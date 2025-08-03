use anyhow::Result;
use clap::Parser;
use crossbeam_channel::{bounded, Receiver, Sender};
use fastq::{each_zipped, OwnedRecord, Parser as FastqParser, Record};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::thread;



#[derive(Parser)]
#[command(name = "fastq_processor")]
#[command(about = "Process R1 and R2 FASTQ files")]
struct Args {
    #[arg(short = '1', long, help = "Input R1 FASTQ file")]
    r1_input: PathBuf,
    
    #[arg(short = '2', long, help = "Input R2 FASTQ file")]
    r2_input: PathBuf,
    
    #[arg(short = 'o', long, help = "Output prefix")]
    output_prefix: String,
    
    #[arg(short = 't', long, default_value = "4", help = "Number of threads")]
    threads: usize,
    
    #[arg(short = 'b', long, default_value = "200000", help = "Batch size for processing")]
    batch_size: usize,
    
    #[arg(short = 'v', long, default_value = "false", help = "Verbose output showing progress")]
    verbose: bool,
    
    #[arg(short = 'c', long, default_value = "false", help = "Compress output files with gzip")]
    compress: bool,
    
    #[arg(short = 'n', long, default_value = "001", help = "Number suffix for output files (e.g., 001, 002)")]
    number_suffix: String,
}

// gzip 或 plain FASTQ 都能自动判断
fn open_fastq<P: AsRef<Path>>(p: P) -> Box<dyn Read + Send> {
    let f = File::open(p.as_ref()).unwrap();
    match p.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gz") => Box::new(MultiGzDecoder::new(f)),
        _          => Box::new(f),
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _    => b'N',
    }).collect()
}

fn extract_base_header(head: &[u8]) -> &[u8] {
    if head.ends_with(b"/1") || head.ends_with(b"/2") { &head[..head.len()-2] } else { head }
}

/// 把两条 FASTQ 读成 batch，发到下游
fn reader_thread(
    r1_path: &Path,
    r2_path: &Path,
    batch_len: usize,
    tx: Sender<(Vec<OwnedRecord>, Vec<OwnedRecord>)>,
) -> Result<()> {
    // 构造两个 parser
    let p1 = FastqParser::new(open_fastq(r1_path));
    let p2 = FastqParser::new(open_fastq(r2_path));

    let mut r1_batch = Vec::with_capacity(batch_len);
    let mut r2_batch = Vec::with_capacity(batch_len);

    // fastq‑rs 原生的"成对遍历"——每回调一次就是一对 read
    each_zipped(p1, p2, |opt1, opt2| {
        match (opt1, opt2) {
            (Some(r1), Some(r2)) => {
                r1_batch.push(r1.to_owned_record()); // OwnedRecord = 结构体版 FASTQ
                r2_batch.push(r2.to_owned_record());
                // 满了就发
                if r1_batch.len() == batch_len {
                    tx.send((r1_batch.split_off(0), r2_batch.split_off(0))).unwrap();
                }
                (true, true) // 两个 parser 都继续
            }
            // 文件长度不一致时提前终止
            _ => (false, false),
        }
    })?;

    if !r1_batch.is_empty() {
        tx.send((r1_batch, r2_batch)).unwrap();
    }
    Ok(())
}

fn create_writer(path: &PathBuf) -> Result<Box<dyn Write + Send>> {
    let file = File::create(path)?;

    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        // ① 更低压缩等级：level 1≈4～5 倍速度
        let encoder = GzEncoder::new(file, Compression::new(1));
        // ② 更大的 BufWriter：1 MiB 而非 8 KiB，减少 sys‑call 次数
        Ok(Box::new(BufWriter::with_capacity(4 << 20, encoder)))
    } else {
        Ok(Box::new(BufWriter::with_capacity(4 << 20, file)))
    }
}

struct ProcessedRecord {
    r1_out: OwnedRecord,
    r2_out: OwnedRecord,
    r3_out: OwnedRecord,
}

fn process_pair(
    r1: OwnedRecord,
    r2: OwnedRecord,
) -> Option<(OwnedRecord, OwnedRecord, OwnedRecord)> {
    if r2.seq().len() != 166 { return None; }

    let id1 = extract_base_header(r1.head());
    let id2 = extract_base_header(r2.head());
    if id1 != id2 { return None; }

    // ---------- R1 ----------
    let id1_vec = id1.to_vec();
    let mut out1 = r1;             // 复用内存；只需截 ID
    out1.head = id1_vec.clone();

    // ---------- R2 ----------
    let (tail_seq, head_seq) = r2.seq().split_at(150); // 0..150, 150..166
    let (tail_qual, head_qual) = r2.qual().split_at(150);

    let out2 = OwnedRecord {
        head : id1_vec.clone(),
        seq  : reverse_complement(head_seq),
        qual : head_qual.iter().rev().cloned().collect(),
        sep  : None,
    };

    // ---------- R3 ----------
    let out3 = OwnedRecord {
        head : id1_vec,
        seq  : tail_seq.to_vec(),
        qual : tail_qual.to_vec(),
        sep  : None,
    };
    Some((out1, out2, out3))
}

fn process_batch(r1_batch: Vec<OwnedRecord>, r2_batch: Vec<OwnedRecord>) -> Vec<ProcessedRecord> {
    let mut results = Vec::new();
    
    for (r1, r2) in r1_batch.into_iter().zip(r2_batch.into_iter()) {
        if let Some((r1_out, r2_out, r3_out)) = process_pair(r1, r2) {
            results.push(ProcessedRecord {
                r1_out,
                r2_out,
                r3_out,
            });
        }
    }
    
    results
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Set up output file paths
    let extension = if args.compress { ".fastq.gz" } else { ".fastq" };
    let r1_output = PathBuf::from(format!("{}_S1_L001_R1_{}{}", args.output_prefix, args.number_suffix, extension));
    let r2_output = PathBuf::from(format!("{}_S1_L001_R2_{}{}", args.output_prefix, args.number_suffix, extension));
    let r3_output = PathBuf::from(format!("{}_S1_L001_R3_{}{}", args.output_prefix, args.number_suffix, extension));
    
    // Clone for printing later
    let r1_output_display = r1_output.clone();
    let r2_output_display = r2_output.clone();
    let r3_output_display = r3_output.clone();
    
    if args.verbose {
        println!("Starting batch processing with batch size: {}", args.batch_size);
    }
    
    // Create channels for batch processing - 增加缓冲区大小
    let (batch_tx, batch_rx): (Sender<(Vec<OwnedRecord>, Vec<OwnedRecord>)>, Receiver<(Vec<OwnedRecord>, Vec<OwnedRecord>)>) = bounded(50);
    let (output_tx, output_rx): (Sender<Vec<ProcessedRecord>>, Receiver<Vec<ProcessedRecord>>) = bounded(50);
    
    // Statistics
    let processed_count = Arc::new(Mutex::new(0usize));
    let filtered_count = Arc::new(Mutex::new(0usize));
    let total_read = Arc::new(Mutex::new(0usize));
    
    // Start reader thread
    let r1_input = args.r1_input.clone();
    let r2_input = args.r2_input.clone();
    let batch_size = args.batch_size;
    let verbose = args.verbose;
    let _read_count = Arc::clone(&total_read);
    let reader_handle = thread::spawn(move || -> Result<()> {
        reader_thread(&r1_input, &r2_input, batch_size, batch_tx)?;
        if verbose {
            println!("Finished reading record pairs");
        }
        Ok(())
    });
    
    // Start processing threads
    let mut processing_handles = Vec::new();
    for _ in 0..args.threads {
        let rx = batch_rx.clone();
        let tx = output_tx.clone();
        let proc_count = Arc::clone(&processed_count);
        let filt_count = Arc::clone(&filtered_count);
        
        let handle = thread::spawn(move || {
            while let Ok((r1_batch, r2_batch)) = rx.recv() {
                let total_in_batch = r2_batch.len();
                let results = process_batch(r1_batch, r2_batch);
                
                let processed_in_batch = results.len();
                let filtered_in_batch = total_in_batch - processed_in_batch;
                
                *proc_count.lock().unwrap() += processed_in_batch;
                *filt_count.lock().unwrap() += filtered_in_batch;
                
                if !results.is_empty() {
                    if tx.send(results).is_err() {
                        break;
                    }
                }
            }
        });
        processing_handles.push(handle);
    }
    
    // Create separate channels for each output file
    let (r1_tx, r1_rx): (Sender<Vec<OwnedRecord>>, Receiver<Vec<OwnedRecord>>) = bounded(50);
    let (r2_tx, r2_rx): (Sender<Vec<OwnedRecord>>, Receiver<Vec<OwnedRecord>>) = bounded(50);
    let (r3_tx, r3_rx): (Sender<Vec<OwnedRecord>>, Receiver<Vec<OwnedRecord>>) = bounded(50);
    
    // Distribution thread - 分发处理结果到各个写入线程
    let verbose_dist = args.verbose;
    let dist_handle = {
        let r1_tx_clone = r1_tx.clone();
        let r2_tx_clone = r2_tx.clone();
        let r3_tx_clone = r3_tx.clone();
        thread::spawn(move || -> Result<()> {
            let mut written_count = 0;
            while let Ok(batch_results) = output_rx.recv() {
                let mut r1_batch = Vec::new();
                let mut r2_batch = Vec::new();
                let mut r3_batch = Vec::new();
                
                for processed in batch_results {
                    r1_batch.push(processed.r1_out);
                    r2_batch.push(processed.r2_out);
                    r3_batch.push(processed.r3_out);
                    written_count += 1;
                }
                
                // 并行发送到各个写入线程
                if !r1_batch.is_empty() {
                    r1_tx_clone.send(r1_batch).map_err(|_| anyhow::anyhow!("Failed to send R1 batch"))?;
                    r2_tx_clone.send(r2_batch).map_err(|_| anyhow::anyhow!("Failed to send R2 batch"))?;
                    r3_tx_clone.send(r3_batch).map_err(|_| anyhow::anyhow!("Failed to send R3 batch"))?;
                }
                
                if verbose_dist && written_count % 100000 == 0 {
                    println!("Written {} records...", written_count);
                }
            }
            if verbose_dist {
                println!("Finished writing {} records", written_count);
            }
            Ok(())
        })
    };
    
    // Start separate writer threads for each output file
    let r1_writer_handle = {
        let r1_output = r1_output.clone();
        thread::spawn(move || -> Result<()> {
            let mut writer = create_writer(&r1_output)?;
            while let Ok(batch) = r1_rx.recv() {
                for record in batch {
                    record.write(&mut writer)?;   // fastq‑rs 一条调用完成
                }
            }
            Ok(())
        })
    };
    
    let r2_writer_handle = {
        let r2_output = r2_output.clone();
        thread::spawn(move || -> Result<()> {
            let mut writer = create_writer(&r2_output)?;
            while let Ok(batch) = r2_rx.recv() {
                for record in batch {
                    record.write(&mut writer)?;   // fastq‑rs 一条调用完成
                }
            }
            Ok(())
        })
    };
    
    let r3_writer_handle = {
        let r3_output = r3_output.clone();
        thread::spawn(move || -> Result<()> {
            let mut writer = create_writer(&r3_output)?;
            while let Ok(batch) = r3_rx.recv() {
                for record in batch {
                    record.write(&mut writer)?;   // fastq‑rs 一条调用完成
                }
            }
            Ok(())
        })
    };
    
    // Wait for reader to finish
    reader_handle.join().unwrap()?;
    
    // Wait for all processing threads to finish
    for handle in processing_handles {
        handle.join().unwrap();
    }
    
    // Close output channel to signal distribution thread to finish
    drop(output_tx);
    
    // Wait for distribution thread to finish
    dist_handle.join().unwrap()?;
    
    // Close writer channels to signal writers to finish
    drop(r1_tx);
    drop(r2_tx);
    drop(r3_tx);
    
    // Wait for all writer threads to finish
    r1_writer_handle.join().unwrap()?;
    r2_writer_handle.join().unwrap()?;
    r3_writer_handle.join().unwrap()?;
    
    let final_processed = *processed_count.lock().unwrap();
    let final_filtered = *filtered_count.lock().unwrap();
    
    println!("Processing complete!");
    println!("Processed records: {}", final_processed);
    println!("Filtered out records: {}", final_filtered);
    println!("Output files:");
    println!("  R1: {}", r1_output_display.display());
    println!("  R2: {}", r2_output_display.display());
    println!("  R3: {}", r3_output_display.display());
    
    Ok(())
}
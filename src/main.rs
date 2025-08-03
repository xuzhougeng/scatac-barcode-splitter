use anyhow::Result;
use clap::Parser;
use crossbeam_channel::{bounded, Receiver, Sender};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
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

#[derive(Debug, Clone)]
struct FastqRecord {
    header: String,
    sequence: String,
    plus: String,
    quality: String,
}

impl FastqRecord {
    fn new(header: String, sequence: String, plus: String, quality: String) -> Self {
        FastqRecord {
            header,
            sequence,
            plus,
            quality,
        }
    }
    
    // 直接写入到buffer的方法
    fn write_to_bytes(&self, buffer: &mut Vec<u8>) {
        buffer.extend_from_slice(self.header.as_bytes());
        buffer.push(b'\n');
        buffer.extend_from_slice(self.sequence.as_bytes());
        buffer.push(b'\n');
        buffer.extend_from_slice(self.plus.as_bytes());
        buffer.push(b'\n');
        buffer.extend_from_slice(self.quality.as_bytes());
        buffer.push(b'\n');
    }
}

fn complement_base(base: char) -> char {
    match base.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        'N' => 'N',
        _ => base,
    }
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(complement_base)
        .collect()
}

fn extract_base_header(header: &str) -> String {
    if header.ends_with("/1") || header.ends_with("/2") {
        header[..header.len() - 2].to_string()
    } else {
        header.to_string()
    }
}

fn read_fastq_record<R: BufRead>(lines: &mut std::io::Lines<R>) -> Result<Option<FastqRecord>> {
    // Read header line
    let header = loop {
        if let Some(line) = lines.next() {
            let line = line?;
            if line.starts_with('@') {
                break line;
            }
            // Skip non-header lines (empty lines, etc.)
        } else {
            return Ok(None); // End of file
        }
    };
    
    // Read sequence lines until we hit a '+' line
    let mut sequence = String::new();
    loop {
        if let Some(line) = lines.next() {
            let line = line?;
            if line.starts_with('+') {
                // This is the plus line, stop reading sequence
                break;
            } else {
                sequence.push_str(&line);
            }
        } else {
            return Err(anyhow::anyhow!("Unexpected end of file while reading sequence"));
        }
    }
    
    // Read quality lines until we have the same length as sequence
    let mut quality = String::new();
    while quality.len() < sequence.len() {
        if let Some(line) = lines.next() {
            let line = line?;
            quality.push_str(&line);
        } else {
            return Err(anyhow::anyhow!("Unexpected end of file while reading quality"));
        }
    }
    
    // Trim quality to exact sequence length (in case we read too much)
    quality.truncate(sequence.len());
    
    Ok(Some(FastqRecord::new(header, sequence, "+".to_string(), quality)))
}

fn read_fastq_batch<R: BufRead>(lines: &mut std::io::Lines<R>, batch_size: usize) -> Result<Vec<FastqRecord>> {
    let mut batch = Vec::with_capacity(batch_size);
    
    for _ in 0..batch_size {
        if let Some(record) = read_fastq_record(lines)? {
            batch.push(record);
        } else {
            break; // End of file
        }
    }
    
    Ok(batch)
}

fn open_reader(path: &PathBuf) -> Result<Box<dyn BufRead + Send>> {
    let file = File::open(path)?;
    
    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        let decoder = MultiGzDecoder::new(file);
        // 增加缓冲区到2MB
        Ok(Box::new(BufReader::with_capacity(2 << 20, decoder)))
    } else {
        // 增加缓冲区到2MB
        Ok(Box::new(BufReader::with_capacity(2 << 20, file)))
    }
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
    r1_out: FastqRecord,
    r2_out: FastqRecord,
    r3_out: FastqRecord,
}

fn process_record_pair(r1: FastqRecord, r2: FastqRecord) -> Option<ProcessedRecord> {
    // Check R2 length is 166bp
    if r2.sequence.len() != 166 {
        return None;
    }
    
    // Check headers match (after removing /1 and /2)
    let r1_base = extract_base_header(&r1.header);
    let r2_base = extract_base_header(&r2.header);
    
    if r1_base != r2_base {
        return None;
    }
    
    // Process R1: remove /1 from header
    let r1_out = FastqRecord::new(
        r1_base.clone(),
        r1.sequence,
        r1.plus,
        r1.quality,
    );
    
    // Process R2: positions 151-166 (16bp), reverse complement
    let r2_seq = &r2.sequence[150..166]; // 0-based indexing, so 150..166 for 151-166
    let r2_qual = &r2.quality[150..166];
    let r2_rc_seq = reverse_complement(r2_seq);
    let r2_rc_qual: String = r2_qual.chars().rev().collect(); // reverse quality scores too
    
    let r2_out = FastqRecord::new(
        r1_base.clone(),
        r2_rc_seq,
        r2.plus.clone(),
        r2_rc_qual,
    );
    
    // Process R3: positions 1-150 (150bp), forward
    let r3_seq = &r2.sequence[0..150];
    let r3_qual = &r2.quality[0..150];
    
    let r3_out = FastqRecord::new(
        r1_base,
        r3_seq.to_string(),
        r2.plus,
        r3_qual.to_string(),
    );
    
    Some(ProcessedRecord {
        r1_out,
        r2_out,
        r3_out,
    })
}

fn process_batch(r1_batch: Vec<FastqRecord>, r2_batch: Vec<FastqRecord>) -> Vec<ProcessedRecord> {
    let mut results = Vec::new();
    
    for (r1, r2) in r1_batch.into_iter().zip(r2_batch.into_iter()) {
        if let Some(processed) = process_record_pair(r1, r2) {
            results.push(processed);
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
    let (batch_tx, batch_rx): (Sender<(Vec<FastqRecord>, Vec<FastqRecord>)>, Receiver<(Vec<FastqRecord>, Vec<FastqRecord>)>) = bounded(50);
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
    let read_count = Arc::clone(&total_read);
    let reader_handle = thread::spawn(move || -> Result<()> {
        let r1_reader = open_reader(&r1_input)?;
        let r2_reader = open_reader(&r2_input)?;
        
        let mut r1_lines = r1_reader.lines();
        let mut r2_lines = r2_reader.lines();
        
        loop {
            let r1_batch = read_fastq_batch(&mut r1_lines, batch_size);
            let r2_batch = read_fastq_batch(&mut r2_lines, batch_size);
            
            match (r1_batch, r2_batch) {
                (Ok(r1_batch), Ok(r2_batch)) => {
                    if r1_batch.is_empty() || r2_batch.is_empty() {
                        if verbose {
                            println!("Reached end of file. R1 batch: {}, R2 batch: {}", r1_batch.len(), r2_batch.len());
                        }
                        break;
                    }
                    
                    let batch_count = r1_batch.len().min(r2_batch.len());
                    *read_count.lock().unwrap() += batch_count;
                    
                    if verbose && *read_count.lock().unwrap() % 1000000 == 0 {
                        println!("Read {} record pairs...", *read_count.lock().unwrap());
                    }
                    
                    if batch_tx.send((r1_batch, r2_batch)).is_err() {
                        println!("Channel send failed, stopping reader");
                        break;
                    }
                }
                (Err(e1), _) => {
                    println!("Error reading R1 batch: {}", e1);
                    return Err(e1);
                }
                (_, Err(e2)) => {
                    println!("Error reading R2 batch: {}", e2);
                    return Err(e2);
                }
            }
        }
        if verbose {
            println!("Finished reading {} record pairs", *read_count.lock().unwrap());
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
    let (r1_tx, r1_rx): (Sender<Vec<FastqRecord>>, Receiver<Vec<FastqRecord>>) = bounded(50);
    let (r2_tx, r2_rx): (Sender<Vec<FastqRecord>>, Receiver<Vec<FastqRecord>>) = bounded(50);
    let (r3_tx, r3_rx): (Sender<Vec<FastqRecord>>, Receiver<Vec<FastqRecord>>) = bounded(50);
    
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
            let mut buffer = Vec::with_capacity(1 << 20); // 1MB buffer
            while let Ok(batch) = r1_rx.recv() {
                buffer.clear();
                for record in batch {
                    record.write_to_bytes(&mut buffer);
                }
                writer.write_all(&buffer)?;
            }
            Ok(())
        })
    };
    
    let r2_writer_handle = {
        let r2_output = r2_output.clone();
        thread::spawn(move || -> Result<()> {
            let mut writer = create_writer(&r2_output)?;
            let mut buffer = Vec::with_capacity(1 << 20); // 1MB buffer
            while let Ok(batch) = r2_rx.recv() {
                buffer.clear();
                for record in batch {
                    record.write_to_bytes(&mut buffer);
                }
                writer.write_all(&buffer)?;
            }
            Ok(())
        })
    };
    
    let r3_writer_handle = {
        let r3_output = r3_output.clone();
        thread::spawn(move || -> Result<()> {
            let mut writer = create_writer(&r3_output)?;
            let mut buffer = Vec::with_capacity(1 << 20); // 1MB buffer
            while let Ok(batch) = r3_rx.recv() {
                buffer.clear();
                for record in batch {
                    record.write_to_bytes(&mut buffer);
                }
                writer.write_all(&buffer)?;
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
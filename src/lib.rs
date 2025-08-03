// lib.rs - 库函数

/// DNA 序列反向互补函数
/// 
/// 将输入的 DNA 序列进行反向互补转换：
/// - A ↔ T
/// - G ↔ C  
/// - 其他字符转为 N
/// - 自动转大写并反向序列
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _    => b'N',
    }).collect()
}

/// 提取 FASTQ header 的基础 ID（移除 /1 或 /2 后缀）
pub fn extract_base_header(head: &[u8]) -> &[u8] {
    if head.ends_with(b"/1") || head.ends_with(b"/2") { &head[..head.len()-2] } else { head }
}
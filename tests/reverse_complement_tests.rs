use scatac_barcode_splitter::reverse_complement;

#[test]
fn test_reverse_complement_basic() {
    // 基本 ATGC 转换
    let input = b"ATGC";
    let expected = b"GCAT";
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}

#[test]
fn test_reverse_complement_lowercase() {
    // 测试小写字母（会被转为大写）
    let input = b"atgc";
    let expected = b"GCAT";
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}

#[test]
fn test_reverse_complement_mixed_case() {
    // 测试大小写混合
    let input = b"AtGc";
    let expected = b"GCAT";
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}

#[test]
fn test_reverse_complement_with_n() {
    // 测试包含 N（未知碱基）
    let input = b"ATGCN";
    let expected = b"NGCAT";
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}

#[test]
fn test_reverse_complement_unknown_bases() {
    // 测试未知碱基（非 ATGC）会被转换为 N
    let input = b"ATXGC";
    let expected = b"GCNAT";  // X 被转换为 N
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}

#[test]
fn test_reverse_complement_empty() {
    // 测试空序列
    let input = b"";
    let expected = b"";
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}

#[test]
fn test_reverse_complement_long_sequence() {
    // 测试较长序列
    let input = b"AAATTTGGGCCC";
    let expected = b"GGGCCCAAATTT";
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}

#[test]
fn test_reverse_complement_palindrome() {
    // 测试回文序列
    let input = b"GAATTC";  // EcoRI 切点
    let expected = b"GAATTC";
    let result = reverse_complement(input);
    assert_eq!(result, expected);
}
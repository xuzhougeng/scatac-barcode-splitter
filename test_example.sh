#!/bin/bash

# 创建测试数据目录
mkdir -p test_data

# 创建示例R1 FASTQ文件
cat > test_data/test_R1.fastq << 'EOF'
@read1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/1
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF

# 创建示例R2 FASTQ文件（166bp）
# 生成166bp的序列：前150bp是ACGT重复，后16bp是TGCA重复
R2_SEQ1=$(python3 -c "print('ACGT' * 37 + 'AC' + 'TGCA' * 4)")
R2_QUAL1=$(python3 -c "print('I' * 166)")
R2_SEQ2=$(python3 -c "print('TGCA' * 37 + 'TG' + 'ACGT' * 4)")
R2_QUAL2=$(python3 -c "print('J' * 166)")

cat > test_data/test_R2.fastq << EOF
@read1/2
$R2_SEQ1
+
$R2_QUAL1
@read2/2
$R2_SEQ2
+
$R2_QUAL2
EOF

# 压缩测试文件
gzip test_data/test_R1.fastq
gzip test_data/test_R2.fastq

echo "测试数据已创建："
ls -la test_data/

echo -e "\n运行测试："
./target/release/scatac-barcode-splitter  \
    -1 test_data/test_R1.fastq.gz \
    -2 test_data/test_R2.fastq.gz \
    -o test_output \
    -t 2

echo -e "\n输出文件："
ls -la test_output_*.fastq.gz 2>/dev/null || echo "未找到输出文件"

echo -e "\n查看输出内容："
if [ -f "test_output_S1_L001_R1_001.fastq.gz" ]; then
    echo "R1输出："
    zcat test_output_S1_L001_R1_001.fastq.gz | head -4
    
    echo -e "\nR2输出（反向互补的151-166bp）："
    zcat test_output_S1_L001_R2_001.fastq.gz | head -4
    
    echo -e "\nR3输出（1-150bp）："
    zcat test_output_S1_L001_R3_001.fastq.gz | head -4
fi
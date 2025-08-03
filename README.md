# 拆分MGI平台的scATAC-seq测序数据

因为华大的拆分方式没办法保留R3， 因此barcode数据并在R2后面，按照如下流程进行数据处理

- 读取R1和R2两个FASTQ文件（支持.gz压缩格式）
- 只处理R2长度为166bp的记录
- 确保R1和R2的header匹配（去掉/1和/2后）
- 输出3个文件：
  - R1: 原始R1序列，删除header中的/1
  - R2: R2序列的151-166bp（16bp），反向互补
  - R3: R2序列的1-150bp（150bp），正向

## 架构设计

程序采用生产者-消费者模式的流式处理架构：

1. **读取线程**: 同时读取R1和R2文件，逐条发送到处理队列
2. **处理线程池**: 多个线程并行处理记录对
3. **写入线程**: 统一写入处理结果到输出文件

这种设计保证了恒定的内存使用量，无论文件多大都不会占用过多内存。

## 编译

```bash
cargo build --release
```

## 使用方法

```bash
./target/release/scatac-barcode-splitter  \
    -1 input_R1.fastq.gz \
    -2 input_R2.fastq.gz \
    -o output_prefix \
    -t 8
```

### 参数说明

- `-1, --r1-input`: 输入R1 FASTQ文件路径
- `-2, --r2-input`: 输入R2 FASTQ文件路径
- `-o, --output-prefix`: 输出文件前缀
- `-t, --threads`: 线程数（默认4）
- `-b, --batch-size`: 批处理大小（默认100000）
- `-n, --number-suffix`: 默认001

### 输出文件

程序会生成以下3个文件：
- `{prefix}_S1_L001_R1_001.fastq.gz`
- `{prefix}_S1_L001_R2_001.fastq.gz`
- `{prefix}_S1_L001_R3_001.fastq.gz`

## 示例

### 基本用法
```bash
./target/release/scatac-barcode-splitter  \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -o processed \
    -t 8 \
    -b 100000
```

### 处理大文件（推荐设置）
```bash
# 大文件优化配置
./target/release/scatac-barcode-splitter  \
    -1 large_R1.fastq.gz \
    -2 large_R2.fastq.gz \
    -o output \
    -t 8 \
    -b 100000
```

### 监控内存使用
```bash
# 后台运行处理程序
./target/release/scatac-barcode-splitter  \
    -1 large_R1.fastq.gz \
    -2 large_R2.fastq.gz \
    -o output \
    -t 8 \
    -b 100000 &

# 监控内存使用
./monitor_memory.sh $(pgrep scatac-barcode-splitter )
```

输出文件：
- `processed_S1_L001_R1_001.fastq.gz`
- `processed_S1_L001_R2_001.fastq.gz`
- `processed_S1_L001_R3_001.fastq.gz`

## 性能特点

- **内存使用恒定**: 无论文件多大，内存使用量都保持在较低水平
- **并行处理**: 读取、处理、写入同时进行，最大化吞吐量
- **实时进度**: 每处理10000条记录显示一次进度
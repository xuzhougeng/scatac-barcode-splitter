#!/bin/bash

# 监控程序内存使用的脚本
# 用法: ./monitor_memory.sh <pid>

if [ $# -eq 0 ]; then
    echo "用法: $0 <进程ID>"
    echo "例如: ./monitor_memory.sh \$(pgrep scatac-barcode-splitter )"
    exit 1
fi

PID=$1

echo "监控进程 $PID 的内存使用..."
echo "时间         RSS(MB)  VSZ(MB)  %MEM   命令"
echo "--------------------------------------------"

while kill -0 $PID 2>/dev/null; do
    ps -p $PID -o pid,rss,vsz,pmem,comm --no-headers | awk '{
        printf "%-12s %-7.1f  %-7.1f  %-5s  %s\n", 
               strftime("%H:%M:%S"), $2/1024, $3/1024, $4, $5
    }'
    sleep 2
done

echo "进程已结束"
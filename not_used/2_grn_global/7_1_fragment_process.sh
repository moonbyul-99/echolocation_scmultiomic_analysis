#!/bin/bash

'''
This script compress the fragment.tsv file to fragment.tsv.gz 
'''

# 定义要处理的父目录
# 请将此路径替换为你的实际目录
PARENT_DIR="/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/fragments"

# 检查父目录是否存在
if [ ! -d "$PARENT_DIR" ]; then
    echo "错误：指定的父目录 '$PARENT_DIR' 不存在。"
    exit 1
fi

echo "开始并行处理 'fragment.tsv' 文件..."
echo "-------------------------------------"

# 使用 find 命令查找所有 fragment.tsv 文件，并对其父目录进行处理
find "$PARENT_DIR" -type f -name "fragment.tsv" | while read -r fragment_file; do
    (
        # 获取 fragment.tsv 所在的目录
        current_dir=$(dirname "$fragment_file")
        file_name=$(basename "$fragment_file")
        output_gz="${current_dir}/${file_name}.gz"
        output_tbi="${current_dir}/${file_name}.gz.tbi"

        echo "正在处理目录: ${current_dir}"

        # 1. 使用 bgzip 压缩 fragment.tsv 文件
        if bgzip -c "$fragment_file" > "$output_gz"; then
            echo "  - ${file_name} 已成功压缩为 ${file_name}.gz"

            # 2. 使用 tabix 创建 .tbi 索引文件
            # 假设 fragment.tsv 的格式是：第一列是序列名，第二列是起始位置，第三列是结束位置
            if tabix -s 1 -b 2 -e 3 "$output_gz"; then
                echo "  - ${file_name}.gz 已成功创建索引 ${file_name}.gz.tbi"
            else
                echo "  - 错误：无法为 ${file_name}.gz 创建索引。请检查 tabix 命令和文件格式。"
            fi
        else
            echo "  - 错误：无法压缩 ${file_name}。请检查 bgzip 命令和文件权限。"
        fi
    ) & # & 符号表示在后台运行子进程，实现并行
done

# 等待所有后台进程完成
wait

echo "-------------------------------------"
echo "所有 'fragment.tsv' 文件的处理已完成。"
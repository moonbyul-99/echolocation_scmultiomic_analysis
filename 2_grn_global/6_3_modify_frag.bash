#!/bin/bash

# Description: 
# This bash script is used to modify the fragment file name

# 定义主目录
MAIN_DIR="/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/fragments" # 请将此替换为您的实际目录路径

# 遍历所有 {sample}_fragments 子目录
for dir in "$MAIN_DIR"/*_fragments/; do
    if [ -d "$dir" ]; then
        SAMPLE_NAME=$(basename "$dir" | sed 's/_fragments//') # 提取 sample 名称

        FRAGMENT_FILE="${dir}fragment.tsv"

        if [ -f "$FRAGMENT_FILE" ]; then
            echo "正在处理文件: $FRAGMENT_FILE"
            # 使用 awk 修改第四列并保存到临时文件，然后替换原文件
            awk -v sample="$SAMPLE_NAME" 'BEGIN{OFS="\t"} {$4=sample"-"$4; print}' "$FRAGMENT_FILE" > "${FRAGMENT_FILE}.tmp" && mv "${FRAGMENT_FILE}.tmp" "$FRAGMENT_FILE"
        else
            echo "警告: 未找到文件 $FRAGMENT_FILE"
        fi
    fi
done

echo "所有文件处理完毕。"
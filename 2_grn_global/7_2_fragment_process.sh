#!/bin/bash
'''
This script is used to rename the cell barcodes in fragment.tsv files to {cell_id}_{sample} format instead of {sample}-{cell_id} format.
'''

# 定义主目录
MAIN_DIR="/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/fragments" #请将此替换为您的实际目录路径

echo "开始处理目录：$MAIN_DIR"

# 遍历所有 {sample}_fragments 子目录
for dir in "$MAIN_DIR"/*_fragments/; do
    if [ -d "$dir" ]; then
        # 提取 sample 名称 (例如 "sample1_fragments" -> "sample1")
        SAMPLE_NAME=$(basename "$dir" | sed 's/_fragments//')

        FRAGMENT_FILE="${dir}fragment.tsv"

        if [ -f "$FRAGMENT_FILE" ]; then
            echo "正在处理文件: $FRAGMENT_FILE"
            # 使用 awk 修改第四列。
            # 首先，将第四列的 sample 部分去除，得到 cell_id。
            # 然后，将 cell_id 和 sample_name 重新组合。
            # AWK_SCRIPT='BEGIN{OFS="\t"} {$4=gensub(/^[^_]*-/, "", "1", $4); $4=$4"_"sample; print}' # 尝试过这种方式，但对于cell_id本身带横线会更复杂
            # 更可靠的做法是先替换掉开头的 {sample}-，然后将剩余部分（即 cell_id）与 sample_name 组合。
            awk -v sample="$SAMPLE_NAME" 'BEGIN{OFS="\t"} {
                # 检查第四列是否以 sample_name- 开头
                if (substr($4, 1, length(sample)+1) == sample"-") {
                    cell_id = substr($4, length(sample)+2); # 提取 cell_id
                    $4 = cell_id "_" sample; # 重新组合为 {cell_id}_{sample}
                } else {
                    # 如果不匹配预期格式，打印警告或按原样保留
                    print "警告: " FILENAME "中第4列的barcode '" $4 "'不符合预期 {sample}-{cell_id} 格式，跳过修改。";
                    # 可以在这里选择保留原样或进行其他处理
                }
                print $0; # 打印整行
            }' "$FRAGMENT_FILE" > "${FRAGMENT_FILE}.tmp" && mv "${FRAGMENT_FILE}.tmp" "$FRAGMENT_FILE"
        else
            echo "警告: 未在目录 $dir 中找到 fragment.tsv 文件。"
        fi
    fi
done

echo "所有文件处理完毕。"
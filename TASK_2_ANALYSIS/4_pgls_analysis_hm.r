library(ape)
library(geomorph)

# ==========================================
# 1. 环境准备与全局数据加载
# ==========================================
# 定义数据所在文件夹
data_dir <- "./pgls_data/"
output_file <- "./pgls_batch_hm_results.csv"

# 载入5物种系统发育树 (请根据您的Newick文本替换或读取文件)
tree_text <- "(((Rhinolophus_sinicus:58.00000000,Cynopterus_sphinx:58.00000000)'14':4.00000000,Miniopterus_fuliginosus:62.00000000)'13':32.00000000,(typhlomys_daloushanensis:48.38724000,Mus_musculus:48.38724000)'25':45.61276000);"
phylo_tree <- read.tree(text = tree_text)

# 将树的物种名规范化为您的简写 (需与CSV中的行名完全一致)
# 映射关系: Rhinolophus_sinicus->JT, Miniopterus_fuliginosus->CY, Cynopterus_sphinx->QF, typhlomys_daloushanensis->T, Mus_musculus->M
phylo_tree$tip.label <- c("JT", "QF", "CY", "T", "M") 

# 载入回声定位表型数据
pheno_data <- read.csv(paste0(data_dir, "phenotype.csv"), row.names = 1)

# 定义要遍历的 8 个脑区名称 (需与您的文件名后缀一致)
regions <- c('OPC', 'astrocyte', 'ependymal', 'excitatory neuron', 'inhibitory neuron', 'microglia', 'oligodendrocyte')

# 初始化一个空的数据框，用于结构化存储最终结果
results_df <- data.frame(
  Region = character(),
  R2 = numeric(),        # 回声定位对转录组差异的解释力度
  F_statistic = numeric(),# F 统计量
  p_value = numeric(),    # 显著性 p 值 (基于置换检验)
  Z_score = numeric(),    # 效应量大小 (Z 分数)
  stringsAsFactors = FALSE
)

# ==========================================
# 2. FOR 循环批量进行 PGLS 检验
# ==========================================
for (region in regions) {
    file_path <- paste0(data_dir, "expr_cell_", region, ".csv")
    
    # 检查文件是否存在，防止个别脑区数据缺失导致中断
    if (!file.exists(file_path)) {
        warning(paste("文件未找到，跳过脑区:", region))
        next
    }
    
    # 读取该脑区的表达矩阵 (第一列物种名作为行名)
    Y_matrix <- as.matrix(read.csv(file_path, row.names = 1))
    
    # 关键步骤：确保转录组、表型、进化树的物种顺序完美严格对齐
    species_order <- phylo_tree$tip.label
    Y_matrix <- Y_matrix[species_order, , drop = FALSE]
    X_vector <- pheno_data[species_order, "Echolocation", drop = FALSE]
    
    # 构建用于 procD.pgls 的高维数据框环境
    gdf <- geomorph.data.frame(
        Y = Y_matrix,
        X = X_vector,
        phy = phylo_tree
    )
    
    # 执行系统发育广义最小二乘法的高维置换检验 (RRPP方式)
    # iter = 999 表示进行 999 次随机置换检验
    pgls_model <- procD.pgls(Y ~ X.Echolocation, 
                             phy = phylo_tree, 
                             data = gdf, 
                            #  lambda = 0.75,
                            #  SS.type = "II",
                             iter = 999, print.progress = FALSE)
    
    # 从 ANOVA 结果表中精确提取结构化统计指标
    anova_table <- pgls_model$aov.table
    
    # 1. 精确提取数值（针对你的 anova_table 结构）
    r2      <- anova_table$Rsq[1]         # 对应表中的 'Rsq'
    f_val   <- anova_table$F[1]           # 对应表中的 'F'
    p_val   <- anova_table$`Pr(>F)`[1]    # 对应表中的 'Pr(>F)'
    z_score <- anova_table$Z[1]           # 对应表中的 'Z'

    # 2. 将当前脑区的结果追加到总表中
    current_res <- data.frame(
    Region = region,
    R2 = r2,
    F_statistic = f_val,
    p_value = p_val,
    Z_score = z_score,
    stringsAsFactors = FALSE
    )
    results_df <- rbind(results_df, current_res)
}

# ==========================================
# 3. 结构化保存与输出
# ==========================================
# 按相关性效应量 (Z_score) 从大到小进行系统排序
results_df <- results_df[order(-results_df$Z_score), ]

# 保存为 CSV 文件，方便您直接阅读或用 Python 的 pandas 载入
write.csv(results_df, file = output_file, row.names = FALSE)

# 在 R 控制台打印预览
print("--- 脑区系统发育相关性检验完成 (按效应量排序) ---")
print(results_df)
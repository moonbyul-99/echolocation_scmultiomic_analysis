library(ape)
library(geomorph)
library(ggplot2)
library(reshape2)
library(viridis)

# ==========================================
# 1. 环境准备与全局数据加载
# ==========================================
data_dir    <- "./pgls_data/"
output_dir  <- "./pgls_data/param_sweep/"

# 参数网格定义
pc_seq    <- c( 30, 40, 50)
lambda_seq <- c(0.25, 0.5, 0.75, 1.0)

# 载入系统发育树
tree_text <- "(((Rhinolophus_sinicus:58.00000000,Cynopterus_sphinx:58.00000000)'14':4.00000000,Miniopterus_fuliginosus:62.00000000)'13':32.00000000,(typhlomys_daloushanensis:48.38724000,Mus_musculus:48.38724000)'25':45.61276000);"
phylo_tree <- read.tree(text = tree_text)
phylo_tree$tip.label <- c("JT", "QF", "CY", "T", "M")

# 载入表型
pheno_data <- read.csv(paste0(data_dir, "phenotype.csv"), row.names = 1)

# 脑区列表
regions <- c('OPC', 'astrocyte', 'ependymal', 'excitatory neuron',
             'inhibitory neuron', 'microglia', 'oligodendrocyte')

# ==========================================
# 2. 参数网格扫描 PGLS
# ==========================================
results_list <- list()

for (region in regions) {
  file_path <- paste0(data_dir, "expr_cell_", region, ".csv")
  if (!file.exists(file_path)) {
    warning(paste("文件未找到，跳过:", region))
    next
  }

  # 读取完整表达矩阵 (50 个 PCs)
  Y_full <- as.matrix(read.csv(file_path, row.names = 1))
  species_order <- phylo_tree$tip.label
  Y_full <- Y_full[species_order, , drop = FALSE]
  X_vector <- pheno_data[species_order, "Echolocation", drop = FALSE]

  cat(sprintf("\n--- 脑区: %s ---\n", region))

  for (n_pc in pc_seq) {
    # 截取前 n_pc 个 PC
    Y_sub <- Y_full[, 1:n_pc, drop = FALSE]

    for (lam in lambda_seq) {
      cat(sprintf("  PCs=%d, lambda=%.2f ... ", n_pc, lam))

      # 构建 geomorph 数据框
      gdf <- geomorph.data.frame(Y = Y_sub, X = X_vector, phy = phylo_tree)

      # 执行 PGLS
      pgls_model <- tryCatch({
        procD.pgls(Y ~ X.Echolocation,
                   phy = phylo_tree,
                   data = gdf,
                   lambda = lam,
                   iter = 999,
                   print.progress = FALSE)
      }, error = function(e) {
        warning(sprintf("PGLS failed for region=%s PCs=%d lambda=%.2f: %s",
                        region, n_pc, lam, e$message))
        return(NULL)
      })

      if (is.null(pgls_model)) {
        cat("FAILED\n")
        # 填入 NA
        results_list[[length(results_list) + 1]] <- data.frame(
          Region    = region,
          PCs       = n_pc,
          Lambda    = lam,
          R2        = NA,
          F_stat    = NA,
          p_value   = NA,
          Z_score   = NA,
          stringsAsFactors = FALSE
        )
        next
      }

      anova_table <- pgls_model$aov.table
      r2      <- anova_table$Rsq[1]
      f_val   <- anova_table$F[1]
      p_val   <- anova_table$`Pr(>F)`[1]
      z_score <- anova_table$Z[1]

      cat(sprintf("Z=%.3f, p=%.2e\n", z_score, p_val))

      results_list[[length(results_list) + 1]] <- data.frame(
        Region    = region,
        PCs       = n_pc,
        Lambda    = lam,
        R2        = r2,
        F_stat    = f_val,
        p_value   = p_val,
        Z_score   = z_score,
        stringsAsFactors = FALSE
      )
    }
  }
}

# 合并结果
results_df <- do.call(rbind, results_list)

# ==========================================
# 3. 保存结构化结果
# ==========================================
write.csv(results_df,
          file = paste0(output_dir, "pgls_param_results.csv"),
          row.names = FALSE)

# ==========================================
# 4. 稳定性汇总：各脑区的统计量均值与变异
# ==========================================
stability <- aggregate(cbind(R2, Z_score) ~ Region, data = results_df,
                       FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                           sd   = sd(x, na.rm = TRUE),
                                           cv   = sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))))
# 展平 aggregate 的多列结果
stability_flat <- data.frame(
  Region         = stability$Region,
  R2_mean        = stability$R2[, "mean"],
  R2_sd          = stability$R2[, "sd"],
  R2_cv          = stability$R2[, "cv"],
  Z_mean         = stability$Z_score[, "mean"],
  Z_sd           = stability$Z_score[, "sd"],
  Z_cv           = stability$Z_score[, "cv"],
  stringsAsFactors = FALSE
)
write.csv(stability_flat,
          file = paste0(output_dir, "stability_summary.csv"),
          row.names = FALSE)

# ==========================================
# 5. 可视化
# ==========================================

# 公共绘图辅助
results_df$PCs    <- as.factor(results_df$PCs)
results_df$Lambda <- as.factor(results_df$Lambda)
results_df$neg_log10_p <- -log10(results_df$p_value)
results_df$neg_log10_p[is.infinite(results_df$neg_log10_p)] <- NA

# -------- 5a. Z-score 热力图 (每个脑区一张子图) --------
pdf(paste0(output_dir, "heatmap_zscore.pdf"), width = 16, height = 12)
p <- ggplot(results_df, aes(x = Lambda, y = PCs, fill = Z_score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Z_score)), size = 3.2) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Z-score") +
  facet_wrap(~ Region, ncol = 3) +
  labs(title = "PGLS Z-score across PCs × Lambda",
       subtitle = "Echolocation trait association (higher |Z| = stronger effect)",
       x = "Lambda (phylogenetic signal)", y = "Number of PCs") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
print(p)
dev.off()
cat("Saved: heatmap_zscore.pdf\n")

# -------- 5b. R² 热力图 --------
pdf(paste0(output_dir, "heatmap_r2.pdf"), width = 16, height = 12)
p <- ggplot(results_df, aes(x = Lambda, y = PCs, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", R2)), size = 3.2) +
  scale_fill_viridis(option = "magma", direction = -1, name = expression(R^2)) +
  facet_wrap(~ Region, ncol = 3) +
  labs(title = "PGLS R² across PCs × Lambda",
       subtitle = "Variance explained by echolocation trait",
       x = "Lambda (phylogenetic signal)", y = "Number of PCs") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
print(p)
dev.off()
cat("Saved: heatmap_r2.pdf\n")

# -------- 5c. -log10(p) 热力图 --------
pdf(paste0(output_dir, "heatmap_pvalue.pdf"), width = 16, height = 12)
p <- ggplot(results_df, aes(x = Lambda, y = PCs, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(!is.na(neg_log10_p),
                               sprintf("%.1f", neg_log10_p), "NA")),
            size = 3.2) +
  scale_fill_viridis(option = "plasma", direction = -1,
                     name = expression(-log[10](p))) +
  facet_wrap(~ Region, ncol = 3) +
  labs(title = "PGLS -log10(p) across PCs × Lambda",
       subtitle = "Higher values = stronger significance",
       x = "Lambda (phylogenetic signal)", y = "Number of PCs") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
print(p)
dev.off()
cat("Saved: heatmap_pvalue.pdf\n")

# -------- 5d. 汇总图：各脑区 Z-score 随 PC 和 Lambda 的变化趋势 --------
pdf(paste0(output_dir, "trend_zscore.pdf"), width = 14, height = 8)
p <- ggplot(results_df, aes(x = PCs, y = Z_score,
                             color = Lambda, group = Lambda)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ Region, ncol = 4, scales = "free_y") +
  scale_color_viridis_d(option = "D", name = "Lambda") +
  labs(title = "Z-score trends: PC dimensions vs Lambda",
       x = "Number of PCs", y = "Z-score") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))
print(p)
dev.off()
cat("Saved: trend_zscore.pdf\n")

# -------- 5e. 排名稳定性热图：各细胞类型在不同参数下的回声相关排名 --------
# 在每个 (PCs, Lambda) 组合内，按 Z-score 降序排名
# rank=1 表示该脑区在该参数设置下回声关联最强
rank_df <- results_df
rank_df$rank <- NA

# 使用原始数值分组赋排名（避免 factor 排序问题）
pc_vals    <- sort(unique(as.numeric(as.character(results_df$PCs))))
lambda_vals <- sort(unique(as.numeric(as.character(results_df$Lambda))))

for (pc in pc_vals) {
  for (lam in lambda_vals) {
    idx <- which(rank_df$PCs == as.character(pc) & rank_df$Lambda == as.character(lam))
    if (length(idx) == 0) next
    z_vals <- rank_df$Z_score[idx]
    # 降序排名：最高 Z-score 获 rank 1
    rank_df$rank[idx] <- rank(-z_vals, ties.method = "min", na.last = "keep")
  }
}
rank_df$rank <- as.integer(rank_df$rank)

# 保存排名数据
write.csv(rank_df[, c("Region", "PCs", "Lambda", "Z_score", "rank")],
          file = paste0(output_dir, "pgls_rank_results.csv"),
          row.names = FALSE)

# 排名热力图：每个脑区一个子图，颜色=排名（1=最红，N=最蓝）
n_regions <- length(unique(rank_df$Region))
pdf(paste0(output_dir, "heatmap_rank.pdf"), width = 16, height = 12)
p <- ggplot(rank_df, aes(x = Lambda, y = PCs, fill = rank)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("#%d", rank)), size = 3.2, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#d73027", "#fc8d59", "#fee090", "#e0f3f8", "#4575b4"),
    limits = c(1, n_regions),
    breaks = 1:n_regions,
    name = "Rank",
    guide = guide_colorbar(reverse = FALSE)
  ) +
  facet_wrap(~ Region, ncol = 3) +
  labs(title = "Region rank stability of echolocation association across PCs × Lambda",
       subtitle = sprintf("Rank #1 = strongest Z-score; Rank #%d = weakest (within each parameter setting)", n_regions),
       x = "Lambda (phylogenetic signal)", y = "Number of PCs") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
print(p)
dev.off()
cat("Saved: heatmap_rank.pdf\n")

# 排名波动总结：每个脑区排名的均值与标准差
rank_summary <- aggregate(rank ~ Region, data = rank_df,
                          FUN = function(x) c(mean_rank = mean(x, na.rm = TRUE),
                                              sd_rank   = sd(x, na.rm = TRUE),
                                              min_rank  = min(x, na.rm = TRUE),
                                              max_rank  = max(x, na.rm = TRUE)))
rank_summary_flat <- data.frame(
  Region    = rank_summary$Region,
  mean_rank = round(rank_summary$rank[, "mean_rank"], 2),
  sd_rank   = round(rank_summary$rank[, "sd_rank"], 2),
  min_rank  = rank_summary$rank[, "min_rank"],
  max_rank  = rank_summary$rank[, "max_rank"],
  stringsAsFactors = FALSE
)
rank_summary_flat <- rank_summary_flat[order(rank_summary_flat$mean_rank), ]
write.csv(rank_summary_flat,
          file = paste0(output_dir, "rank_summary.csv"),
          row.names = FALSE)

cat(sprintf("\n--- Rank Stability Summary ---\n"))
print(rank_summary_flat)

# -------- 5f. 稳定性气泡图：Z_mean vs CV --------
pdf(paste0(output_dir, "stability_bubble.pdf"), width = 10, height = 7)
p <- ggplot(stability_flat, aes(x = Z_mean, y = Z_cv,
                                 size = abs(Z_mean), label = Region)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_text(vjust = -1, size = 3.5) +
  scale_size_continuous(range = c(3, 10), guide = "none") +
  labs(title = "Stability of PGLS results across parameter sweep",
       subtitle = "Higher |Z_mean| = stronger effect; Lower CV = more stable",
       x = "Mean Z-score", y = "Coefficient of Variation (CV) of Z-score") +
  theme_minimal(base_size = 13)
print(p)
dev.off()
cat("Saved: stability_bubble.pdf\n")

# ==========================================
# 6. 控制台摘要
# ==========================================
cat("\n===== Parameter Sweep Complete =====\n")
cat(sprintf("Total combinations evaluated: %d\n", nrow(results_df)))
cat(sprintf("Unique regions: %d\n", length(unique(results_df$Region))))
cat(sprintf("Results saved to: %s\n", output_dir))
cat("Files:\n")
cat("  pgls_param_results.csv  - 完整扫描结果\n")
cat("  stability_summary.csv   - 稳定性汇总 (mean/SD/CV)\n")
cat("  heatmap_zscore.pdf      - Z-score 热力图\n")
cat("  heatmap_r2.pdf          - R² 热力图\n")
cat("  heatmap_pvalue.pdf      - -log10(p) 热力图\n")
cat("  trend_zscore.pdf        - Z-score 趋势线图\n")
cat("  heatmap_rank.pdf        - 回声关联排名热力图 (跨参数稳定性)\n")
cat("  pgls_rank_results.csv   - 排名数据\n")
cat("  rank_summary.csv        - 排名稳定性汇总 (mean/SD/min/max)\n")
cat("  stability_bubble.pdf    - 稳定性气泡图\n")

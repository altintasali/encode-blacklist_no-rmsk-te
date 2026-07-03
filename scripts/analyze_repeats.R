#==============================================================================#
# Project info ---- 
#==============================================================================#
# Project    : RepeatMasker stats
# Data       : RepeatMasker and TE-removed blacklisted regions (mm10)
# Task       : Basic descriptive stats + Fisher Enrichment
# Author     : Ali Altintas 
# Date       : 03.07.2026
#==============================================================================#

rm(list = ls(all = TRUE)) 

#==============================================================================#
# Libraries & Directories ----
#==============================================================================#
if (!require("pacman")) install.packages("pacman")
pacman::p_load_gh("altintasali/aamisc")
pacman::p_load("aamisc", "data.table", "magrittr", "ggplot2", "scales")

# Set paths relative to repo root folder
wd <- aamisc::get_script_path() %>% dirname() %>% dirname() # Steps out of 'scripts/' to repo root
setwd(wd)

# Establish required output directory structures dynamically
dir.create("results/enrichment", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

# Establish console mirroring pipeline into logs folder
sInfo <- file("logs/analyze_repeats.log", open = "wt")
sink(sInfo, type = "output", split = TRUE) # Mirrors code evaluations to console + file log
sink(sInfo, type = "message")

# Open a null device to suppress unwanted Rplots.pdf generation
pdf(NULL)

#==============================================================================#
# Read data ----
#==============================================================================#
rmsk_colnames <- c(
  "bin", "swScore", "milliDiv", "milliDel", "milliIns", 
  "chrom", "start", "end", "genoLeft", "strand", 
  "repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft", "id"
)

rmsk <- data.table::fread(
  "data/raw/rmsk.txt.gz",
  header = FALSE,
  col.names = rmsk_colnames,
  sep = "\t",
  showProgress = TRUE
)
rmsk$repClass <- gsub("\\?", "", rmsk$repClass)

bed_colnames <- c(
  "chrom", "start", "end", "repName", "score", "strand", "repClass", "repFamily"
)

trapped_rmsk <- fread("data/process_temp/mm10_all-rmsk_trapped_in_blacklist.bed", 
                      header = FALSE, col.names = bed_colnames, sep = "\t")
trapped_rmsk$repClass <- gsub("\\?", "", trapped_rmsk$repClass)

trapped_tes <- fread("data/process_temp/mm10_TEs_trapped_in_blacklist.bed", 
                     header = FALSE, col.names = bed_colnames, sep = "\t")
trapped_tes$repClass <- gsub("\\?", "", trapped_tes$repClass)

tes <- fread("data/process_temp/mm10_only-TEs.bed", 
             header = FALSE, col.names = bed_colnames, sep = "\t")
tes$repClass <- gsub("\\?", "", tes$repClass)
te_classes   <- unique(tes$repClass)

#==============================================================================#
# Plot Composition Stats ----
#==============================================================================#
dt1 <- rmsk$repClass %>% table %>% data.table
colnames(dt1) <- c("Class", "Count")
dt1[, Type := "non-TE"]; dt1[Class %in% te_classes, Type := "TE"]; dt1[, Pool := "all"]
total_te_pool <- sum(dt1$Count)

dt2 <- trapped_rmsk$repClass %>% table %>% data.table
colnames(dt2) <- c("Class", "Count")
dt2[, Type := "non-TE"]; dt2[Class %in% te_classes, Type := "TE"]; dt2[, Pool := "black listed"]

dat <- rbind(dt1, dt2)
dat[, Pct := Count/sum(Count)*100, by = Pool]
dat[, Label := paste0(comma(Count), " (", round(Pct, 2), "%)")]

sorted_classes <- dat[Pool == "all"][order(Pct), Class]
dat[, Class := factor(Class, levels = sorted_classes)]

ggplot(dat, aes(x = Pct, y = Class, fill = Type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, show.legend = TRUE) +
  geom_text(aes(label = Label, color = Type), hjust = -0.05, size = 3.2, fontface = "bold") +
  guides(color = "none") +
  scale_x_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.35))) +
  labs(title = "Composition of the mm10 Transposable Element Pool",
       subtitle = paste0("Total Elements Analyzed: ", comma(total_te_pool)),
       x = "Percentage of Total Pool inside Context", y = "Repeat Class", fill = "TE Classification") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
        panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "bold"),
        legend.position = "top", strip.text = element_text(face = "bold", size = 12)) + 
  facet_grid(~Pool) 

ggsave("results/enrichment/mm10_TE_pool_horizontal.png", width = 11, height = 7, dpi = 300)

#==============================================================================#
# Enrichment Analysis ----
#==============================================================================#
stat_dt <- merge(dt1[, .(Class, Genome_Count = Count)], dt2[, .(Class, Blacklist_Count = Count)], by = "Class", all.x = TRUE)
stat_dt[is.na(Blacklist_Count), Blacklist_Count := 0]
stat_dt[, Not_Blacklist_Count := Genome_Count - Blacklist_Count]

total_blacklisted     <- sum(stat_dt$Blacklist_Count)
total_not_blacklisted <- sum(stat_dt$Not_Blacklist_Count)

results <- list()
for (i in 1:nrow(stat_dt)) {
  current_class <- stat_dt$Class[i]
  a <- stat_dt$Blacklist_Count[i]; b <- stat_dt$Not_Blacklist_Count[i]
  c <- total_blacklisted - a; d <- total_not_blacklisted - b
  
  contingency_matrix <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  fisher_result      <- fisher.test(contingency_matrix)
  
  results[[current_class]] <- data.table(
    Class = current_class, Odds_Ratio = as.numeric(fisher_result$estimate),
    Lower_CI = fisher_result$conf.int[1], Upper_CI = fisher_result$conf.int[2], P_Value = fisher_result$p.value
  )
}

enrichment_results <- rbindlist(results)
enrichment_results[, Adj_P_Value := p.adjust(P_Value, method = "BH")]

plot_df <- enrichment_results[Odds_Ratio > 0 & is.finite(Upper_CI)]
plot_df[, Status := "Not Significant"]
plot_df[Odds_Ratio > 1 & Adj_P_Value < 0.05, Status := "Enriched"]
plot_df[Odds_Ratio < 1 & Adj_P_Value < 0.05, Status := "Depleted"]

plot_df[, Class := factor(Class, levels = plot_df[order(Odds_Ratio), Class])]
plot_df[, OR_Label := paste0("OR = ", round(Odds_Ratio, 2))]

message("\n--- Enrichment Statistics Breakdown ---")
print(plot_df[order(-Odds_Ratio)])

#==============================================================================#
# Export Statistical Results ----
#==============================================================================#
message("-> Exporting enrichment statistical tables...")

data.table::fwrite(
  plot_df[order(-Odds_Ratio)], 
  file = "results/enrichment/enrichment_results.tsv.gz", 
  sep = "\t", 
  compress = "gzip"
)

ggplot(plot_df, aes(x = Odds_Ratio, y = Class, color = Status)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI), size = 0.8, fatten = 2.5) +
  #geom_text(aes(x = Upper_CI, label = OR_Label), hjust = -0.2, size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), expand = expansion(mult = c(0.05, 0.25))) + 
  scale_color_manual(values = c("Enriched" = "#D0021B", "Depleted" = "#4A90E2", "Not Significant" = "gray60")) +
  labs(title = "Blacklist Enrichment Across RepeatMasker Classes",
       subtitle = "Fisher's Exact Test Odds Ratios with 95% Confidence Intervals",
       x = "Odds Ratio (Log10 Scale)", y = "Repeat Class", color = "Genomic Status") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
        panel.grid.minor = element_blank(), axis.text.y = element_text(face = "bold"), legend.position = "top")

ggsave("results/enrichment/repeat_blacklist_forest_plot.png", width = 10, height = 6, dpi = 300)

#==============================================================================#
# Session Info ----
#==============================================================================#
message("\n--- Session Environment Info ---")
sessionInfo()

# Close the sink channels safely
sink(type = "message")
sink(type = "output")

library(data.table)
library(qqman)
library(tidyverse)
library(ggplot2)
library(CMplot)
library(dplyr)




# Load the data
# adjust filepath to where you saved the summary stats
diabetes <- read.csv("Gwas Diabetes 2 data.csv")

# Inspect column names
dim(diabetes)
colnames(diabetes)
head(diabetes)
summary(diabetes)

# Histogram of P-values
hist(diabetes$P)

# Log-transform P-values
diabetes$logP <- -log10(diabetes$P)

# Histogram of -log10(P)
ggplot(diabetes, aes(x = logP)) +
  geom_histogram(fill = "steelblue", color = "black", bins = 50) +
  labs(
    title = "Histogram of -log10(P) values (Diabetes Disease)",
    x = expression(-log[10](P)),
    y = "Frequency"
  ) +
  theme_minimal()


# Histogram of SNPs/Chromosomes
# Ensure chromosome order is numeric or natural (1–22, X, Y)
diabetes$CHR <- factor(diabetes$CHR, 
                     levels = c(as.character(1:22), "X", "Y"))

# Histogram of SNPs per Chromosome
# ---- Order chromosomes (1–22, then X) ----

diabetes$CHR <- factor(diabetes$CHR,
                       levels = c(as.character(1:22), "X"))

ggplot(diabetes, aes(x = CHR)) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(
    title = "Histogram of SNPs per Chromosome",
    x = "Chromosome",
    y = "Frequency"
  ) +
  theme_minimal()

# Normalize per chromosome size:
# ============================================================
# SNP Density per Chromosome (Normalized by Chromosome Size)
# ============================================================


# ---- Chromosome sizes (GRCh38, in base pairs) ----
chr_sizes <- data.frame(
  CHR = as.character(c(1:22, "X")),
  size_bp = c(
    248956422, 242193529, 198295559, 190214555, 181538259,
    170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189,
    90338345, 83257441, 80373285, 58617616, 64444167,
    46709983, 50818468, 156040895  # X
  )
)

# ---- Compute SNP count and density ----
snp_density <- diabetes %>%
  filter(!is.na(CHR)) %>%
  mutate(CHR = as.character(CHR)) %>%
  group_by(CHR) %>%
  summarise(
    SNP_count = n(),
    Significant = sum(P < 1e-5),
    .groups = "drop"
  ) %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(
    SNPs_per_Mb = SNP_count / (size_bp / 1e6),
    Sig_density = Significant / (size_bp / 1e6)
  )

# ---- Order chromosomes (1–22, then X) ----
snp_density$CHR <- factor(snp_density$CHR, 
                          levels = c(as.character(1:22), "X"))

# ---- Plot SNP density ----
ggplot(snp_density, aes(x = CHR, y = SNPs_per_Mb, fill = Sig_density > 0)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "blue"),
    name = "Contains Significant SNP (P < 1e-5)"
  ) +
  labs(
    title = "SNP Density per Chromosome (Normalized by Chromosome Size)",
    x = "Chromosome",
    y = "SNPs per Megabase (Mb)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



# ============================================================
#  Manhattan Plot
# ============================================================

# CMplot requires columns in the order: SNP, CHR, BP, P
# So if needed, reorder or subset:

# Clean P-values
diabetes_clean <- diabetes %>%
  mutate(P = as.numeric(P)) %>%   # ensure numeric
  filter(!is.na(P) & P > 0 & P < 1)  # keep valid range

diabetes_cm <- diabetes_clean[, c("SNP", "CHR", "BP", "P")]

# Manhattan plot
CMplot(
  diabetes_cm,
  plot.type = "m",
  col = c("grey30", "skyblue3"),      # alternating colors
  LOG10 = TRUE,
  threshold = c(1e-5, 5e-8),          # suggestive & genome-wide
  threshold.col = c("orange", "red"),
  threshold.lty = c(2, 1),
  main = "Crohn's Disease",
  cex = 0.6,                          # point size
  file = "jpg",                       # output format
  dpi = 300,
  file.output = FALSE, # show in RStudio, don't save
  verbose = TRUE
)
# ============================================================
#  QQ Plot
# ============================================================

# Compute expected and observed
n <- length(diabetes$P)
expected <- -log10(ppoints(n))
observed <- -log10(sort(diabetes$P))

qq_df <- data.frame(expected, observed)

ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(color = "dodgerblue3", size = 1) +
  labs(
    title = "QQ Plot: Diabetes Disease",
    x = "Expected -log10(P)",
    y = "Observed -log10(P)"
  ) +
  theme_minimal()

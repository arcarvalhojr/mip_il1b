# ==============================================================================
# 
# Import raw count data, tidy, run QC and save colData and filtered dds object
# ==============================================================================


# import libraries --------------------------------------------------------

# if not, install DESeq2 package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")


# require libs
library(tidyverse)
library(readxl)
library(janitor)
library(DESeq2)
library(pheatmap)

# import data -------------------------------------------------------------
seq_data <- readr::read_tsv("Data/raw/rna_seq/counts_hisat_RefSeq.tsv")

exp_design <- readxl::read_xlsx("Data/raw/rna_seq/exp_design.xlsx")


# data tyding -------------------------------------------------------------

# transform the seq data into matrix with genes as row names
seq_data <- seq_data %>%
  dplyr::select(-width) %>% 
  tibble::column_to_rownames("...1") %>% 
  as.matrix()

colnames(seq_data) <- sub(".*_", "", colnames(seq_data))

# create col sample_id to match with the id of the seq_data cols
exp_design <- exp_design %>%
  janitor::clean_names() %>% 
  dplyr::mutate(
    sample_id = paste0("S", id)
  )

# re-order seq_data cols for match with sample_id col of the exp_design
seq_data <- seq_data[, exp_design$sample_id]

# total reads per sample
colSums(seq_data)

# create dds object -------------------------------------------------------

# create the colData for
colData <- exp_design %>%
  dplyr::mutate(group = paste0(gestational_day, "_", infection)) %>% 
  dplyr::select(sample_id, gestational_day, infection, group) %>% 
  tibble::column_to_rownames("sample_id")

# check if sample names match in both files
all(rownames(colData) == colnames(seq_data))

# define the reference conditions for each factor
colData$gestational_day <- factor(colData$gestational_day,
                                     levels = c("Gd16.5", "Gd17.5", "Gd18.5"))

colData$infection <- factor(colData$infection, levels = c("no", "yes"))

# dds object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = seq_data,
                                      colData = colData,
                                      design = ~ gestational_day+infection+gestational_day:infection)


# exploratory data analysis -----------------------------------------------

# filter out low count genes
keep <- rowSums(counts(dds) > 20) >= 4
dds <- dds[keep, ]

# rlog transformation
rld <- DESeq2::rlog(dds, blind = TRUE)

# extract vsd count matrix
rld_mat <- assay(rld) 

# calculate variance
rv <- rowVars(rld_mat)

# select the top500 variable genes
top500 <- order(rv, decreasing = TRUE)[1:500]


# ------------------------------------- PCA
# perform PCA on the transposed matrix of data
pca <- prcomp(t(rld_mat[top500, ]))  

# get the "importance" of each pc, stored in the second row
pca_var <- round(summary(pca)$importance[2,] * 100, 2) 

# create df with metadata
df <- cbind(exp_design, pca$x) 

# before plot, lets see which variables drive each pc
anova(lm(PC1 ~ gestational_day + infection, data = df))
anova(lm(PC2 ~ gestational_day + infection, data = df))
anova(lm(PC3 ~ gestational_day + infection, data = df))
anova(lm(PC4 ~ gestational_day + infection, data = df))

# pc1 vs pc4 for gd effect
ggplot(df, aes(x = PC1, y = PC2, 
               color = gestational_day, shape = infection)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 3) +
  labs(
    x = paste0("PC1: ", pca_var["PC1"], "%"),
    y = paste0("PC2: ", pca_var["PC2"], "%")
  ) +
  theme_bw()

# pc1 vs pc3 faceted plot to better visualize the effect of infection
ggplot(df, aes(x = PC1, y = PC3, color = infection)) +
  geom_point(size = 3) +
  facet_wrap(~ gestational_day) +
  labs(
    x = paste0("PC1: ", pca_var["PC1"], "%"),
    y = paste0("PC3: ", pca_var["PC3"], "%")
  ) +
  theme_bw()


# ------------------------------------ distance plot
# get sample dists on the transposed count matrix
sample_dists <- dist(t(rld_mat[top500, ]))

# transform into a matrix
sample_dist_matrix <- as.matrix(sample_dists)

pheatmap(
  sample_dist_matrix,
  annotation_col = colData %>% dplyr::select(-group),
  color = colorRampPalette(
    rev(RColorBrewer::brewer.pal(9, "Blues"))
  )(255)
)


# save data ---------------------------------------------------------------

saveRDS(dds, "Data/processed/dds_filtered.rds")
saveRDS(colData, "Data/processed/colData.rds")


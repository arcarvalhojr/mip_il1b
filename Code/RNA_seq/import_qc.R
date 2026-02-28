#### Exploratory analysis (qc) for the bulk rna-seq data  ####


# import libraries --------------------------------------------------------

# if not, install DESeq2 package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")


# require libs
library(tidyverse)
library(DESeq2)
library(readxl)
library(pheatmap)

# import data -------------------------------------------------------------
seq_data <- read_tsv("Data/raw/rna_seq/counts_hisat_RefSeq.tsv")

exp_design <- read_xlsx("Data/raw/rna_seq/exp_design.xlsx")


# data tyding -------------------------------------------------------------

# transform the seq data into matrix with genes as row names
seq_data <- seq_data %>%
  dplyr::select(-width) %>% 
  dplyr::rename(gene_symbol = ...1) %>%
  tibble::column_to_rownames("gene_symbol") %>% 
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
  dplyr::select(sample_id, gestational_day, infection) %>% 
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
dds <- dds[rowSums(counts(dds)) >= 10, ]

# perform a variance stabilizing transformation
vsd <- DESeq2::vst(dds, blind = T)

# rlog transform
#rld <- DESeq2::rlog(dds, blind = FALSE)


# ------------------------------------- PCA
# extract vsd count matrix
vsd_mat <- assay(vsd) 

# calculate variance
rv <- rowVars(vsd_mat)

# select the top500 variable genes
top500 <- order(rv, decreasing = TRUE)[1:500]

# perform PCA on the transposed matrix of data
pca <- prcomp(t(vsd_mat[top500, ]))  

# create df with metadata
df <- cbind(exp_design, pca$x) 

# before plot, lets see what variables are driven the each pc
anova(lm(PC1 ~ gestational_day + infection, data = df))
anova(lm(PC2 ~ gestational_day + infection, data = df))
anova(lm(PC3 ~ gestational_day + infection, data = df))
anova(lm(PC4 ~ gestational_day + infection, data = df))

# get the "importance" of each pc, stored in the second row
pca_var <- round(summary(pca)$importance[2,] * 100, 2) 

ggplot(df, aes(x = PC1, y = PC4, 
               color = gestational_day, shape = infection)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 3) +
  labs(
    x = paste0("PC1: ", pca_var["PC1"], "%"),
    y = paste0("PC4: ", pca_var["PC4"], "%")
  ) +
  theme_bw()

# facet plot by gd to better visualize the infection
ggplot(df, aes(x = PC1, y = PC3, color = infection)) +
  geom_point(size = 3) +
  facet_wrap(~ gestational_day) +
  labs(
    x = paste0("PC1: ", pca_var["PC1"], "%"),
    y = paste0("PC3: ", pca_var["PC3"], "%")
  ) +
  theme_bw()


colSums(seq_data >= 10)
boxplot(log2(seq_data + 1), outline = FALSE)


# ------------------------------------ distance plot
# get sample dists on the transposed count matrix
sample_dists <- dist(t(assay(vsd)))

# transform into a matrix
sample_dist_matrix <- as.matrix(sample_dists)

pheatmap(
  sample_dist_matrix,
  annotation_col = colData,
  color = colorRampPalette(
    rev(RColorBrewer::brewer.pal(9, "Blues"))
  )(255)
)


# ---------------------------------- distance plot (genes)
# select 40 genes
top_var_mat <- vsd_mat[top500[1:40], ]

# subtract out the means from each row, leaving the variances for each gene
top_var_mat <- top_var_mat - rowMeans(top_var_mat)

pheatmap(
  top_var_mat,
  annotation_col = colData,
  color = colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )(255),
  scale = "none",
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = NA
)


# save data ---------------------------------------------------------------

saveRDS(dds, "Data/processed/dds_filtered.rds")

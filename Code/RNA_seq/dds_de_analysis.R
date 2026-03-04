#### Differential expression analysis for the bulk rna-seq data  ####


# import libraries --------------------------------------------------------

# require libs
library(tidyverse)
library(DESeq2)


# import data -------------------------------------------------------------
dds <- readRDS("Data/processed/dds_filtered.rds")

# differential expression -------------------------------------------------

# run deseq2
dds <- DESeq2::DESeq(dds)

# normalized counts
colSums(DESeq2::counts(dds, normalized=T))

# plot dispersion estimates
plotDispEsts(dds)

# get the names of the comparissons
resultsNames(dds)

# extract the results of the model
# results for the effect of pregnacy for gds vs gd16.5
res_preg_gd17 <- results(dds, name = "gestational_day_Gd17.5_vs_Gd16.5",
                         alpha = 0.05)
res_preg_gd18 <- results(dds, name = "gestational_day_Gd18.5_vs_Gd16.5",
                         alpha = 0.05)

# results for the effect of infection in each gd vs its control
res_inf_gd16 <- results(dds, name = "infection_yes_vs_no", alpha = 0.05)
res_inf_gd17 <- results(dds, contrast = list(c("infection_yes_vs_no",
                                               "gestational_dayGd17.5.infectionyes")),
                        alpha = 0.05)
res_inf_gd18 <- results(dds, contrast = list(c("infection_yes_vs_no",
                                               "gestational_dayGd18.5.infectionyes")),
                        alpha = 0.05)
# for the interaction
res_interaction_gd17 <- results(dds, name = "gestational_dayGd17.5.infectionyes",
                                alpha = 0.05)
res_interaction_gd18 <- results(dds, name = "gestational_dayGd18.5.infectionyes",
                                alpha = 0.05)


teste <- DESeq2::lfcShrink(dds, 
                           contrast = list(c("infection_yes_vs_no",
                                                  "gestational_dayGd18.5.infectionyes")),
                          res = res_inf_gd18, type = "ashr")
teste2 <- DESeq2::lfcShrink(dds, contrast = "infection_yes_vs_no",
                            res = res_inf_gd16, type = "ashr")

# check the results summary
DESeq2::summary(res_preg_gd17)
DESeq2::summary(res_preg_gd18)
DESeq2::summary(res_inf_gd16)
DESeq2::summary(res_inf_gd17)
DESeq2::summary(res_inf_gd18)
DESeq2::summary(res_interaction_gd17)
DESeq2::summary(res_interaction_gd18)

DESeq2::summary(teste)

# create and save results tables ------------------------------------------


# function to create tbl for each result with all genes
make_res_tbl <- function(res) {
  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    dplyr::filter(!is.na(padj)) %>% 
    as_tibble()
}

df_preg_gd17_gd16 <- make_res_tbl(res_preg_gd17)
df_preg_gd18_gd16 <- make_res_tbl(res_preg_gd18)
df_inf_gd16 <- make_res_tbl(res_inf_gd16)
df_inf_gd17_gd16 <- make_res_tbl(res_inf_gd17)
df_inf_gd18_gd16 <- make_res_tbl(res_inf_gd18)
df_inter_gd17 <- make_res_tbl(res_interaction_gd17)
df_inter_gd18 <- make_res_tbl(res_interaction_gd18)

df_teste <- make_res_tbl(teste)

# function to create a filtered tbl for the results
make_res_filt_tbl <- function(tbl) {
  tbl %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>% 
    dplyr::arrange(desc(padj))
}

df_inf_gd16_filt <- make_res_filt_tbl(df_inf_gd16)
df_inf_gd18_gd16_filt <- make_res_filt_tbl(df_inf_gd18_gd16)



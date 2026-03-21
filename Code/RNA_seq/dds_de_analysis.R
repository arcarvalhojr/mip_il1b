# ==============================================================================
# 
# Perform Differential Expression analysis.
#
# ==============================================================================


# import libraries --------------------------------------------------------

# require libs
library(tidyverse)
library(DESeq2)
library(purrr)
library(ashr)
library(writexl)

# import data -------------------------------------------------------------
dds <- readRDS("Data/processed/dds_filtered.rds")
colData <- readRDS("Data/processed/colData.rds")


# DE analysis -------------------------------------------------------------

# run deseq function 
dds <- DESeq2::DESeq(dds)

# plot dispersion estimates
DESeq2::plotDispEsts(dds)

# get comparisons name
DESeq2::resultsNames(dds)


# extract results of the model
# results for the effect of pregnacy for gds vs gd16.5
res_gd17_vs_gd16_ctrl <- DESeq2::results(
  dds, name = "gestational_day_Gd17.5_vs_Gd16.5",
  alpha    = 0.05
)

res_gd18_vs_gd16_ctrl <- DESeq2::results(
  dds,
  name = "gestational_day_Gd18.5_vs_Gd16.5",
  alpha    = 0.05
)

# results for the effect of infection in each gd vs its control
res_gd16_inf <- DESeq2::results(
  dds,
  name = "infection_yes_vs_no",
  alpha    = 0.05
)

res_gd17_inf <- DESeq2::results(
  dds,
  contrast = list(c("infection_yes_vs_no",
                    "gestational_dayGd17.5.infectionyes")),
  alpha    = 0.05
)

res_gd18_inf <- DESeq2::results(
  dds,
  contrast = list(c("infection_yes_vs_no",
                    "gestational_dayGd18.5.infectionyes")),
  alpha    = 0.05
)

# results for the interaction
res_gd17_vs_gd16_int <- DESeq2::results(
  dds,
  name = "gestational_dayGd17.5.infectionyes",
  alpha = 0.05
)

res_gd18_vs_gd16_int <- DESeq2::results(
  dds,
  name = "gestational_dayGd18.5.infectionyes",
  alpha = 0.05
)

# results list
results_list <- list(
  gd16_inf = res_gd16_inf,
  gd17_inf = res_gd17_inf,
  gd18_inf = res_gd18_inf,
  gd17_vs_gd16_int = res_gd17_vs_gd16_int,
  gd18_vs_gd16_int = res_gd18_vs_gd16_int,
  gd17_vs_gd16_ctrl = res_gd17_vs_gd16_ctrl,
  gd18_vs_gd16_ctrl = res_gd18_vs_gd16_ctrl
)

# results dummary
purrr::walk2(results_list, names(results_list), ~ {
  cat("\n══", .y, "══\n")
  print(DESeq2::summary(.x))
})

saveRDS(results_list, "Data/processed/deseq2_results.rds")

# shrinkage ---------------------------------------------------------------

# shrink results for latter visualizations
shrink_result <- function(contrast, res) {
  DESeq2::lfcShrink(
    dds,
    contrast = contrast,
    type     = "ashr",
    res      = res
  )
}

res_gd16_inf_shr <- shrink_result(c("infection_yes_vs_no"), res_gd16_inf)

res_gd18_inf_shr <- shrink_result(
  c("infection_yes_vs_no", "gestational_dayGd17.5.infectionyes"),
  res_gd18_inf
)

res_gd17_ctrl_shr <- shrink_result(
  c("gestational_day_Gd17.5_vs_Gd16.5"), res_gd17_vs_gd16_ctrl
)

res_gd18_ctrl_shr <- shrink_result(
  c("gestational_day_Gd18.5_vs_Gd16.5"), res_gd18_vs_gd16_ctrl
)

res_gd18_int_shr <- shrink_result(
  c("gestational_dayGd18.5.infectionyes"), res_gd18_vs_gd16_int
)

# shrank results list
results_shr <- list(
  gd16_inf = res_gd16_inf_shr,
  gd18_inf = res_gd18_int_shr,
  gd18_vs_gd16_int = res_gd18_int_shr,
  gd18_vs_gd16_ctrl = res_gd18_ctrl_shr,
  gd17_vs_gd16_ctrl = res_gd17_ctrl_shr
)

saveRDS(results_shr, "Data/processed/deseq2_results_shrunk.rds")

# save results ------------------------------------------------------------

# create results table
format_results <- function(res, name) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(padj)) %>% 
    dplyr::mutate(
      contrast  = name,
      sig = padj < 0.05 & abs(log2FoldChange) > 0.58,
      direction = dplyr::case_when(
        sig & log2FoldChange > 0 ~ "up",
        sig & log2FoldChange < 0 ~ "down",
        TRUE                     ~ "ns"
      )
    ) %>%
    dplyr::arrange(padj)
}

tables_list <- purrr::imap(results_list, format_results)
degs_list <- purrr::map(tables_list, ~ dplyr::filter(.x, sig))


writexl::write_xlsx(
  tables_list,
  "Code/RNA_seq/results/deseq2_results_all_contrasts.xlsx"
)

writexl::write_xlsx(
  degs_list,
  "Code/RNA_seq/results/deseq2_DEGs.xlsx"
)

# packages + data

library(tidyverse)

bulk_ccrcc <- readxl::read_excel("data/ccRCC/41591_2020_839_MOESM2_ESM.xlsx", 
                                 sheet = "S4A_RNA_Expression",
                                 skip = 1) 

metadata <- readxl::read_excel("data/ccRCC/41591_2020_839_MOESM2_ESM.xlsx",
                               sheet = "S1_Clinical_and_Immune_Data",
                               skip = 1)


# lets semi-randomly select 

# RCC25-1048 (mtor) - RNA ID: EA639156

# and

# RCC25-516 (nivo) - RNA ID: P66425-08F-Run1_S7_L001
# these are both ~60yo men

bulk_ccrcc_two <- bulk_ccrcc |>
  select(gene_name, `EA639156`, `P66425-08F-Run1_S7_L001`)

genes_i_need <- read_csv("data/ccRCC/test_seq.csv") |>
  select(Gene)


bulk_ccrcc_genes_i_need <- genes_i_need |>
  left_join(bulk_ccrcc_two, by = c("Gene" = "gene_name")) |>
  rename("nivo" = "P66425-08F-Run1_S7_L001",
         "mtor_inhib" = "EA639156") |>
  drop_na()

write_csv(bulk_ccrcc_genes_i_need, file = "data/ccRCC/ccrcc_mtorvnivo_mapletest.csv")
  



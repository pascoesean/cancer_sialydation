# load packages
################################################################################

library(tidyverse)
library(Seurat)


################################################################################
################################################################################
# KRISHNA ANALYSIS 
################################################################################
################################################################################

#kris_data <- readRDS("data/ccRCC/krishna/ccRCC_6pat_Seurat.rds")

#kris_data_new <- SeuratObject::UpdateSeuratObject(kris_data)
#SeuratObject::SaveSeuratRds(kris_data_new, "data/ccRCC/krishna/kris_data_new.rds")

kris_data_new <- SeuratObject::LoadSeuratRds("data/ccRCC/krishna/kris_data_new.rds")

kris_annot <- read_delim("data/ccRCC/krishna/ccRCC_6pat_cell_annotations.txt", delim = "\t")

# need to get from ensembl ids

DimPlot(kris_data_new, reduction = "umap")

# maybe st3gal4 in pDCs???
VlnPlot(kris_data_new, features = c("ENSG00000110080", # st3gal4
                                    "ENSG00000126091", # st3gal3
                                    "ENSG00000064225", # st3gal6
                                    "ENSG00000177575", # cd163
                                    "ENSG00000026508", # cd44
                                    "ENSG00000129450"  # siglec9
                                    ))


# using which cells to look at those

# look in sample2; need tumor. can compare to healthy. im going to bed.
sialcd44 <- WhichCells(kris_data_new, expression = ENSG00000026508 > 1 &  ENSG00000107159 > .1 & (ENSG00000110080 > .1 | ENSG00000126091 > .1 | ENSG00000064225 > .1))

sig9tams <- WhichCells(kris_data_new, expression = ENSG00000177575 > 1 &  ENSG00000129450 > .1)
sig9s <- WhichCells(kris_data_new, expression = ENSG00000129450 > .1)

# lets literally just look at numbers of these in each condition
total_cells <- FetchData(kris_data_new, vars = c("Sample2", "Sample")) |>
  group_by(Sample2) |>
  summarize(n_total = n())

sialcd44_anno <- FetchData(kris_data_new, vars = c("Sample2"), cells = sialcd44)|>
  group_by(Sample2) |>
  summarize(n_sialcd44s = n())

sig9s_anno <- FetchData(kris_data_new, vars = c("Sample2"), cells = sig9s)|>
  group_by(Sample2) |>
  summarize(n_sig9s = n())

sig9tams_anno <- FetchData(kris_data_new, vars = c("Sample2"), cells = sig9tams)|>
  group_by(Sample2) |>
  summarize(n_sig9tams = n())

# need ntotals by patient number 
totals <- total_cells |>
  group_by(Sample2) |>
  summarize(total_bycat = sum(n_total))

total_cells_joined <- total_cells |>
  left_join(y = sialcd44_anno) |>
  left_join(y = sig9s_anno) |>
  left_join(y = sig9tams_anno) |>
  left_join(y = totals) |>
  mutate(prop_sialcd44 = n_sialcd44s/total_bycat,
         prop_sig9s = n_sig9s/total_bycat,
         prop_sig9tams = n_sig9tams/total_bycat) |>
  separate_wider_delim(Sample2, delim = "_", names = c("patient", "sample_area")) |>
  mutate(sample_area = case_when(
    sample_area == "Upper" ~ "Near",
    sample_area == "Lower" ~ "Far",
    sample_area == "SupraLateral" ~ "Near",
    sample_area == "LowerLateral" ~ "Far",
    sample_area == "LowerMedial" ~ "Center",
    sample_area == "Medial" ~ "Near", 
    sample_area == "Lateral" ~ "Far",
    TRUE ~ sample_area
  ))

# interesting....
total_cells_joined |>
  pivot_longer(cols = starts_with("prop_"), names_to = "subset", values_to = "proportion", names_prefix = "prop_") |>
  ggplot(aes(x = patient, y = proportion)) +
  geom_col() +
  facet_grid(rows = vars(subset), cols = vars(sample_area), scales = "free_y") +
  #scale_fill_brewer(palette = "Set3") +
  ggpubr::theme_pubr() +
  labs(y = "proportion of total cells in this category") + 
  theme(axis.text.x = element_text(angle = 45, hjust=  1, vjust = 1))

total_cells_joined |>
  filter(patient != "UT1" & patient != "UT2") |>
  filter(sample_area != "PBMC" & sample_area != "LymphNode") |>
  #pivot_longer(cols = starts_with("prop_"), names_to = "subset", values_to = "proportion", names_prefix = "prop_") |>
  ggplot(aes(x = prop_sialcd44, y = prop_sig9s, color = patient, shape = sample_area)) +
  geom_point(size = 4) +
  ggpubr::theme_pubr() +
  labs(x = "Proportion with Sialylated CD44", y = "Proportion with Siglec9") 
  #facet_wrap(~patient)
  



################################################################################
################################################################################
# BI ANALYSIS ----
################################################################################
################################################################################

normalized_counts <- data.table::fread("data/ccRCC/bi/ccRCC_scRNASeq_NormalizedCounts.txt.gz", data.table = FALSE)

meta_data <- read_delim(file = "data/ccRCC/bi/Final_SCP_Metadata.txt", delim = "\t") |>
  filter(NAME != "TYPE") |> # not sure what this row is but i know i don't want it
  column_to_rownames("NAME") |>
  as.data.frame(row.names = rownames(~x))

# If I just want to repeat above analysis; I don't actually have to fuck around 
# with the seurat objects

bi_combo <- normalized_counts |>
  filter(GENE %in% c("SIGLEC9", "CD44", "ST3GAL4", "ST3GAL3", "ST3GAL6")) |>
  column_to_rownames("GENE") |>
  t() |>
  as_tibble(rownames = "NAME") |>
  left_join(y = meta_data) 

saveRDS(bi_combo, file = "data/ccRCC/bi/bi_combo.rds")
# for each individual, want to plot % of tumor cells expressing CD44 + sialtransferases
# as well as % other cells expressing siglec9


## NO LINEAGE ----
################################################################################

# might keep both tumor and not tumor here and then subset on later
sialcd44 <- bi_combo |>
  filter(CD44 > 1 & (ST3GAL3 > .1 | ST3GAL4 > .1 | ST3GAL6 > .1))

sig9s <- bi_combo |>
  filter(SIGLEC9 > .1)

sialcd44_counts <- sialcd44 |>
  group_by(donor_id, ICB_Response, organ__ontology_label) |>
  summarize(num_sialcd44 = n())

sig9_counts <- sig9s |>
  group_by(donor_id, ICB_Response, organ__ontology_label) |>
  summarize(num_sig9 = n())

total_counts <- bi_combo |>
  group_by(donor_id, ICB_Response, organ__ontology_label) |>
  summarize(total_bycat = n())


bi_props <- total_counts |>
  full_join(y = sialcd44_counts) |>
  full_join(y = sig9_counts) |>
  # introduced NA values here are just zeros (there were no counts of that)
  mutate_all(~replace(., is.na(.), 0)) |> 
  mutate(prop_sig9 = num_sig9/total_bycat,
         prop_sialcd44 = num_sialcd44/total_bycat) 

bi_props |>
  mutate(ICB_Response = case_when(
    ICB_Response == "ICB_PR" ~ "Partial Response",
    ICB_Response == "ICB_SD" ~ "Stable Disease",
    ICB_Response == "ICB_PD" ~ "Progressive Disease",
    ICB_Response == "ICB_NE" ~ "Not Evaluable",
    ICB_Response == "NoICB" ~ "No ICB Treatment"
  )) |>
  ggplot(aes(x = prop_sialcd44, y = prop_sig9, color = ICB_Response, shape = organ__ontology_label)) +
  geom_point(size = 4) +
  ggpubr::theme_pubr() +
  labs(x = "Proportion with Sialylated CD44", y = "Proportion with Siglec9") 

## WITH LINEAGE ----
################################################################################

sialcd44 <- bi_combo |>
  filter(CD44 > 1 & (ST3GAL3 > .1 | ST3GAL4 > .1 | ST3GAL6 > .1) & Lineage == "Putative Tumor")

sig9s <- bi_combo |>
  filter(SIGLEC9 > .1 & Lineage == "Myeloid")

sialcd44_counts <- sialcd44 |>
  group_by(donor_id, ICB_Response, organ__ontology_label) |>
  summarize(num_sialcd44 = n())

sig9_counts <- sig9s |>
  group_by(donor_id, ICB_Response, organ__ontology_label) |>
  summarize(num_sig9 = n())

total_counts_myeloid <- bi_combo |>
  filter(Lineage == "Myeloid") |>
  group_by(donor_id, ICB_Response, organ__ontology_label) |>
  summarize(total_mye_bycat = n())

total_counts_tumor <- bi_combo |>
  filter(Lineage == "Putative Tumor") |>
  group_by(donor_id, ICB_Response, organ__ontology_label) |>
  summarize(total_tum_bycat = n())

bi_props <- total_counts_myeloid |>
  left_join(y = total_counts_tumor) |>
  full_join(y = sialcd44_counts) |>
  full_join(y = sig9_counts) |>
  # introduced NA values here are just zeros (there were no counts of that)
  mutate_all(~replace(., is.na(.), 0)) |> 
  mutate(prop_mye_sig9 = num_sig9/total_mye_bycat,
         prop_tum_sialcd44 = num_sialcd44/total_tum_bycat) 

bi_props |>
  mutate(ICB_Response = case_when(
    ICB_Response == "ICB_PR" ~ "Partial Response",
    ICB_Response == "ICB_SD" ~ "Stable Disease",
    ICB_Response == "ICB_PD" ~ "Progressive Disease",
    ICB_Response == "ICB_NE" ~ "Not Evaluable",
    ICB_Response == "NoICB" ~ "No ICB Treatment"
  )) |>
  ggplot(aes(x = prop_tum_sialcd44, y = prop_mye_sig9, color = ICB_Response, shape = organ__ontology_label)) +
  geom_point(size = 4) +
  ggpubr::theme_pubr() +
  labs(x = "Proportion of Tumor cells with Sialylated CD44", y = "Proportion of Myeloid cells with Siglec9") 





## Seurat Stuff ----
################################################################################

# none of this really works. and that doesn't bother me anymore #free

ncounts_assay <- CreateAssayObject(data = normalized_counts |> 
                                     column_to_rownames("GENE") |> 
                                     as.matrix()
                                   )

bi_seurat <- CreateSeuratObject(ncounts_assay, meta.data = meta_data)

bi_seurat <- readRDS("data/ccRCC/bi/bi_seurat.rds")

bi.genes <- rownames(bi_seurat)
bi_seurat <- ScaleData(bi_seurat, features = bi.genes)

bi_seurat <- FindVariableFeatures(object = bi_seurat)
bi_seurat <- RunPCA(bi_seurat, features = VariableFeatures(object = bi_seurat))
ElbowPlot(bi_seurat, ndims = 50)

DimHeatmap(bi_seurat, dims = 20:40, cells = 500, balanced = TRUE)

bi_seurat <- FindNeighbors(bi_seurat, dims = 1:40)
bi_seurat <- FindClusters(bi_seurat, resolution = 0.05)

bi_seurat <- RunUMAP(bi_seurat, dims = 1:40)


#SeuratObject::SaveSeuratRds(bi_seurat, "data/ccRCC/bi/bi_seurat.rds")



# using which cells to look at those
sialcd44 <- WhichCells(bi_seurat, expression = CD44 > 1 &  ("ST3GAL5" > 1 | "ST3GAL4" > 1 | "ST3GAL6" > 1))

sig9tams <- WhichCells(bi_seurat, expression = CD163 > 1 &  SIGLEC9 > 1)
sig9s <- WhichCells(bi_seurat, expression = SIGLEC9 > 1)

# lets literally just look at numbers of these in each condition
total_cells <- FetchData(bi_seurat, vars = c("ICB_Response", "donor_id")) |>
  group_by(ICB_Response, donor_id) |>
  summarize(n_total = n())

sialcd44_anno <- FetchData(bi_seurat, vars = c("ICB_Response", "donor_id"), cells = sialcd44)|>
  group_by(ICB_Response, donor_id) |>
  summarize(n_sialcd44s = n())

sig9s_anno <- FetchData(bi_seurat, vars = c("ICB_Response", "donor_id"), cells = sig9s)|>
  group_by(ICB_Response, donor_id) |>
  summarize(n_sig9s = n())

sig9tams_anno <- FetchData(bi_seurat, vars = c("ICB_Response", "donor_id"), cells = sig9tams)|>
  group_by(ICB_Response, donor_id) |>
  summarize(n_sig9tams = n())

# need ntotals by patient number 
totals <- total_cells |>
  group_by(ICB_Response) |>
  summarize(total_bycat = sum(n_total))

total_cells_joined <- total_cells |>
  left_join(y = sialcd44_anno) |>
  left_join(y = sig9s_anno) |>
  left_join(y = sig9tams_anno) |>
  mutate(total_forcat = case_when(
    ICB_Response == "ICB_NE" ~ 280,
    ICB_Response == "ICB_PD" ~ 2470,
    ICB_Response == "ICB_PR" ~ 11178,
    ICB_Response == "ICB_SD" ~ 3744,
    ICB_Response == "NoICB" ~ 16654
  )) |>
  mutate(prop_sialcd44 = n_sialcd44s/total_forcat,
         prop_sig9s = n_sig9s/total_forcat,
         prop_sig9tams = n_sig9tams/total_forcat)

# interesting....
total_cells_joined |>
  filter(ICB_Response != "ICB_NE") |>
  select(ICB_Response, donor_id, starts_with("prop_")) |>
  pivot_longer(cols = starts_with("prop_"), names_to = "subset", values_to = "proportion", names_prefix = "prop_") |>
  ggplot(aes(x = ICB_Response, y = proportion, fill = donor_id)) +
  geom_col() +
  facet_wrap(~subset, scales = "free_y") +
  #scale_fill_brewer(palette = "Set3") +
  ggpubr::theme_pubr() +
  labs(y = "proportion of total cells in this category") + 
  theme(axis.text.x = element_text(angle = 45, hjust=  1, vjust = 1))

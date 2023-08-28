# UCLA-ASD_KCNI
```{r}
library(dplyr)
library(tidyverse)
library(edgeR)
library(markerGeneProfile) 
library(matrixStats)
library(cowplot)
library(broom)
library(knitr)
library(ggpubr)
library(biomaRt)
theme_set(theme_cowplot())
```
Combining results from UCLA
```{r}
#Create vectors with gene and isoform filenames
csv_dir <- "/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/UCLA_ASD/RNA/"
gene_csv <- list.files(path = csv_dir, pattern = "\\genes.results$")
isoform_csv <- list.files(path = csv_dir, pattern = "\\isoforms.results$")

#Loop through gene files
setwd(csv_dir)
library(dplyr)
library(tools)

output_dir <- "/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_count_matrix <- NULL

for (i in 1:length(gene_csv)) {
  curr_file <- gene_csv[i] # get current file 
  sample_matrix <- read.delim(file.path(csv_dir, curr_file)) # intake current file table
  sample_name <- gsub("\\..*", "", curr_file) # get sample name (approximate)
  colnames(sample_matrix)[5] <- sample_name # set sample name over expected counts column
  sample_matrix <- sample_matrix[, c("gene_id", sample_name)] # take columns of interest
  if (is.null(gene_count_matrix)) {
    gene_count_matrix <- sample_matrix
  } else {
    gene_count_matrix <- left_join(gene_count_matrix, sample_matrix, by = "gene_id")
  }
}

row.names(gene_count_matrix) <- gene_count_matrix$gene_id # set row names as gene names
gene_count_matrix <- gene_count_matrix[-1] # remove gene name column to make object pure matrix

output_file <- file.path(output_dir, "UCLA_count_matrix.csv")
write.csv(gene_count_matrix, file = output_file, row.names = TRUE)

```
Load and format saved matrices and format for ease of access
```{r}
setwd("/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/")
UCLA_matrix = read.csv("UCLA_count_matrix.csv")

#names(UCLA_matrix) = gsub("X", '', names(UCLA_matrix))
names(UCLA_matrix)[1] = "gene_symbol"

#Rename Ensembl id to be compatible with gene names
UCLA_matrix$gene_symbol = gsub("\\..*", '', UCLA_matrix$gene_symbol)

#Make EnsemblIDs row names
row.names(UCLA_matrix) = UCLA_matrix$gene_symbol
UCLA_matrix = UCLA_matrix[,-1]
```
Import Metadata
```{r}
UCLA_metadata = read.csv("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/UCLA_ASD/Metadata/rnaSeq_UCLA-ASD_metadata.csv")
psychencode_metadata = read.csv(("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))
psychencode_metadata<- psychencode_metadata %>% rename(specimenID= individualID)

merge(x = UCLA_metadata, y = psychencode_metadata, by = "specimenID")
#Remove Sample from specimenID
UCLA_metadata$specimenID <- gsub("^Sample_", "", UCLA_metadata$specimenID)

# Separate the 'specimenID' column into two separate columns 'part1' and 'part2'
UCLA_metadata <- separate(UCLA_metadata,specimenID, into= c("part1", "part2"), sep="_", extra= "merge")

# Rename 'part1' to 'specimenID'
UCLA_metadata<- UCLA_metadata %>%
  rename(specimenID = part1)

#Rename "part2" as 'Brainarea'
UCLA_metadata<- UCLA_metadata %>%
  rename(Brainarea = part2)

#Split the brain regions into the Brodmann Areas and the following part:
UCLA_metadata <- separate(UCLA_metadata, Brainarea, into = c("part1", "part2"), sep="_|-", extra= "merge")

# Rename 'part1' to 'Brain_region'
UCLA_metadata<- UCLA_metadata %>%
  rename(Brain_region = part1)

write.csv(UCLA_metadata, "/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/UCLA_metadata.csv")

#psychencode_metadata = read.csv(("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))
UCLA_merged_metadata = left_join(x = UCLA_metadata, y = psychencode_metadata, by = "specimenID")

# use these brain regions for analysis: BA9, BA41
UCLA_merged_metadata_regions_filtered = UCLA_merged_metadata %>% filter(Brain_region == 'ba9' |
                                                                        Brain_region == 'ba41')
write.csv(UCLA_merged_metadata_regions_filtered,"/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/UCLA_merged_metadata_regions_filtered.csv")

```
Normalize and Process 
```{r}
library(edgeR)
gene_ids <- UCLA_matrix[, 1]
expression_matrix <- as.matrix(UCLA_matrix[, -1])
#Convert matrices to counts per million
UCLA_cpm = cpm(expression_matrix, log = TRUE, prior.count = 0.1)


#Remove genes with low standard deviations
UCLA_sds = rowSds(UCLA_cpm, na.rm = T)

UCLA_matrix = UCLA_cpm[UCLA_sds > 0.1, ] %>% as.data.frame() %>% rownames_to_column( var = "gene_symbol")
```
Converting Ensembl ID to Gene Name
```{r}
library(biomaRt)
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = mart) 
ensembl_to_gene = (data.frame(ensembl$ensembl_gene_id, ensembl$hgnc_symbol))
names(ensembl_to_gene) = c("gene_symbol", "gene_name")

#remove duplicates
ensembl_to_gene = ensembl_to_gene[!duplicated(ensembl_to_gene[,1]),]

UCLA_matrix = merge(x=UCLA_matrix, y=ensembl_to_gene, by = "gene_symbol", all.x = T)
```
Cell-Type Proportion Estimation:
```{r}
sonny_markers = read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
colnames(sonny_markers) = colnames(sonny_markers) %>% make.names() %>% tolower()
```

```{r}
hgnc_mapping = read_tsv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt')
sonny_hgnc_merged_markers = left_join(sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id), 
                                      hgnc_mapping %>% distinct(entrez_id, .keep_all = T)%>% 
                                        dplyr::select(entrez_id, ensembl_gene_id) %>% 
                                        dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>% 
  group_by(subclass) %>% 
  arrange(subclass, -average.log.fold.change) %>% 
  ungroup()
```


```{r}
# get ensembl list of markers
new_markers = sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types = new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique
new_marker_list <- lapply(new_cell_types, function(cell_type){
  return(new_markers %>%
           filter(subclass == cell_type, ensembl_id %in% unique(c(UCLA_matrix$gene_symbol))) %>%
           pull(ensembl_id)
         )
})

names(new_marker_list) = c('Astrocyte', 'Endothelial', 'Exc_IT', 'Exc_L4_IT', 'Exc_L5_ET', 'Exc_L5/6_IT_Car3', 'Exc_L5/6_NP', 'Exc_L6_CT', 'Exc_L6b', 'Inh_LAMP5', 'Microglia', 'Oligodendrocyte', 'OPC', 'Inh_PAX6', 'Pericyte', 'Inh_PVALB', 'Inh_SST', 'Inh_VIP', 'VLMC')
print(new_cell_types)
```
UCLA:  mgpEstimate to get cell type proportions
```{r}
#Estimation counts NULL as a gene - so remove all NULL gene names to avoid duplicates
UCLA_genes_only = UCLA_matrix %>% subset(gene_symbol != "")
if(length(which(duplicated(UCLA_genes_only$gene_symbol))) != 0){
  UCLA_genes_only = UCLA_genes_only[-which(duplicated(UCLA_genes_only$gene_symbol)),]
}

#Remove ensembl_ID and move gene names to first column 
UCLA_genes_only[,1] = UCLA_genes_only$gene_symbol
colnames(UCLA_genes_only)[1] = "gene_name"
UCLA_genes_only = UCLA_genes_only[, -(ncol(UCLA_genes_only))]

UCLA_estimations =  mgpEstimate(
  exprData = UCLA_genes_only,
  genes = new_marker_list,
  geneColName = 'gene_name',
  outlierSampleRemove = F, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = FALSE)

#Coerce estimations list into data frame 
UCLA_estimations_scaled = UCLA_estimations$estimates %>% as.data.frame() %>% scale() %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")

# Separate the 'specimenID' column into two separate columns 'part1' and 'part2'
UCLA_estimations_scaled <- separate(UCLA_estimations_scaled,specimenID, into= c("part1", "part2"), sep="_", extra= "merge")

# Rename 'part1' to 'specimenID'
UCLA_estimations_scaled<- UCLA_estimations_scaled %>%
  rename(specimenID = part1)

#Rename "part2" as 'region'
UCLA_estimations_scaled<- UCLA_estimations_scaled %>%
  rename(region = part2)

UCLA_estimations_scaled <- separate(UCLA_estimations_scaled, region, into = c("part1", "part2"), sep="_", extra= "merge")

#Rename "part1" as 'Brain_region'
UCLA_estimations_scaled<- UCLA_estimations_scaled %>%
  rename(Brain_region = part1)

# use these brain regions for analysis: BA9, BA41
UCLA_estimations_regions_filtered = UCLA_estimations_scaled %>% filter(Brain_region == 'ba41' |
                                                                         Brain_region == 'ba9')

#merge these only after finding the above brain regions
#Merge cell type proportions with sample metadata
UCLA_estimations_metadata = inner_join(UCLA_merged_metadata_regions_filtered %>% mutate(specimenID = make.names(specimenID)),
                                      UCLA_estimations_regions_filtered %>% mutate(specimenID = make.names(specimenID)))

csv_dir <- "/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/UCLA_ASD/Metadata"
#Check for individuals with the 15q duplication mutation:
individual_UCLA_ASD_metadata = read.csv("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/UCLA_ASD/Metadata/individual_UCLA-ASD_metadata.csv")

UCLA_estimations_metadata = left_join(UCLA_estimations_metadata, individual_UCLA_ASD_metadata %>% dplyr::select(individualID, primaryDiagnosisDetails), by = c("specimenID" = "individualID")) 
                                     
#Remove '+' from ageDeath for modelling
UCLA_estimations_metadata$ageDeath = as.numeric(gsub("[+]", "", UCLA_estimations_metadata$ageDeath))

write.csv(UCLA_estimations_metadata, "/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/UCLA_cell_prop_estimates.csv")

#Merge cell type proportions with sample metadata
UCLA_estimations_metadata = inner_join(UCLA_metadata %>% mutate(specimenID = make.names(specimenID)) %>% distinct(specimenID, .keep_all =T), 
                                       UCLA_estimations_scaled %>% mutate(specimenID = make.names(specimenID)))
names(UCLA_estimations_metadata)[1] = "gene_symbol"

# use these brain regions for analysis: BA9, BA41
UCLA_estimations_metadata = UCLA_estimations_metadata %>% filter(Brain_region == 'ba41' |
                                                                         Brain_region == 'ba9')

write.csv(UCLA_estimations_metadata, "/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/UCLA_cell_prop_estimates.csv")
```

#Data Visualization:
 Loading Packages:
```{r}
library(dplyr)
library(tidyverse)
library(edgeR)
library(markerGeneProfile) 
library(matrixStats)
library(cowplot)
library(broom)
library(knitr)
library(ggpubr)
library(biomaRt)
library(ggrepel)
library(patchwork)
library(ggsignif)
library(modelr)
library(ggbeeswarm)
theme_set(theme_classic2())
#Colour palette
cbPalette = c("#56B4E9", "#009E73","#E69F00", "#0072B2", "#D55E00", "#CC79A7","#000000","#F0E442")

setwd("/external/rprshnas01/kcni/aaleem/KCNI_data_analyses")

```
Loading data:
```{r}

UCLA_matrix = read.csv("UCLA_count_matrix.csv")
UCLA_metadata = read.csv("/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/UCLA_metadata.csv")
psychencode_metadata = read.csv(("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))
UCLA_metadata = UCLA_metadata %>% cross_join(psychencode_metadata)
UCLA_merged_metadata_regions_filtered = read.csv("/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/UCLA_merged_metadata_regions_filtered.csv")
individual_UCLA_ASD_metadata = read.csv("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/UCLA_ASD/Metadata/individual_UCLA-ASD_metadata.csv")
UCLA_estimations_metadata = read.csv(("/external/rprshnas01/kcni/aaleem/KCNI_data_analyses/UCLA_cell_prop_estimates.csv"))
UCLA_estimations_metadata$primaryDiagnosis = factor(UCLA_estimations_metadata$primaryDiagnosis, levels =c("control", "Autism Spectrum Disorder"))
UCLA_estimations_metadata = UCLA_estimations_metadata  %>% filter(!is.na(individualIdSource), !is.na(primaryDiagnosis), ageDeath >= 20) 
UCLA_estimations_metadata_long = UCLA_estimations_metadata %>% pivot_longer(cols = Astrocyte:VLMC, names_to = 'cell_type', values_to = 'rel_prop')

UCLA_estimations_metadata_long$ageDeath = UCLA_estimations_metadata_long$ageDeath %>% as.numeric()
#Remove '+' from ageDeath for modelling
UCLA_estimations_metadata_long$ageDeath = as.numeric(gsub("[+]", "", UCLA_estimations_metadata_long$ageDeath))
UCLA_estimations_metadata_long = UCLA_estimations_metadata_long %>% filter(primaryDiagnosis %in% c("Autism Spectrum Disorder", "control"))
UCLA_estimations_metadata_long$primaryDiagnosis = UCLA_estimations_metadata_long$primaryDiagnosis %>% factor(levels = c("control", "Autism Spectrum Disorder"))
```
Generating figures:
```{r}
#Generating plots to see cell-type estimations in BA9 and BA41
UCLA_estimations_metadata %>% ggplot(aes(x = primaryDiagnosis, y = Inh_PVALB)) + geom_boxplot() + facet_wrap(~Brain_region)

UCLA_estimations_metadata %>% ggplot(aes(x = primaryDiagnosis, y = Inh_SST)) + geom_boxplot() + facet_wrap(~Brain_region)

UCLA_estimations_metadata %>% ggplot(aes(x = primaryDiagnosis, y = Exc_L4_IT)) + geom_boxplot() + facet_wrap(~Brain_region)

UCLA_estimations_metadata %>% ggplot(aes(x = primaryDiagnosis, y = OPC)) + geom_boxplot() + facet_wrap(~Brain_region)

UCLA_estimations_metadata %>% ggplot(aes(x = primaryDiagnosis, y = Endothelial)) + geom_boxplot() + facet_wrap(~Brain_region)

#combined linear models
combined_lms = UCLA_estimations_metadata_long %>%
  # group stacked data by cell_type
  group_by(cell_type, Brain_region) %>%
  # fit all the cell_type_prop data according to the model
  # using the broom package to tidy the results
  do(tidy(lm(scale(rel_prop) ~ scale(RIN) + scale(ageDeath)+ reportedGender + primaryDiagnosis, data = .))) %>% 
  
ungroup() %>%
  mutate(padj = p.adjust(p.value, method = 'BH')) %>%
  # add cell class labels
  mutate(class = case_when(
    str_detect(cell_type, "Inh") ~ "Inhibitory",
    str_detect(cell_type, "Exc") ~ "Excitatory",
    TRUE ~ "Non-Neuronal")) %>%
   # clean up the names in the term column
  mutate(term = recode(term,
                       `(Intercept)` = "Intercept",
                       `reportedGendermale` = "gender:Male",
                       `primaryDiagnosisAutism Spectrum Disorder` = "ASD",
                       `primaryDiagnosisControl` = "C",
                       `scale(ageDeath)` = "Age",
                       `scale(RIN)` = "RIN"))
 
```
Generating figures:
```{r}
#Cell-type estimations in BA9 and BA41

beta_plot_UCLA_ASD = combined_lms %>% 
  filter(term %in% 'ASD') %>% 
  mutate(cell_type = fct_reorder(cell_type, estimate)) %>% 
  ggplot(aes(x = cell_type, y = estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ggtitle("UCLA: Autism Spectrum Disorder vs. Controls") +
  ylab('Beta Coefficient') + 
  xlab('Cell Type Proportions') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(rows = vars(Brain_region), cols = vars(class), drop = T, scale = "free")
beta_plot_UCLA_ASD
```

```{r}
#Plotting the distribution of SST cells across all the primary diagnosis details given in the metadata
UCLA_estimations_metadata %>% ggplot(aes(x = primaryDiagnosisDetails, y = Inh_SST)) + geom_boxplot() + geom_jitter() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

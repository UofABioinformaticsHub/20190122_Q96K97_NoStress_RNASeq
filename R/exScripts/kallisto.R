library(ngsReports)
library(tidyverse)
library()
 
fdl <- list.files("/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/0_rawData/FastQC/", pattern = "zip", full.names = TRUE) %>%
  FastqcDataList()

plotGcContent(fdl, plotType = "line", theoreticalGC = TRUE, species = "Drerio", gcType = "Transc", usePlotly = TRUE)

counts <- list.dirs("/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/3_kallisto/") %>% 
  .[-1] %>%
  catchKallisto()

genesGR <- genes(ensDb)
trans <- transcripts(ensDb)

minCpm <- 1
minSamples <- 4
transDge <- counts$counts %>%
  set_colnames(basename(colnames(.))) %>%
  magrittr::extract(rowSums(cpm(.) > minCpm) >= minSamples,) %>%
  divide_by(counts$annotation[rownames(.),"Overdispersion"]) %>% # This is the key step that Gordon showed at a talk last year
  set_rownames(str_remove(rownames(.), "\\.[0-9]+")) %>%
  DGEList(genes = trans[rownames(.)])
transVoom <- transDge %>%
  calcNormFactors() %>%
  voomWithQualityWeights(plot = TRUE, design = matrix(1, nrow = ncol(.), ncol = 1))


geneDge <- transDge$counts %>%
  as.data.frame() %>%
  rownames_to_column("tx_id") %>%
  as_tibble() %>%
  left_join(transDge$genes) %>%
  dplyr::select(tx_id, gene_id, starts_with("Q"), starts_with("W"), -width) %>%
  gather(key = "sample", value = "count", -ends_with("id")) %>%
  group_by(gene_id, sample) %>%
  summarise(count = sum(count)) %>%
  spread(key = sample, value = count) %>%
  as.data.frame() %>%
  column_to_rownames("gene_id") %>%
  DGEList(genes = genesGR[rownames(.)]) %>%
  calcNormFactors()



pca <- counts$counts %>%
  set_colnames(basename(colnames(.))) %>%
  .[rowSums(cpm(.) > 1) > 4,] %>%
  cpm(log = TRUE) %>%
  t() %>%
  prcomp()

plotly::ggplotly(
  pca$x %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    as_tibble() %>%
    mutate(WT = grepl("W", sample)) %>%
    ggplot(aes(PC1, PC2, colour = WT, label = sample)) +
    geom_point(size = 3) +
    theme_bw() 
)

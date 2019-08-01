library(rtracklayer)
library(limma)
library(edgeR)
library(AnnotationHub)
library(tidyverse)
library(magrittr)
library(scales)
library(pander)
library(ggrepel)
library(here)

#Load files in with catchKallisto
counts <- list.files(
  path = "/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/3_kallisto/",
  full.names = TRUE) %>%
  catchKallisto

#Rename colnames
colnames(counts$counts) <- colnames(counts$counts) %>%
  basename() %>%
  paste("test", sep = "_")

#Import ensembl 96 gtf
gr <- import.gff("ftp://ftp.ensembl.org/pub/release-96/gtf/danio_rerio/Danio_rerio.GRCz11.96.gtf.gz")

#Generating gene names linked to transcript names
tr2gn <- gr %>% 
  subset(!is.na(transcript_id)) %>%
  subset(type == "transcript") %>%
  mcols() %>%
  as.data.frame() %>%
  dplyr::select(gene_id, transcript_id) %>%
  as_tibble()

#Create DGE list
#geneDGE <- 

counts$counts %>%
  as.data.frame() %>%
  rownames_to_column("transcript_id") %>%
  as_tibble() %>%
  mutate(transcript_id = str_remove_all(transcript_id, "\\.[0-9]+")) %>%
  mutate(transcript_id = str_remove_all(transcript_id, "_unspliced")) #%>%
  gather(key = "Sample", value = "Counts", ends_with("test")) %>%
  left_join(tr2gn) %>%
  group_by(Sample, gene_id) %>%
  summarise(Counts = sum(Counts)) %>%
  spread(key = "Sample", value = "Counts") %>%
  dplyr::filter(!is.na(gene_id)) %>%
  column_to_rownames("gene_id") %>%
  DGEList() %>%
  calcNormFactors()

#geneDGE <- counts$counts %>%
#  as.data.frame() %>%
#  rownames_to_column("transcript_id") %>%
#  dplyr::filter(!grepl("unspliced", transcript_id)) %>%
# as_tibble() %>%
#  mutate(transcript_id = str_remove_all(transcript_id, "\\.[0-9]+")) %>%
#  gather(key = "Sample", value = "Counts", ends_with("_test")) %>%
#  left_join(tr2gn) %>%
#  group_by(Sample, gene_id) %>%
#  summarise(Counts = sum(Counts)) %>%
#  spread(key = "Sample", value = "Counts") %>%
#  column_to_rownames("gene_id") #%>%
#  DGEList() %>%
#  calcNormFactors()
  
#Rename columns and rows
colnames(geneDGE$counts) <- str_remove(colnames(geneDGE$counts), "_test")
rownames(geneDGE$samples) <- str_remove(rownames(geneDGE$samples), "_test")

#Set group variable
geneDGE$samples$group <- colnames(geneDGE) %>%
  str_extract("(W|Q)") %>%
  factor(levels = c("W", "Q"))

#Add gene information
ah <- AnnotationHub()

#Subset Annotation Hub to search for zebrafish
ah %>%
  subset(species == "Danio rerio") %>%
  subset(dataprovider == "Ensembl") %>%
  subset(rdataclass == "EnsDb")

#Select correct Ensembl release
ensDb <- ah[["AH64906"]]

#Extract GenomicRanges object from ensDb
genesGR <- genes(ensDb)

#Remove redundant columns from mcols
mcols(genesGR) <- mcols(genesGR)[c("gene_id", "gene_name", "gene_biotype", "entrezid")]

#Use rownames of dgeList to reorder the genes object
geneDGE$genes <- genesGR[rownames(geneDGE),]

#Perform logical test to see how many genes were not detected in dataset
geneDGE$counts %>%
  rowSums() %>%
  is_greater_than(0) %>%
  table()

#Check for genes having > 4 samples with cpm > 1
geneDGE %>%
  cpm() %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_weakly_greater_than(4) %>%
  table()

#Create logical vector of genes to keep for further analysis
genes2keep <- geneDGE %>%
  cpm() %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_weakly_greater_than(4)

#Create new DGE list of filtered genes
dgeFilt <- geneDGE[genes2keep,, keep.lib.sizes = FALSE] %>%
  calcNormFactors()

#Compare distributions of the two datasets
par(mfrow = c(1,2))
geneDGE %>%
  cpm(log = TRUE) %>%
  plotDensities(legend = FALSE, main = "Before Filtering")
dgeFilt %>%
  cpm(log = TRUE) %>%
  plotDensities(legend = FALSE, main = "After Filtering")
par(mfrow = c(1,1))

#Create box plot to check library sizes
dgeFilt$samples %>%
  ggplot(aes(group, lib.size, fill = group)) +
  geom_boxplot() +
  scale_y_continuous(labels = comma) +
  labs(x = "Genotype", y = "Library Size") +
  scale_fill_discrete(name ="Genotype", labels = c("Wildtype", "Mutant")) +
  scale_x_discrete(labels=c("W" = "Wildtype", "Q" = "Mutant")) +
  theme_bw()

#Assess cpm values to make sure PCA results are not heavily skewed by highly expressed genes
pca <- dgeFilt %>%
  cpm(log = TRUE) %>%
  t() %>%
  prcomp()

#Quick inspection to check whether first two PCA components capture most of the variability
summary(pca)$importance %>% pander(split.tables = Inf)

design <- model.matrix(~group, data = dgeFilt$samples)

#Create object for voom analysis
voomData <- voomWithQualityWeights(dgeFilt, design = design, plot = TRUE)

#Plot PCA
#plotly::ggplotly(
pca$x %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  as_tibble() %>%
  dplyr::select(sample, PC1, PC2) %>%
  left_join(rownames_to_column(voomData$targets, "sample")) %>%
  ggplot(aes(PC1, PC2, colour = group, label = sample)) +
  geom_point() +
  geom_text_repel() +
  theme_bw()
#  )
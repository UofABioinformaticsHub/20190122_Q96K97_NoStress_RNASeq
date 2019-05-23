#Load libraries
library(limma)
library(edgeR)
library(AnnotationHub)
library(tidyverse)
library(magrittr)
library(scales)
library(pander)
library(ggrepel)
library(here)

#Load counts
counts <- read_tsv(here("../2_alignedData/featureCounts/genes.out")) %>%
  set_colnames(basename(colnames(.))) %>%
  set_colnames(str_remove(colnames(.), "Aligned.sortedByCoord.out.bam"))  %>%
  dplyr::select(Geneid, starts_with("W"), starts_with("Q"))

#Create DGE list and calculate normalisaton factors
dgeList <- counts %>%
  as.data.frame() %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors()

#Set group variable
#Why do we need to do this?
dgeList$samples$group <- colnames(dgeList) %>%
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
#How do we know which release matches the the data we have generated?
ensDb <- ah[["AH64906"]]

#Extract GenomicRanges object from ensDb
genesGR <- genes(ensDb)

#The following can be used to return information about the GRanges object:
#granges(genesGR)
#seqinfo(genesGR)
#mcols(genesGR)

#Remove redundant columns from mcols
mcols(genesGR) <- mcols(genesGR)[c("gene_id", "gene_name", "gene_biotype", "entrezid")]

#Use rownames of dgeList to reorder the genes object
dgeList$genes <- genesGR[rownames(dgeList),]

#Perform logical test to see how many genes were not detected in dataset
dgeList$counts %>%
  rowSums() %>%
  is_greater_than(0) %>%
  table()

#Check for genes having > 4 samples with cpm > 1
#Why use > smallest group and not >=? Is this summing the cpm's or the number of samples > 1?
dgeList %>%
  cpm() %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_weakly_greater_than(4) %>%
  table()
#Approximately 40% of genes do not fit this criteria

#Create logical vector of genes to keep for further analysis
genes2keep <- dgeList %>%
  cpm() %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_weakly_greater_than(4)

#Create new DGE list of filtered genes
dgeFilt <- dgeList[genes2keep,, keep.lib.sizes = FALSE] %>%
  calcNormFactors()

#Compare distributions of the two datasets
par(mfrow = c(1,2))
dgeList %>%
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

#What does groupQ do?
topTable <- voomData %>% 
  lmFit(design = design) %>%
  eBayes %>%
  topTable(coef = "groupQ", n = Inf) %>%
  as_tibble() %>%
  unite("Range", ID.start, ID.end, sep = "-") %>%
  unite("Location", ID.seqnames, Range, ID.strand, sep = ":") %>%
  dplyr::select(Geneid = ID.gene_id, 
                Symbol = ID.gene_name,
                AveExpr, logFC, t, P.Value, 
                FDR = adj.P.Val, 
                Location, 
                Entrez = ID.entrezid) %>%
  mutate(DE = FDR < 0.05)

topTable %>%
  mutate(DE = FDR < 0.05) %>%
  ggplot(aes(logFC, -log10(P.Value), colour = DE)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = . %>% 
                    dplyr::filter(DE) %>%
                    dplyr::filter(-log10(P.Value) > 4 | abs(logFC) > 2.5),
                  aes(label = Symbol)) + 
  scale_colour_manual(values = c("grey", "red")) +
  theme_bw() +
  theme(legend.position = "none")

  topTable %>%
  mutate(DE = FDR < 0.05) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(AveExpr, logFC, colour = DE)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = . %>% 
                    dplyr::filter(DE) %>%
                    dplyr::filter(abs(logFC) > 2 | AveExpr > 14),
                  aes(label = Symbol)) + 
  scale_colour_manual(values = c("grey", "red")) +
  labs(x = "Average Expression (log2 CPM)",
       y = "log Fold-Change") +
  theme_bw() +
  theme(legend.position = "none")
  
topTable %>%
   dplyr::filter(FDR < 0.05, abs(logFC) > 1) %>%
   dplyr::select(ID = Geneid, Symbol, AveExpr, logFC, P.Value, FDR) %>%
   pander(caption = paste("The", nrow(.), "most DE genes when ranked by p-value, and filtered on a logFC beyond the range $\\pm1$"))
  
ens2Entrez <- file.path("https://uofabioinformaticshub.github.io/Intro-NGS-fib", "data", "ens2Entrez.tsv") %>% 
  url() %>%
  read_tsv()

de <- topTable %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Geneid) %>%
  left_join(ens2Entrez) %>%
  dplyr::filter(!is.na(Entrez)) %>%
  .[["Entrez"]] %>%
  unique()
uv <- topTable %>%
  dplyr::select(Geneid) %>%
  left_join(ens2Entrez) %>%
  dplyr::filter(!is.na(Entrez)) %>%
  .[["Entrez"]] %>%
  unique()

goResults <- goana(de = de, universe = uv, species = "Hs")

goResults %>% 
  rownames_to_column("GO ID") %>%
  as_tibble() %>%
  dplyr::filter(DE > 1) %>%
  arrange(P.DE) %>%
  mutate(FDR = p.adjust(P.DE, "fdr")) %>%
  dplyr::filter(FDR < 0.05) %>%
  mutate(`GO ID` = str_replace(`GO ID`, ":", "\\\\:")) %>%
  pander(caption = "GO Terms potentially enriched in the set of differentially expressed genes")

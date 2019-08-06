### 1. 2019 Analysis
  - [x] DE analysis
  - [x] fgsea
  - [ ] fry
    - [ ] Compare p-values from fgsea and fry (using -log10) just for fun
  - [x] Network visualistion of GSEA results
    - [ ] Select best visualisation and zoom in on OXPHOS (`clusterProfiler::cnetplot()`)
  - [ ] IRE
  - Put aside for now:
    - Community detection using Ville's turnip & burdock
      - Only use significant pathways
      - Only use genes up to the leading edge
      - Find common genes between pathways that c luster together

### 2. Integration with 2017
  - [x] Summary of 2017
    - [x] DGE Analysis
    - [x] GSEA Analysis
  - [x] Meta analysis of GSEA/fry results
    - Use -log10(p) & Fisher's Method
    - No easy visualisation...
  - Integration of counts
    - [x] LogFC comparison
    - [x] Put counts in same dge list
    - [x] MDS/PCA
    - [x] RUVSeq
      - [ ] Compare with Dream results (logFC vs logFC; -log10(p) vs -log10(p))
      - [ ] fgsea/fry on RUVSeq results
    - [ ] Dream
      - [ ] DGE
      - [ ] fgsea/fry on Dream results
      
### 3. Analysis of Wildtypes
  - Compare WT's at 6 months between datasets
    - [ ] Plot WTs on PCA and find genes which are biggest contributors
    - [ ] WGCNA to identify modules clearly identified between families
      - 20181113_MorganLardelli_mRNASeq: psen2 (samples w/o "FAD" or "FS")
      - 2017_Lardelli_Q96K97del ("non-mutant")
      - 2017_Lardelli_K97Gfs ("non-mutant")
      - 20170906_Morgan_Hypoxia ("6P_PN")
      - 170327MichaelLardelli_NextSeq ("WT_6month")
      - Karissa's
    - Glucocorticoid receptor target genes?
  
### 4. Integration using alternative approaches
  - [ ] Singscore
    - Initially within single datasets
  - [ ] Concordance (Blalock) (Unlikely enough time)
  
### 5. Analysis using other mutants (Psen2, Sorl1, Q96K97)
  - [ ] WGCNA or Singscore
    - EOfAD genotypes: Sorl1l- ; Psen2^S4Ter^; Psen2^T141_L142delinsMISLISV^; psen1^Q96K97^
    - non-EOfAD genotypes: psen1^K97fs^; psen2^N140fs^
  - [ ] Meta-analysis using Fisher's method & Stouffer's method (use all genesets)
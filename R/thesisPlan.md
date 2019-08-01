### 1. 2019 Analysis
  - [x] DE analysis
  - [x] GSEA
    - [x] Wiki pathways
    - [x] 50 Hallmark
    - [x] KEGG
  - [ ] IRE
  - [ ] Compare p-values from GSEA and fry (using -log10) just for fun
  - [ ] Select best visualisation and zoom in on OXPHOS (`clusterProfiler::cnetplot()`)
  - Put aside for now:
    - Community detection using Ville's turnip & burdock
      - Only use significant pathways
      - Only use genes up to the leading edge
      - Find common genes between pathways that cluster together

### 2. Integration with 2017
  - [x] Summary of 2017
    - [x] DGE Analysis
    - [x] GSEA Analysis
  - [x] Meta analysis of GSEA/fry results
    - Use -log10(p) & Fisher's Method
    - No easy visualisation...
  - [ ] Integration of counts
    - [x] LogFC comparison
    - [x] Put counts in same dge list
    - [x] MDS/PCA
    - [x] RUVSeq
      - [ ] Compare with Dream results (logFC vs logFC; -log10(p) vs -log10(p))
      - [ ] GSEA (or fry)    
    - [ ] Dream
      - [ ] DGE
      - [ ] GSEA/fry on Dream Results

### 3. Integration using alternative approaches
  - [ ] Singscore
  - [ ] Concordance (Blalock)
  
### 4. Analysis using other mutants (Psen2, Sorl1, Q96K97)
    - Do we include non-EOfAD's? (K97fs, Psen null)
    - Meta-analysis using Fisher's method & Stouffer's method (use all genesets)
  - [ ] WGCNA or Singscore
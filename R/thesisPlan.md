1. 2019
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

2. Integration with 2017
  1. [x] Summary of 2017
    1. [x] DGE Analysis
    2. [x] GSEA Analysis
  2. [x] Meta analysis of GSEA/fry results
    - Use -log10(p) & Fisher's Method
    - No easy visualisation...
  3. Integration of counts
    1. [x] LogFC comparison
    2. [x] Put counts in same dge list
    3. [x] MDS/PCA
    4. RUVSeq
      1. Compare with Dream results (logFC vs logFC; -log10(p) vs -log10(p))
      2. GSEA (or fry)    
    5. Dream
      1. DGE
      2. GSEA/fry on Dream Results

3. Integration using alternative approaches
  1. Singscore
  2. Concordance (Blalock)
  
4. Analysis using other mutants (Psen2, Sorl1, Q96K97)
    - Do we include non-EOfAD's? (K97fs, Psen null)
    - Meta-analysis using Fisher's method & Stouffer's method (use all genesets)
  1. WGCNA or Singscore
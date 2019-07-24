1. 2019
  - [x] DE analysis
  - [x] GSEA
    - [x] Wiki pathways
    - [x] 50 Hallmark
    - [x] KEGG
    - IRE
  - Compare p-values from GSEA and fry (using -log10) just for fun
  - Select best visualisation and zoom in on OXPHOS (`clusterProfiler::cnetplot()`)
  x Community detection using Ville's turnip & burdock
    x Only use significant pathways
    x Only use genes up to the leading edge
    x Find common genes between pathways that cluster together

2. Integration with 2017
  - Summary of 2017
    - [x] DGE Analysis
    - [x] GSEA Analysis
  - Meta analysis of GSEA/fry results
    - Use -log10(p) & Fisher's Method
    - No easy visualisation...
  - Integration of counts
    1. Put counts in same dge list
    2. MDS/PCA
    3. Dream
      - DGE
      - GSEA/fry on Dream Results
    4. RUVSeq
      - Compare with Dream results (logFC vs logFC; -log10(p) vs -log10(p))
      - GSEA (or fry)

  
3. Integration using alternative approaches
  - Singscore
  - Concordance (Blalock)
  
4. Analysis using other mutants (Psen2, Sorl1, Q96K97)
    - Do we include non-EOfAD's? (K97fs, Psen null)
    - Meta-analysis using Fisher's method & Stouffer's method (use all genesets)
  - WGCNA or Singscore
  
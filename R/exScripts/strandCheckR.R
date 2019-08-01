library(strandCheckR)

files <- list.files(
  "/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/2_alignedData/bam/",pattern = ".+\\.bam$",full.names = TRUE
)

win <- getWinFromBamFile(files[1], sequences = "9")
# shorten the file name
win$File <- basename(as.character(win$File))
win

plotHist(
  windows = win, group_by = c("File","OverlapTranscript"), 
  normalize_by = "File", scales = "free_y"
)

filterDNA(
  file = "/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/2_alignedData/bam/Q1Aligned.sortedByCoord.out.bam", destination = "/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/2_alignedData/bam_filtered/Q1Aligned.sortedByCoord.out.bam", 
  statfile = "statfile", 
  threshold = 0.8)

postwin <- getWinFromBamFile("2_alignedData/bam_filtered/Q1Aligned.sortedByCoord.out.bam", sequences = 9)

plotHist(
  windows = postwin, group_by = c("File","OverlapTranscript"), 
  normalize_by = "File", scales = "free_y"
)


katPoint <- function(data, sample = "sample", min.mut = 6, max.dis = 1000, 
                     txdb = NULL) {
  build = data$build[1]
  genome.opts = c("mm10", "mm9")
  if (!build %in% genome.opts) {
    stop("Available reference builds: mm10, mm9")
  }
  if (build == "mm10") {
    chr.arm = c(1.957e+08, 1.821e+08, 1.6e+08, 1.565e+08, 1.518e+08, 
                1.497e+08, 1.454e+08, 1.294e+08, 1.246e+08, 1.307e+08, 
                1.221e+08, 1.201e+08, 1.204e+08, 1.249e+08, 1.04e+08, 
                9.82e+07, 9.498e+07, 9.071e+07, 6.143e+07, 1.71e+08, 
                9.174e+07)
  } else if (build == "mm9") {
    chr.arm = c(1.972e+08, 1.817e+08, 1.596e+08, 1.556e+08, 1.525e+08, 
                1.495e+08, 1.525e+08, 1.317e+08, 1.24e+08, 1.3e+08, 
                1.218e+08, 1.213e+08, 1.202e+08, 1.252e+08, 1.034e+08, 
                9.832e+07, 9.527e+07, 9.077e+07, 6.134e+07, 1.666e+08, 
                9.174e+07)
  } else {
    stop("Available reference builds: mm10, mm9")
  }
  num = dim(data)[1] - 5
  katPoint <- matrix(nrow = num, ncol = 8)
  i = 1
  mutnum = 1
  Cmutnum = 0
  for (i in 1:num) {
    if (data$ref[i] %in% c("C", "G")){
      Cmutnum = Cmutnum + 1
    }
    if (data$dis[i + 1] <= max.dis) {
      mutnum = mutnum + 1
    } else {
      if (mutnum >= min.mut) {
        len = data$pos[i] - data$pos[i - mutnum + 1] + 1
        chr.n = gsub(pattern = "chr", replacement = "", x = data$chr[i], 
                     fixed = TRUE)
        chr.n = gsub(pattern = "X", replacement = "20", x = chr.n, 
                     fixed = TRUE)
        chr.n = gsub(pattern = "Y", replacement = "21", x = chr.n, 
                     fixed = TRUE)
        chr.n = as.numeric(chr.n)
        if (data$pos[i] <= chr.arm[chr.n]) {
          arm = paste(chr.n, "p", sep = "")
        } else if (data$pos[i - mutnum + 1] >= chr.arm[chr.n]) {
          arm = paste(chr.n, "q", sep = "")
        } else {
          arm = paste(chr.n, "p, ", chr.n, "q", sep = "")
        }
        katPoint[i, 1:8] = c(sample, data$chr[i], data$pos[i - mutnum + 
                                                             1], data$pos[i], arm, len, mutnum, round(Cmutnum/mutnum,3))
      }
      mutnum = 1
      Cmutnum = 0
    }
  }
  katPoint.out = data.frame(na.omit(katPoint))
  names(katPoint.out) = c("sample", "chrom", "start", "end", "chrom.arm", "length", "number.mut", 
                          "weight.C>X")
  for (i in 1:dim(katPoint.out)[1]) {
    if (as.numeric(as.character(katPoint.out$"weight.C>X"[i])) < 0.8) {
      katPoint.out$confidence[i] = 0
    } else {
      katPoint.out$confidence[i] <- length(which(subset(katPoint.out,
                                                        as.numeric(as.character(katPoint.out$"weight.C>X")) >= 0.8)$chrom == katPoint.out$chrom[i]))
      if (katPoint.out$confidence[i] > 3) {
        katPoint.out$confidence[i] = 3
      }
    }
  }
  if (!is.null(txdb)) {
    gr <- GRanges(seqnames = Rle(katPoint.out$chrom), ranges=IRanges(start = 
                                                                       as.numeric(as.character(katPoint.out$start)), end =as.numeric(as.character(katPoint.out$end))))
    peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
    katPoint.out$annotation <- peakAnno@anno$annotation
    katPoint.out$distanceToTSS <- peakAnno@anno$distanceToTSS
    katPoint.out$geneName <- peakAnno@anno$SYMBOL
    katPoint.out$geneID <- peakAnno@anno$geneId
  } 
  message(paste(dim(katPoint.out)[1], "potential kataegis events identified", 
                sep = " "))
  return(katPoint.out)
}

mutDis.plot <- function(plot.data, sample = "sample", chr = NULL, color = NULL, 
                        min.mut = 6, max.dis = 1000) {
  if (is.null(color)) {
    col = RColorBrewer::brewer.pal(n = 6, name = "Set1")
  } else {
    col = color
  }
  names(col) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  build = plot.data$build[1]
  genome.opts = c("mm10", "mm9")
  if (!build %in% genome.opts) {
    stop("Available reference builds: mm10, mm9")
  }
  if (build == "mm10") {
    chr.lens = c(195471971, 182113224, 160039680, 156508116, 151834684, 
                 149736546, 145441459, 129401213, 124595110, 130694993, 
                 122082543, 120129022, 120421639, 124902244, 104043685, 
                 98207768, 94987271, 90702639, 61431566, 171031299, 91744698)
    chr.arm = c(1.97e+08, 1.82e+08, 1.6e+08, 1.56e+08, 1.51e+08, 
                1.49e+08, 1.45e+08, 1.29e+08, 1.24e+08, 1.3e+08, 
                1.22e+08, 1.2e+08, 1.2e+08, 1.25e+08, 1.04e+08, 
                9.8e+07, 9.5e+07, 9.1e+07, 6.1e+07, 1.71e+08, 9.17e+07)
  } else if (build == "mm9") {
    chr.lens = c(197195432, 181748087, 159599783, 155630120, 152537259, 
                 149517037, 152524553, 131738871, 124076172, 129993255, 
                 121843856, 121257530, 120284312, 125194864, 103494974, 
                 98319150, 95272651, 90772031, 61342430, 166650296, 91744698)
    chr.arm = c(1.97e+08, 1.82e+08, 1.6e+08, 1.56e+08, 1.51e+08, 
                1.49e+08, 1.52e+08, 1.31e+08, 1.24e+08, 1.3e+08, 
                1.22e+08, 1.21e+08, 1.2e+08, 1.25e+08, 1.03e+08, 
                9.8e+07, 9.5e+07, 9.1e+07, 6.1e+07, 1.67e+08, 9.17e+07)
  } else {
    stop("Available reference builds: mm10, mm9")
  }
  if (is.null(chr)) {
    seq = c(1:21)
    seq0 = c(1:19, "X", "Y")
  } else {
    plot.data$pos.updated = plot.data$pos
    seq0 = gsub(pattern = "chr", replacement = "", x = chr, fixed = TRUE)
    seq = gsub(pattern = "X", replacement = "20", x = seq0, fixed = TRUE)
    seq = gsub(pattern = "Y", replacement = "21", x = seq, fixed = TRUE)
    seq = as.numeric(seq)
    plot.data = plot.data[which(plot.data$seq %in% seq), ]
    for (i in 1:length(chr.lens)) {
      if (!i %in% seq) {
        chr.lens[i] = 0
      }
    }
  }
  xlim = c(0, sum(chr.lens))
  chr.lens.sum = cumsum(chr.lens)
  chr.lens.sum = c(0, chr.lens.sum)
  plot.data$pos.updated = plot.data$pos + chr.lens.sum[plot.data$seq]
  chr.abline = chr.lens[seq]
  chr.abline.sum = cumsum(chr.abline)
  chr.abline.sum = c(0, chr.abline.sum)
  chr.axis = seq
  for (i in 1:length(seq)) {
    chr.axis[i] = (chr.abline.sum[i] + chr.abline.sum[i + 1])/2
  }
  par(mai = c(1, 1, 1, 1.5))
  plot.data$dis.log10 = log10(plot.data$dis)
  plot(plot.data$pos.updated, plot.data$dis.log10, col = col[plot.data$alteration], 
       pch = 16, cex = 0.5, xlim = xlim, axes = F, box(), mgp = c(2, 
                                                                  1, 0), main = paste("Rainfall plot for", sample, sep = " "), 
       xlab = "Chromosomes", ylab = "Intermutation distance (bp, log10)")
  abline(h = 2, col = "red", lty = 2)
  abline(v = chr.abline.sum, col = "grey")
  axis(2, las = 1, lwd.tick = 0.5)
  axis(1, at = chr.axis, labels = c(seq0), tick = FALSE, mgp = c(1, 
                                                                 0.5, 0))
  legend("right", names(col), col = col, pch = 16, inset = c(-0.15, 
                                                             0), bty = "n", xpd = TRUE)
  arrows.point = c()
  
  num = dim(plot.data)[1] - 5
  katPoint <- matrix(nrow = num, ncol = 8)
  i = 1
  mutnum = 1
  Cmutnum = 0
  for (i in 1:num) {
    if (plot.data$ref[i] %in% c("C", "G")){
      Cmutnum = Cmutnum + 1
    }
    if (plot.data$dis[i + 1] <= max.dis) {
      mutnum = mutnum + 1
    } else {
      if (mutnum >= min.mut) {
        len = plot.data$pos[i] - plot.data$pos[i - mutnum + 
                                                 1] + 1
        chr.n = gsub(pattern = "chr", replacement = "", x = plot.data$chr[i], 
                     fixed = TRUE)
        chr.n = gsub(pattern = "X", replacement = "20", x = chr.n, 
                     fixed = TRUE)
        chr.n = gsub(pattern = "Y", replacement = "21", x = chr.n, 
                     fixed = TRUE)
        chr.n = as.numeric(chr.n)
        if (plot.data$pos[i] <= chr.arm[chr.n]) {
          arm = paste(chr.n, "p", sep = "")
        } else if (plot.data$pos[i - mutnum + 1] >= chr.arm[chr.n]) {
          arm = paste(chr.n, "q", sep = "")
        } else {
          arm = paste(chr.n, "p, ", chr.n, "q", sep = "")
        }
        wei = round(Cmutnum/mutnum,3)
        katPoint[i, 1:8] = c(sample, plot.data$chr[i], plot.data$pos[i - 
                                                                       mutnum + 1], plot.data$pos[i], arm, len, mutnum, wei)
        arrows.point = (plot.data$pos.updated[i] + plot.data$pos.updated[i - 
                                                                           mutnum + 1])/2
        if (wei >= 0.8) {
          arrows(arrows.point, 0, arrows.point, 0.5, length = 0.1, 
                 angle = 30, code = 2, col="black")
        } 
      }
      mutnum = 1
      Cmutnum = 0
    }
  }
  message(paste(dim(na.omit(katPoint))[1], "potential kataegis events identified. Use katPoint() function to see detail.", 
                sep = " "))
}

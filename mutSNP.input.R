mutSNP.input <- function(mut.data, chr = "chr", pos = "pos", ref = "ref", 
                         alt = "alt", build = NULL, k = 10) {
  if (exists("mut.data", mode = "list")) {
    mut.data <- mut.data
  } else {
    if (file.exists(mut.data)) {
      mut.data <- utils::read.table(mut.data, sep = "\t", header = TRUE, 
                                    as.is = FALSE, check.names = FALSE)
    } else {
      stop("mut.data is neither a file nor a loaded data frame")
    }
  }
  mut.data <- mut.data[, c(chr, pos, ref, alt)]
  genome.opts = c("mm10", "mm9")
  if (!build %in% genome.opts || is.null(build)) {
    stop("Available reference builds: mm10, mm9")
  }
  if (build == "mm10") {
    chr.lens = c(195471971, 182113224, 160039680, 156508116, 151834684, 
                 149736546, 145441459, 129401213, 124595110, 130694993, 
                 122082543, 120129022, 120421639, 124902244, 104043685, 
                 98207768, 94987271, 90702639, 61431566, 171031299, 91744698)
    bsg = BSgenome.Mmusculus.UCSC.mm10
  } else if (build == "mm9") {
    chr.lens = c(197195432, 181748087, 159599783, 155630120, 152537259, 
                 149517037, 152524553, 131738871, 124076172, 129993255, 
                 121843856, 121257530, 120284312, 125194864, 103494974, 
                 98319150, 95272651, 90772031, 61342430, 166650296, 91744698)
    bsg = BSgenome.Mmusculus.UCSC.mm9
  } else {
    stop("Available reference builds: mm10, mm9")
  }
  mut.data$build = build
  if (!all(mut.data$ref %in% DNA_BASES & mut.data$alt %in% DNA_BASES)) {
    stop("Only SNV substitutions are currently supported.")
  }
  ref_base = DNAStringSet(mut.data$ref)
  alt_base = DNAStringSet(mut.data$alt)
  # Ensure start and end are within chromosome bounds
  conv.start = pmax(mut.data$pos - k, 1)  # Ensure start is not less than 1
  conv.end = pmin(mut.data$pos + k, sapply(mut.data$chr, function(chr) seqlengths(bsg)[chr]))  # Ensure end is within the chromosome length
  
  # Filter out any invalid ranges where start > end
  valid_ranges <- conv.start <= conv.end
  
  if (any(!valid_ranges)) {
    warning(paste("Skipping", sum(!valid_ranges), "invalid ranges where start > end"))
  }
  
  ref_base = ref_base[valid_ranges]
  alt_base = alt_base[valid_ranges]
  conv.start = conv.start[valid_ranges]
  conv.end = conv.end[valid_ranges]
  mut.data = mut.data[valid_ranges, ]
  
  # Create a GRanges object to retrieve sequences
  gr <- GRanges(seqnames = mut.data$chr, ranges = IRanges(start = conv.start, end = conv.end))
  
  context = getSeq(bsg, gr)
  if (TRUE) {
    idx = mut.data$ref %in% c("A", "G")
    context[idx] = reverseComplement(context[idx])
    ref_base[idx] = reverseComplement(ref_base[idx])
    alt_base[idx] = reverseComplement(alt_base[idx])
  }
  mut.data$alteration = paste(ref_base, alt_base, sep = ">")
  mut.data$context = context
  # Replace chr X and Y with numeric value (20 and 21) for better ordering
  seq = gsub(pattern = "chr", replacement = "", x = mut.data$chr, fixed = TRUE)
  seq = gsub(pattern = "X", replacement = "20", x = seq, fixed = TRUE)
  seq = gsub(pattern = "Y", replacement = "21", x = seq, fixed = TRUE)
  mut.data$seq = as.numeric(seq)
  mut.data = mut.data[order(mut.data$seq, mut.data$pos), ]
  chr.lens.sum = cumsum(chr.lens)
  chr.lens.sum = c(0, chr.lens.sum)
  mut.data$dis = c(mut.data$pos[1], diff(mut.data$pos + chr.lens.sum[mut.data$seq]))
  return(mut.data)
}

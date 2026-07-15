
read_plink_out <- function(path, ...) {
  if (!file.exists(path)) {
    stop(sprintf("Expected PLINK output file not found: %s", path))
  }
  read.table(path, ...)
}

save_hist <- function(x, file, title, xlab) {
  jpeg(file, 900, 650)
  hist(x, col = "steelblue", breaks = 40, main = title, xlab = xlab)
  dev.off()
}

frq   <- read_plink_out("outputs/SNPfreq.frq", header = TRUE)

cat("Minimum MAF :", min(frq$MAF), "\n")
cat("Maximum MAF :", max(frq$MAF), "\n")


save_hist(frq$MAF,       "outputs/MAF_Histogram.jpg",      "Minor Allele Frequency", "MAF")
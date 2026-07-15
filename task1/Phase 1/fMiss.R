
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

lmiss <- read_plink_out("outputs/miss.lmiss", header = TRUE)
imiss <- read_plink_out("outputs/miss.imiss", header = TRUE)

save_hist(lmiss$F_MISS,  "outputs/SNP_Missingness.jpg",    "SNP Missingness",        "Missing Rate")
save_hist(imiss$F_MISS,  "outputs/Sample_Missingness.jpg", "Sample Missingness",     "Missing Rate")
run_plink <- function(cmd, step_name = NULL) {
  full_cmd <- paste("plink1.9", cmd)
  status <- system(full_cmd)
  if (status != 0) {
    stop(sprintf(
      "PLINK step failed%s (exit status %d).\nCommand: %s",
      if (!is.null(step_name)) paste0(" [", step_name, "]") else "",
      status, full_cmd
    ))
  }
  invisible(status)
}

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

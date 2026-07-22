require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R package(s): ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

analysis_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#17324D"),
      plot.subtitle = ggplot2::element_text(colour = "#4E6476"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )
}

save_plot <- function(plot, filename, width = 8, height = 5, dpi = 320) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = "white")
}

write_session_info <- function(filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), filename)
}

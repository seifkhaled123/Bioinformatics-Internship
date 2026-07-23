source("professional_helpers.R")

source("functions.R")
require_packages(c("readxl", "dplyr", "tidyr", "ggplot2", "glmnet", "randomForest", "pROC", "caret", "igraph", "ggrepel"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

task2_dir <- normalizePath(".")
root <- normalizePath("..")
phase_path <- function(phase, ...) file.path(task2_dir, paste("Phase", phase), ...)
out_dir <- function(phase) {
  x <- phase_path(phase, "outputs")
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}
write_csv <- function(x, phase, name, row.names = FALSE) write.csv(x, file.path(out_dir(phase), name), row.names = row.names)

read_qc_data <- function(scaled = TRUE) {
  name <- if (scaled) "Qatari_metabolomics_Cleaned_Scaled.csv" else "Qatari_metabolomics_Cleaned_Imputed.csv"
  path <- phase_path(11, "outputs", name)
  if (!file.exists(path)) stop("Run Phase 11 first; missing ", path)
  read.csv(path, check.names = FALSE)
}

phase11 <- function() {
  raw <- as.data.frame(readxl::read_excel(file.path(task2_dir, "data", "Qatari_metabolomics.xlsx")))
  mapping <- as.data.frame(readxl::read_excel(file.path(task2_dir, "data", "mapping.xlsx")))
  stopifnot(!anyDuplicated(mapping$mapped_id), !anyDuplicated(mapping$main_id))
  data <- merge(mapping, raw, by = "mapped_id", all.x = TRUE, sort = FALSE)
  if (anyNA(data$Diabetes)) stop("Some mapped metabolomics records have no diabetes status.")
  features <- setdiff(names(raw), c("mapped_id", "Diabetes"))
  x <- data[, features, drop = FALSE]
  if (!all(vapply(x, is.numeric, logical(1)))) stop("All metabolite columns must be numeric.")
  rownames(x) <- data$main_id

  metabolite_missing <- colMeans(is.na(x))
  keep_features <- names(metabolite_missing)[metabolite_missing <= .20]
  sample_missing <- rowMeans(is.na(x[, keep_features, drop = FALSE]))
  keep_samples <- names(sample_missing)[sample_missing <= .20]
  x <- x[keep_samples, keep_features, drop = FALSE]
  variances <- vapply(x, var, numeric(1), na.rm = TRUE)
  keep_features <- names(variances)[is.finite(variances) & variances > 1e-6]
  x <- x[, keep_features, drop = FALSE]

  # Median imputation is performed only after the documented exclusion thresholds.
  imputed_cells <- colSums(is.na(x))
  for (met in names(x)) x[[met]][is.na(x[[met]])] <- median(x[[met]], na.rm = TRUE)
  x_imputed <- x
  x_scaled <- as.data.frame(scale(x_imputed))
  metadata <- data.frame(main_id = rownames(x), Diabetes = data$Diabetes[match(rownames(x), data$main_id)])
  imputed_output <- cbind(metadata, x_imputed)
  scaled_output <- cbind(metadata, x_scaled)

  write_csv(data.frame(Metabolite = names(metabolite_missing), Missing_Proportion = metabolite_missing, Passed_20pct = metabolite_missing <= .20), 11, "metabolite_missingness.csv")
  write_csv(data.frame(main_id = names(sample_missing), Missing_Proportion = sample_missing, Passed_20pct = sample_missing <= .20), 11, "sample_missingness.csv")
  write_csv(data.frame(Metabolite = names(variances), Variance = variances, Passed_variance_filter = names(variances) %in% keep_features, Imputed_values = imputed_cells[names(variances)]), 11, "metabolite_qc_metrics.csv")
  write.table(setdiff(features, keep_features), file.path(out_dir(11), "metabolites_removed.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(keep_features, file.path(out_dir(11), "metabolites_passed.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write_csv(imputed_output, 11, "Qatari_metabolomics_Cleaned_Imputed.csv")
  write_csv(scaled_output, 11, "Qatari_metabolomics_Cleaned_Scaled.csv")
  write_csv(data.frame(Metric = c("samples_before", "samples_after", "metabolites_before", "metabolites_after", "samples_removed", "metabolites_removed", "imputed_values"), Value = c(nrow(raw), nrow(x), length(features), ncol(x), nrow(raw) - nrow(x), length(features) - ncol(x), sum(imputed_cells))), 11, "qc_summary.csv")

  plot_data <- bind_rows(
    data.frame(kind = "Metabolites", missingness = metabolite_missing),
    data.frame(kind = "Samples", missingness = sample_missing)
  )
  p <- ggplot(plot_data, aes(missingness, fill = kind)) + geom_histogram(bins = 25, alpha = .75, position = "identity") + geom_vline(xintercept = .20, linetype = 2, colour = "#B22222") + facet_wrap(~kind, scales = "free_y") + labs(title = "Metabolomics missingness before QC", subtitle = "Dashed line: 20% exclusion threshold", x = "Missing-value proportion", y = "Count") + analysis_theme() + theme(legend.position = "none")
  save_plot(p, file.path(out_dir(11), "missingness_distribution.png"), 8, 4.5)
  write_session_info(file.path(out_dir(11), "sessionInfo.txt"))
}

load_analysis_data <- function() {
  metabolites <- read_qc_data(scaled = TRUE)
  covariates <- read.table(file.path(root, "task1", "Phase 4", "outputs", "covariates.txt"), header = TRUE, check.names = FALSE)
  dat <- merge(metabolites, covariates, by.x = "main_id", by.y = "IID", all.x = TRUE, sort = FALSE)
  pcs <- paste0("PC", 1:5)
  needed <- c("main_id", "Diabetes", "Sex", pcs)
  if (anyNA(dat[, needed])) stop("Missing matched covariates. Check Task 1 Phase 4 output and IDs.")
  dat$Sex <- factor(dat$Sex, levels = c(1, 2), labels = c("Male", "Female"))
  dat
}

phase12 <- function() {
  dat <- load_analysis_data()
  feature_names <- names(read_qc_data(scaled = TRUE))[-c(1, 2)]
  model_covariates <- c("Diabetes", "Sex", paste0("PC", 1:5))
  result <- lapply(feature_names, function(met) {
    fit <- lm(reformulate(model_covariates, response = met), data = dat)
    coef <- summary(fit)$coefficients["Diabetes", ]
    data.frame(Metabolite = met, N = stats::nobs(fit), Beta = coef[["Estimate"]], SE = coef[["Std. Error"]], CI_low = coef[["Estimate"]] - 1.96 * coef[["Std. Error"]], CI_high = coef[["Estimate"]] + 1.96 * coef[["Std. Error"]], P_value = coef[["Pr(>|t|)"]])
  }) |> bind_rows() |> mutate(FDR_q_value = p.adjust(P_value, "BH"), Bonferroni_significant = P_value < .05 / n(), FDR_significant = FDR_q_value < .05, Direction = if_else(Beta > 0, "Higher in diabetes", "Lower in diabetes")) |> arrange(P_value)
  write_csv(result, 12, "metabolite_diabetes_associations.csv")
  write.table(result$Metabolite[result$Bonferroni_significant], file.path(out_dir(12), "bonferroni_significant_metabolites.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write_csv(data.frame(Covariate = model_covariates, Included = TRUE), 12, "model_specification.csv")

  threshold <- -log10(.05 / nrow(result))
  labels <- head(result, 12)
  volcano <- ggplot(result, aes(Beta, -log10(P_value), colour = Bonferroni_significant)) + geom_point(alpha = .8, size = 2) + geom_hline(yintercept = threshold, linetype = 2, colour = "#B22222") + geom_vline(xintercept = 0, colour = "grey45") + geom_text(data = labels, aes(label = Metabolite), check_overlap = TRUE, hjust = 0, nudge_x = .02, size = 3, show.legend = FALSE) + scale_colour_manual(values = c(`FALSE` = "#78909C", `TRUE` = "#006D77")) + labs(title = "Metabolite associations with diabetes", subtitle = "Linear models adjusted for sex and PCs 1--5; dashed line is Bonferroni α = 0.05/136", x = "Adjusted mean difference (diabetes vs non-diabetes; SD units)", y = expression(-log[10](p))) + analysis_theme()
  save_plot(volcano, file.path(out_dir(12), "volcano_plot.png"), 9, 6)

  forest_data <- head(result, 15) |> mutate(Metabolite = factor(Metabolite, levels = rev(Metabolite)))
  forest <- ggplot(forest_data, aes(Beta, Metabolite, colour = Direction)) + geom_vline(xintercept = 0, colour = "grey50") + geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", height = .16) + geom_point(size = 2.5) + scale_colour_manual(values = c("Higher in diabetes" = "#B22222", "Lower in diabetes" = "#0072B2")) + labs(title = "Top adjusted metabolite associations", subtitle = "Point: β; line: 95% confidence interval", x = "Mean difference (SD units)", y = NULL) + analysis_theme()
  save_plot(forest, file.path(out_dir(12), "top_associations_forest.png"), 9, 6)

  top_features <- head(result$Metabolite, 12)
  long <- dat |> select(Diabetes, all_of(top_features)) |> pivot_longer(-Diabetes, names_to = "Metabolite", values_to = "Abundance") |> mutate(Group = factor(Diabetes, levels = c(0, 1), labels = c("Non-diabetic", "Diabetic")))
  comparison <- ggplot(long, aes(Group, Abundance, colour = Group)) + geom_boxplot(outlier.shape = NA, alpha = .12) + geom_jitter(width = .13, alpha = .55, size = 1.2) + facet_wrap(~Metabolite, scales = "free_y", ncol = 4) + labs(title = "Individual metabolite values by diabetes status", subtitle = "Top 12 adjusted associations; points are individual participants", x = NULL, y = "Z-scored abundance") + scale_colour_manual(values = c("Non-diabetic" = "#6C757D", "Diabetic" = "#D55E00")) + analysis_theme() + theme(legend.position = "none", strip.text = element_text(size = 8))
  save_plot(comparison, file.path(out_dir(12), "top_metabolites_by_diabetes.png"), 10, 8)
  write_session_info(file.path(out_dir(12), "sessionInfo.txt"))
}

classification_metrics <- function(truth, probability, threshold = .5) {
  truth <- factor(truth, levels = c(0, 1))
  predicted <- factor(ifelse(probability >= threshold, 1, 0), levels = c(0, 1))
  cm <- caret::confusionMatrix(predicted, truth, positive = "1")
  data.frame(Accuracy = unname(cm$overall[["Accuracy"]]), Sensitivity = unname(cm$byClass[["Sensitivity"]]), Specificity = unname(cm$byClass[["Specificity"]]), Balanced_accuracy = unname(cm$byClass[["Balanced Accuracy"]]), AUC = as.numeric(pROC::auc(pROC::roc(truth, probability, levels = c("0", "1"), direction = "<", quiet = TRUE))))
}

phase13 <- function() {
  dat <- read_qc_data(scaled = FALSE)
  features <- names(dat)[-(1:2)]
  y <- factor(dat$Diabetes, levels = c(0, 1))
  set.seed(20260722)
  train_index <- caret::createDataPartition(y, p = .80, list = FALSE)[, 1]
  train <- dat[train_index, ]; test <- dat[-train_index, ]
  pre <- caret::preProcess(train[, features], method = c("center", "scale"))
  train_x <- as.matrix(predict(pre, train[, features])); test_x <- as.matrix(predict(pre, test[, features]))
  train_y <- factor(train$Diabetes, levels = c(0, 1))
  test_y <- factor(test$Diabetes, levels = c(0, 1))

  set.seed(20260722)
  elastic <- glmnet::cv.glmnet(train_x, as.numeric(as.character(train_y)), family = "binomial", alpha = 1, type.measure = "auc", nfolds = 5)
  elastic_prob <- as.numeric(predict(elastic, test_x, s = "lambda.1se", type = "response"))
  set.seed(20260722)
  rf <- randomForest::randomForest(x = train_x, y = train_y, ntree = 1500, mtry = max(1, floor(sqrt(ncol(train_x)))), importance = TRUE)
  rf_prob <- predict(rf, test_x, type = "prob")[, "1"]
  metrics <- bind_rows(cbind(Model = "Elastic-net logistic regression", classification_metrics(test_y, elastic_prob)), cbind(Model = "Random forest", classification_metrics(test_y, rf_prob)))
  write_csv(metrics, 13, "classification_metrics.csv")
  prediction <- data.frame(main_id = test$main_id, Actual = test$Diabetes, Elastic_net_probability = elastic_prob, Elastic_net_prediction = as.integer(elastic_prob >= .5), Random_forest_probability = rf_prob, Random_forest_prediction = as.integer(rf_prob >= .5))
  write_csv(prediction, 13, "test_set_predictions.csv")
  write.table(features, file.path(out_dir(13), "classification_predictors.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  coef <- as.matrix(coef(elastic, s = "lambda.1se"))[-1, , drop = FALSE]
  importance <- full_join(data.frame(Metabolite = rownames(coef), Elastic_net_coefficient = coef[, 1]), data.frame(Metabolite = rownames(randomForest::importance(rf)), Random_forest_MeanDecreaseGini = randomForest::importance(rf)[, "MeanDecreaseGini"]), by = "Metabolite") |> mutate(Abs_elastic_net_coefficient = abs(Elastic_net_coefficient)) |> arrange(desc(Random_forest_MeanDecreaseGini))
  write_csv(importance, 13, "metabolite_importance.csv")
  importance_plot <- head(importance, 20) |> mutate(Metabolite = factor(Metabolite, levels = rev(Metabolite))) |> ggplot(aes(Random_forest_MeanDecreaseGini, Metabolite)) + geom_col(fill = "#006D77") + labs(title = "Random-forest metabolite importance", subtitle = "Importance estimated on the training set; validate on future external data", x = "Mean decrease in Gini", y = NULL) + analysis_theme()
  save_plot(importance_plot, file.path(out_dir(13), "metabolite_importance.png"), 8, 7)
  roc_data <- bind_rows(data.frame(FPR = 1 - pROC::roc(test_y, elastic_prob, levels = c("0", "1"), direction = "<", quiet = TRUE)$specificities, TPR = pROC::roc(test_y, elastic_prob, levels = c("0", "1"), direction = "<", quiet = TRUE)$sensitivities, Model = "Elastic-net"), data.frame(FPR = 1 - pROC::roc(test_y, rf_prob, levels = c("0", "1"), direction = "<", quiet = TRUE)$specificities, TPR = pROC::roc(test_y, rf_prob, levels = c("0", "1"), direction = "<", quiet = TRUE)$sensitivities, Model = "Random forest"))
  roc_plot <- ggplot(roc_data, aes(FPR, TPR, colour = Model)) + geom_line(linewidth = 1) + geom_abline(linetype = 2, colour = "grey50") + coord_equal() + labs(title = "Held-out test-set ROC curves", subtitle = "Small test set: estimates are imprecise", x = "False-positive rate", y = "True-positive rate") + scale_colour_manual(values = c("Elastic-net" = "#0072B2", "Random forest" = "#D55E00")) + analysis_theme()
  save_plot(roc_plot, file.path(out_dir(13), "test_set_roc.png"), 6, 6)
  write_session_info(file.path(out_dir(13), "sessionInfo.txt"))
}

phase14 <- function() {
  dat <- read_qc_data(scaled = FALSE)
  x <- scale(dat[, -(1:2), drop = FALSE])
  n <- nrow(x); p <- ncol(x)
  edge_threshold <- .20

  # Ridge regularisation limits the covariance condition number to approximately 100.
  # The resulting graph is exploratory; the edge threshold is not a formal significance test.
  s <- cov(x)
  eig <- eigen(s, symmetric = TRUE, only.values = TRUE)$values
  ridge <- max(0, (max(eig) - 100 * min(eig)) / 99) + 1e-6
  precision <- solve(s + diag(ridge, p))
  partial <- -cov2cor(precision)
  diag(partial) <- 1

  edge_index <- which(upper.tri(partial) & abs(partial) >= edge_threshold, arr.ind = TRUE)
  edges <- data.frame(
    from = colnames(x)[edge_index[, 1]],
    to = colnames(x)[edge_index[, 2]],
    partial_correlation = partial[edge_index],
    abs_partial_correlation = abs(partial[edge_index]),
    sign = ifelse(partial[edge_index] > 0, "Positive", "Negative")
  ) |>
    arrange(desc(abs_partial_correlation))
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = colnames(x)))
  graph_degree <- igraph::degree(graph)
  graph_strength <- igraph::strength(graph, weights = igraph::E(graph)$abs_partial_correlation)
  graph_betweenness <- igraph::betweenness(
    graph, directed = FALSE,
    weights = 1 / pmax(igraph::E(graph)$abs_partial_correlation, 1e-6),
    normalized = TRUE
  )
  graph_components <- igraph::components(graph)

  active_graph <- igraph::induced_subgraph(graph, vids = igraph::V(graph)[graph_degree > 0])
  communities <- igraph::cluster_louvain(active_graph, weights = igraph::E(active_graph)$abs_partial_correlation)
  raw_membership <- igraph::membership(communities)
  community_sizes <- sort(table(raw_membership), decreasing = TRUE)
  community_map <- setNames(seq_along(community_sizes), names(community_sizes))
  ordered_membership <- as.integer(community_map[as.character(raw_membership)])
  names(ordered_membership) <- names(raw_membership)

  nodes <- data.frame(
    Metabolite = igraph::V(graph)$name,
    Display_name = gsub("_", " ", igraph::V(graph)$name),
    Degree = as.integer(graph_degree),
    Strength = as.numeric(graph_strength),
    Betweenness = as.numeric(graph_betweenness),
    Component = as.integer(graph_components$membership),
    Community = 0L,
    Isolate = graph_degree == 0,
    stringsAsFactors = FALSE
  )
  nodes$Community[match(names(ordered_membership), nodes$Metabolite)] <- ordered_membership
  hub_order <- order(-nodes$Degree, -nodes$Strength, -nodes$Betweenness)
  hub_names_for_labels <- nodes$Metabolite[head(hub_order, 10)]
  nodes$Hub_label <- ifelse(nodes$Metabolite %in% hub_names_for_labels, nodes$Display_name, "")

  set.seed(20260722)
  layout <- igraph::layout_with_fr(
    active_graph,
    weights = igraph::E(active_graph)$abs_partial_correlation,
    niter = 4000
  )
  layout <- igraph::norm_coords(layout, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  node_plot <- data.frame(Metabolite = igraph::V(active_graph)$name, x = layout[, 1], y = layout[, 2]) |>
    left_join(nodes, by = "Metabolite") |>
    mutate(Community = factor(Community))
  active_edges <- igraph::as_data_frame(active_graph, what = "edges")
  edge_plot <- active_edges |>
    left_join(transmute(node_plot, from = Metabolite, x_from = x, y_from = y), by = "from") |>
    left_join(transmute(node_plot, to = Metabolite, x_to = x, y_to = y), by = "to")
  label_nodes <- node_plot |>
    arrange(desc(Degree), desc(Strength)) |>
    slice_head(n = 10)

  network_plot <- ggplot() +
    geom_segment(
      data = edge_plot,
      aes(x_from, y_from, xend = x_to, yend = y_to, colour = sign, linewidth = abs_partial_correlation, alpha = abs_partial_correlation),
      lineend = "round"
    ) +
    geom_point(
      data = node_plot,
      aes(x, y, size = Degree, fill = Community),
      shape = 21, colour = "white", stroke = .6
    ) +
    ggrepel::geom_text_repel(
      data = label_nodes,
      aes(x, y, label = Display_name),
      seed = 42, size = 3.25, fontface = "bold", colour = "#17324D",
      box.padding = .55, point.padding = .35, min.segment.length = 0,
      max.overlaps = Inf, max.time = 10, max.iter = 20000, force = 5, force_pull = .1, segment.colour = "#8A98A8", segment.alpha = .7
    ) +
    scale_colour_manual(values = c(Positive = "#087E8B", Negative = "#D1495B"), name = "Edge sign") +
    scale_fill_viridis_d(option = "C", begin = .08, end = .9, guide = "none") +
    scale_size_continuous(range = c(2.5, 9), breaks = c(2, 5, 8, 11), name = "Node degree") +
    scale_linewidth_continuous(range = c(.25, 1.7), breaks = c(.20, .25, .30, .35), name = "|Partial r|") +
    scale_alpha_continuous(range = c(.18, .72), guide = "none") +
    coord_equal(clip = "off") +
    labs(
      title = "Metabolite partial-correlation network",
      subtitle = sprintf("%d connected metabolites · %d edges at |partial r| ≥ %.2f · labels show 10 topology hubs", igraph::vcount(active_graph), igraph::ecount(active_graph), edge_threshold),
      caption = "Ridge-regularised estimate. Communities reflect graph topology, not curated pathways. Isolates remain in the exported node table."
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 20, colour = "#17324D"),
      plot.subtitle = element_text(size = 11.5, colour = "#526777", margin = margin(b = 12)),
      plot.caption = element_text(size = 9, colour = "#66788A", hjust = 0),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = margin(18, 28, 14, 28)
    ) +
    guides(
      size = guide_legend(order = 1),
      colour = guide_legend(order = 2, override.aes = list(linewidth = 1.4)),
      linewidth = guide_legend(order = 3)
    )
  save_plot(network_plot, file.path(out_dir(14), "partial_correlation_network.png"), 13, 10)

  strongest_edges <- slice_head(edges, n = min(30, nrow(edges)))
  strong_graph <- igraph::graph_from_data_frame(strongest_edges, directed = FALSE)
  set.seed(20260722)
  strong_layout <- igraph::layout_with_kk(strong_graph, weights = 1 / igraph::E(strong_graph)$abs_partial_correlation)
  strong_layout <- igraph::norm_coords(strong_layout, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  strong_nodes <- data.frame(Metabolite = igraph::V(strong_graph)$name, x = strong_layout[, 1], y = strong_layout[, 2]) |>
    left_join(nodes, by = "Metabolite") |>
    mutate(Community = factor(Community))
  strong_edge_plot <- igraph::as_data_frame(strong_graph, what = "edges") |>
    left_join(transmute(strong_nodes, from = Metabolite, x_from = x, y_from = y), by = "from") |>
    left_join(transmute(strong_nodes, to = Metabolite, x_to = x, y_to = y), by = "to")
  strong_plot <- ggplot() +
    geom_segment(
      data = strong_edge_plot,
      aes(x_from, y_from, xend = x_to, yend = y_to, colour = sign, linewidth = abs_partial_correlation),
      alpha = .65, lineend = "round"
    ) +
    geom_point(data = strong_nodes, aes(x, y, fill = Community), shape = 21, size = 5.2, colour = "white", stroke = .7) +
    ggrepel::geom_text_repel(
      data = strong_nodes, aes(x, y, label = Display_name), seed = 42,
      size = 3.15, colour = "#17324D", box.padding = .5, point.padding = .3,
      min.segment.length = 0, max.overlaps = Inf, max.time = 5, force = 1.5, segment.colour = "#9AA7B5"
    ) +
    scale_colour_manual(values = c(Positive = "#087E8B", Negative = "#D1495B"), name = "Edge sign") +
    scale_fill_viridis_d(option = "C", begin = .08, end = .9, guide = "none") +
    scale_linewidth_continuous(range = c(.7, 2.3), breaks = c(.28, .30, .32, .34, .36), name = "|Partial r|") +
    coord_equal(clip = "off") +
    labs(
      title = "Thirty strongest conditional metabolite relationships",
      subtitle = "Focused view for readable metabolite labels; edge colour shows sign and width shows absolute partial correlation",
      caption = "This focused view is a visual subset of the complete |partial r| ≥ 0.20 network."
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 20, colour = "#17324D"),
      plot.subtitle = element_text(size = 11.5, colour = "#526777", margin = margin(b = 12)),
      plot.caption = element_text(size = 9, colour = "#66788A", hjust = 0),
      legend.position = "bottom",
      plot.margin = margin(20, 35, 16, 35)
    )
  save_plot(strong_plot, file.path(out_dir(14), "partial_correlation_network_strongest_edges.png"), 13, 10)

  sensitivity <- bind_rows(lapply(c(.15, .20, .25, .30), function(threshold) {
    idx <- which(upper.tri(partial) & abs(partial) >= threshold, arr.ind = TRUE)
    threshold_edges <- data.frame(from = colnames(x)[idx[, 1]], to = colnames(x)[idx[, 2]])
    threshold_graph <- igraph::graph_from_data_frame(threshold_edges, directed = FALSE, vertices = data.frame(name = colnames(x)))
    threshold_degree <- igraph::degree(threshold_graph)
    threshold_components <- igraph::components(threshold_graph)
    data.frame(
      Threshold = threshold,
      Edges = igraph::ecount(threshold_graph),
      Connected_metabolites = sum(threshold_degree > 0),
      Isolates = sum(threshold_degree == 0),
      Components = threshold_components$no,
      Largest_component = max(threshold_components$csize)
    )
  }))
  hub_table <- nodes |>
    filter(!Isolate) |>
    arrange(desc(Degree), desc(Strength), desc(Betweenness)) |>
    slice_head(n = 20)
  hub_names <- paste(head(hub_table$Display_name, 5), collapse = ", ")
  interpretation <- paste0(
    "At |partial r| ≥ 0.20, the ridge-regularised network contains ", nrow(edges),
    " edges among ", igraph::vcount(active_graph), " connected metabolites, with ", sum(nodes$Isolate),
    " isolate. The highest-degree hubs are ", hub_names,
    ". Several prominent hubs are lipid-related metabolites, while amino-acid and nucleotide intermediates also occur among central nodes, a qualitatively plausible pattern for coordinated metabolism. However, communities are topology-derived rather than curated pathway assignments, the 0.20 threshold is heuristic, and no bootstrap stability or external replication was performed; hub and pathway interpretations therefore remain exploratory."
  )

  write.csv(partial, file.path(out_dir(14), "ridge_partial_correlation_matrix.csv"))
  igraph::write_graph(graph, file.path(out_dir(14), "metabolite_network.graphml"), format = "graphml")
  write_csv(edges, 14, "cytoscape_edge_list.csv")
  write_csv(nodes, 14, "cytoscape_node_table.csv")
  write_csv(hub_table, 14, "network_hubs.csv")
  write_csv(sensitivity, 14, "network_threshold_sensitivity.csv")
  writeLines(interpretation, file.path(out_dir(14), "network_interpretation.txt"))
  write_csv(data.frame(
    Samples = n,
    Metabolites = p,
    Ridge_penalty = ridge,
    Edge_threshold = edge_threshold,
    Edges = nrow(edges),
    Connected_metabolites = igraph::vcount(active_graph),
    Isolates = sum(nodes$Isolate),
    Positive_edges = sum(edges$sign == "Positive"),
    Negative_edges = sum(edges$sign == "Negative"),
    Topology_communities = length(unique(ordered_membership)),
    Method = "ridge-regularised partial correlation"
  ), 14, "network_parameters.csv")
  write_session_info(file.path(out_dir(14), "sessionInfo.txt"))
}

phase <- Sys.getenv("BIOINF_PHASE", unset = if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else "")
switch(phase, `11` = phase11(), `12` = phase12(), `13` = phase13(), `14` = phase14(), stop("Usage: Rscript professional_pipeline.R {11|12|13|14}"))

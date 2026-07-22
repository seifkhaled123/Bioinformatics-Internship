source("professional_helpers.R")

source("functions.R")
require_packages(c("readxl", "dplyr", "tidyr", "ggplot2", "glmnet", "randomForest", "pROC", "caret", "igraph"))
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
  # Ridge-regularised precision matrix: stable for p close to n. This is an exploratory network.
  s <- cov(x)
  eig <- eigen(s, symmetric = TRUE, only.values = TRUE)$values
  ridge <- max(0, (max(eig) - 100 * min(eig)) / 99) + 1e-6
  precision <- solve(s + diag(ridge, p))
  partial <- -cov2cor(precision); diag(partial) <- 1
  edge_index <- which(upper.tri(partial) & abs(partial) >= .20, arr.ind = TRUE)
  edges <- data.frame(from = colnames(x)[edge_index[, 1]], to = colnames(x)[edge_index[, 2]], partial_correlation = partial[edge_index], sign = ifelse(partial[edge_index] > 0, "positive", "negative")) |> arrange(desc(abs(partial_correlation)))
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = colnames(x)))
  nodes <- data.frame(Metabolite = igraph::V(graph)$name, Degree = igraph::degree(graph), Betweenness = igraph::betweenness(graph, normalized = TRUE), Component = igraph::components(graph)$membership)
  write.csv(partial, file.path(out_dir(14), "ridge_partial_correlation_matrix.csv"))
  write_csv(edges, 14, "cytoscape_edge_list.csv")
  write_csv(nodes, 14, "cytoscape_node_table.csv")
  write_csv(data.frame(Samples = n, Metabolites = p, Ridge_penalty = ridge, Edge_threshold = .20, Edges = nrow(edges), Method = "ridge-regularised partial correlation"), 14, "network_parameters.csv")
  set.seed(20260722)
  png(file.path(out_dir(14), "partial_correlation_network.png"), width = 2400, height = 2100, res = 260)
  plot(graph, layout = igraph::layout_with_fr(graph), vertex.size = 4 + sqrt(igraph::degree(graph)) * 2, vertex.label.cex = .55, vertex.color = "#78B7C5", edge.width = 1 + 4 * abs(igraph::E(graph)$partial_correlation), edge.color = ifelse(igraph::E(graph)$sign == "positive", "#0072B2", "#B22222"), main = "Ridge partial-correlation network (|r| ≥ 0.20)")
  dev.off()
  write_session_info(file.path(out_dir(14), "sessionInfo.txt"))
}

phase <- Sys.getenv("BIOINF_PHASE", unset = if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else "")
switch(phase, `11` = phase11(), `12` = phase12(), `13` = phase13(), `14` = phase14(), stop("Usage: Rscript professional_pipeline.R {11|12|13|14}"))

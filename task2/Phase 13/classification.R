library(caret)
library(randomForest)
library(pROC)

# 1. Load your Data
data <- read.csv("../Phase 11/outputs/Qatari_metabolomics_Cleaned_Scaled.csv") 

# 2. Format the target variable using the correct column name
data$Diabetes <- as.factor(data$Diabetes)

# 3. CRITICAL: Remove the ID column so the model only trains on metabolites
data <- data[, -which(names(data) == "main_id")]

# 4. Data Splitting (80% Train, 20% Test)
set.seed(42) 
trainIndex <- createDataPartition(data$Diabetes, p = 0.8, list = FALSE)
train_data <- data[trainIndex, ]
test_data  <- data[-trainIndex, ]

# 5. Train Models
# Model A: Logistic Regression
log_model <- glm(Diabetes ~ ., data = train_data, family = "binomial")

# Model B: Random Forest
rf_model <- randomForest(Diabetes ~ ., data = train_data, importance = TRUE, ntree = 500)

# 6. Evaluate Performance on the TEST SET
# Get predictions and probabilities
rf_preds <- predict(rf_model, test_data)
rf_probs <- predict(rf_model, test_data, type = "prob")[,2] 

# Confusion Matrix (Accuracy, Sensitivity, Specificity)
conf_matrix <- confusionMatrix(rf_preds, test_data$Diabetes)
print(conf_matrix)

# Calculate AUC
roc_curve <- roc(test_data$Diabetes, rf_probs, quiet = TRUE)
auc_value <- auc(roc_curve)
cat("\n==============================\n")
cat("Random Forest AUC:", auc_value, "\n")
cat("==============================\n")

# 7. Generate Deliverables
# Variable Importance Plot
varImpPlot(rf_model, main="Random Forest: Top Predictive Metabolites", n.var = 15)

# Save test predictions table
test_results <- data.frame(
  Actual = test_data$Diabetes,
  Predicted = rf_preds,
  Probability = rf_probs
)
write.csv(test_results, "outputs/Phase13_Test_Set_Predictions.csv", row.names = FALSE)
source('../functions.R')
library(ggrepel)


run_plink("--bfile ../data/dataQC --pca 10 --out outputs/PCA", step_name = "PCA")

pcs <- read_plink_out("outputs/PCA.eigenvec", header = FALSE)
eig <- read_plink_out("outputs/PCA.eigenval", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

var_exp <- eig$V1 / sum(eig$V1) * 100
jpeg("outputs/ScreePlot.jpg", 900, 650)
plot(var_exp, type = "b", pch = 19,
     xlab = "Principal Component", ylab = "Variance Explained (%)",
     main = "Scree Plot")
dev.off()

p1 <- ggplot(pcs, aes(PC1, PC2)) +
  geom_point(size = 2, color = "steelblue") +
  theme_minimal() +
  labs(title = "Population Structure")
ggsave("outputs/PCA_PC1_PC2.jpg", p1, width = 7, height = 5)
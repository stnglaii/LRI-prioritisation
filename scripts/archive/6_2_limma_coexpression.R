library(limma)

# Assuming your data frame 'data' has columns: 
# interaction, sampleID, condition (Cancer/Healthy), mean

# Pivot data to wide format
library(tidyr)
expression_data <- data %>%
  select(interaction, sampleID, mean) %>%
  spread(key = sampleID, value = mean)

# Set row names and remove the interaction column
rownames(expression_data) <- expression_data$interaction
expression_data <- expression_data[,-1]

# Check data
head(expression_data)

# Log-transform if necessary
# Adding a small constant to avoid log(0) if needed
expression_data_log <- log2(expression_data + 0.1)

# Create a design matrix
sample_info <- unique(data[, c("sampleID", "condition")])
condition <- factor(sample_info$condition)
design <- model.matrix(~0 + condition)  # Use model without intercept
colnames(design) <- levels(condition)

# Check design matrix
design

# Fit the linear model
fit <- lmFit(expression_data_log, design)

# Define contrasts (Cancer vs Healthy)
contrast_matrix <- makeContrasts(Cancer_vs_Healthy = Cancer - Healthy, levels=design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Examine top differentially expressed interactions
results <- topTable(fit2, coef="Cancer_vs_Healthy", number=Inf, adjust.method="BH")

# View results
head(results)
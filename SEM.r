# Load packages
library(lavaan)
library(semPlot)
library(dplyr)
# Data loading
data <- read.csv("data_file.csv", header = TRUE, row.names = 1)

all.scaled <- data.frame(scale(data, center = FALSE)) # Standardization

# 1. Model specification
# AvaiK~ ph + TOC + TN
# AvaiZn~ ph + TOC + TN
# TOC~ WaterContent + RealTemp + Latitude + Longitude + ph
# TN~ WaterContent + RealTemp + Latitude + Longitude + ph
# top3_bac_sum~ ph + C.N + AvaiP + AvaiCu
model <- '
ph~ MAT + MAP + Latitude
AvaiP~ Latitude + ph 
Ca~ ph + AvaiP + Latitude
total_ARG~ ph + AvaiP + Ca + Thermoproteota + integrase + Zinc.Cd
'

model1 <- '
ph~ WaterContent + RealTemp + Latitude
AvaiP~ WaterContent + Latitude + ph + RealTemp
Ca~ ph + Latitude + WaterContent + AvaiP + RealTemp
Thermoproteota~ ph + WaterContent + RealTemp + AvaiP + Ca + Latitude
Zinc.Cd~ Ca + ph + Thermoproteota + WaterContent + RealTemp + Latitude
total_ARG~ ph + Latitude + Ca + Zinc.Cd + AvaiP + WaterContent + RealTemp + Thermoproteota
'

# [Other model definitions (model2 to model10) remain unchanged with same structure...]

# 2. Fit the model (SEM)
fit_model <- sem(model1, data = all.scaled, se = 'bootstrap')

# 3. Display summary output
summary(fit_model, standardize = TRUE, rsq = TRUE, modindices = TRUE)
regressions <- parameterEstimates(fit_model, standardized = TRUE)
total_ARG_regressions <- subset(regressions, lhs == "total_ARG" & op == "~")
total_ARG_regressions <- total_ARG_regressions[, c("lhs", "op", "rhs", "pvalue", "std.all")]
total_ARG_regressions

# Iterative selection and comparison of metrics
fitmeasures(fit_model, c("chisq", "df", "pvalue", "gfi", "cfi", "rmsea", "rmr", "srmr", "AIC"))

# Visualization
semPaths(fit_model, "std", edge.label.cex = 1.5, sizeMan = 10, nCharNodes = 5,
         fade = FALSE,
         layout = "spring",
         optimizeLatRes = FALSE, residuals = FALSE)

# Check modification indices
mf <- modificationIndices(fit_model)
mf <- mf[order(mf$mi, decreasing = TRUE), ]
head(mf)

# Key interpretation guidelines:
# - Results require p-value > 0.05 for good fit
# - Chi-square/df < 2 indicates acceptable model-data match
# - GFI > 0.90, CFI > 0.90, RMSEA < 0.06 suggest good fit
# - MI values help identify potential model improvements

# Created by: Vladislav Vlasov
# Created on: 05/05/2021

load("./cortex.rdata")

# Study and describe the data. Do you see indications of potential issues when sta- tistically modeling the data?

# Train and compare Ridge and LASSO models to separate the C/S from S/C (the Behavior variable) samples based on protein
# expression. Interpret the results of the optimizations:
# 1 — Do correlations between variables influence the results? How?
# 2 — Explain the shape of the trend in the cross-validation plot for selecting the op- timal lambda value.
# 3 —Can a reduced set of variables predict the Behavior variable?

# Build a boosting model and compare it to the ridge and lasso models. Are the same variables important for the
# predictions? Do you see evidence for non-linear effects or interactions between the most important predictor
# variables?

# Does Memantine injection influence protein values when controlling for genotype and treatment?

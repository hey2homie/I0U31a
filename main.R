# Created by: Vladislav Vlasov
# Created on: 05/05/2021

library(dplyr)
library(tidyr)
library(expss)		# Better table
library(glmnet)		# Ridge and Lasso regressions
library(caret)		# Dataset partition
library(corrplot)	# Correlation plot
library(tree)		# Decision Trees
library(gbm)		# Gradient Boosting
library(dismo)		# Cross-Validation of number of trees

options(warn = 0)

# Loading provided dataset
load("./cortex.rdata")
cortex$Behavior <- ifelse(cortex$Behavior == "C/S", 1, 0)

# Setting seed for reproduction
set.seed(10)

# ----------------------------------------------------------------------------------------------------------------------
# Study and describe the data. Do you see indications of potential issues when statistically modeling the data?
# ----------------------------------------------------------------------------------------------------------------------

# Dimensions
dim(cortex)

# ----------------------------------------------------------------------------------------------------------------------
# We can see that it's high dimensional data — 70 observations and 73 parameters. This will result in difficulties with
# the fitting, hence, will limit prediction power of the linear regression models. Issue like this can be fixed by
# variable selection methods.
# ----------------------------------------------------------------------------------------------------------------------


# Classes distributions
cortex %>%
	tab_cells(Genotype) %>%
	tab_cols(Treatment, Behavior) %>%
	tab_stat_cases() %>%
	tab_pivot()

# Histogram
cortex[1:70] %>%
	as.matrix() %>%
	hist(main = "Histrogram of gene expression across the dataset", xlab = "Value")

# ----------------------------------------------------------------------------------------------------------------------
# We can see that the distrubition of the data is not near close to the normal distribution. We need to transform it
# with, for example, logarithmic transformation.
# ----------------------------------------------------------------------------------------------------------------------

# Histrogram with logarithmic transformation
cortex[1:70] %>%
	as.matrix() %>%
	log() %>%
	hist(main = "Histrogram of gene expression across the dataset after logarithmic transformation", xlab = "Value")

# ----------------------------------------------------------------------------------------------------------------------
# Now the distribution is more close to the normal distribution, hence, we can apply it to the data.
# ----------------------------------------------------------------------------------------------------------------------

cortex[1:70] <- log(cortex[1:70])

# ----------------------------------------------------------------------------------------------------------------------

# Correlation plot
cor.plot <- cortex[1:70] %>%
	cor() %>%
	corrplot()

# ----------------------------------------------------------------------------------------------------------------------
# Train and compare Ridge and LASSO models to separate the C/S from S/C (the Behavior variable) samples based on protein
# expression. Interpret the results of the optimizations:
# 1 — Do correlations between variables influence the results? How?
# 2 — Explain the shape of the trend in the cross-validation plot for selecting the optimal lambda value.
# 3 — Can a reduced set of variables predict the Behavior variable?
# ----------------------------------------------------------------------------------------------------------------------

# Splitting initial dataset into training and test data
training.samples <- cortex$Behavior %>%
  createDataPartition(p = 0.7, list = FALSE)
train.data  <- cortex %>%
	{ subset(., select = -c(71, 72)) } %>%
	{.[training.samples, ]}
test.data <- cortex %>%
	{ subset(., select = -c(71, 72)) } %>%
	{.[-training.samples, ]}

# Making matrix of predictors and converting factor variable into numerical one
predictors.train <- model.matrix(Behavior ~., train.data)
responce.train <- train.data$Behavior

# Choosing appropriate value of lambda by 10 fold cross-validation and type measure deviance, which properly suites
# logistric regressions
calculate.lambda <- function(alpha, predictors = predictors.train, responce = responce.train) {
	return(cv.glmnet(predictors, responce,
					 alpha = alpha, family = "binomial", type.measure = "deviance"))
}

# Making either Ridge or Lasso logistic regression model for classification of Behavior. Predictors are scaled by
# default in glmnet function.
make.model <- function(alpha, lambda.min, predictors = predictors.train, responce = responce.train) {
	return(glmnet(predictors, responce,
				  alpha = alpha, family = "binomial", lambda = lambda.min))
}

# Testing model accuracy with type "response" since it's logistic regression
test.model <- function(model, test = test.data) {
	return(
		model %>%
		{ predict(., newx = model.matrix(Behavior ~., test), type = "response") }  %>%
		{ ifelse(. > 0.5, 1, 0) } %>%
		{ mean(. == test$Behavior) }
	)
}

# ----------------------------------------------------------------------------------------------------------------------
ridge.lamba <- calculate.lambda(0)
lasso.lamba <- calculate.lambda(1)

ridge.model <- make.model(0, ridge.lamba$lambda.min)
lasso.model <- make.model(1, lasso.lamba$lambda.min)

ridge.accuracy <- test.model(ridge.model)
lasso.accuracy <- test.model(lasso.model)

# Printing used lambda value for each model
cat("Lambda used in Ridge regression: ", ridge.lamba$lambda.min, "\n",
    "Lambda used in Lasso regression: ", lasso.lamba$lambda.min, "\n",
	sep = "")

# Printing amount of predictors used by both models
cat("Predictors used in Ridge regression: ", length(summary(coef(ridge.lamba, ridge.lamba$lambda.min))$x), "\n",
    "Predictors used in Lasso regression: ", length(summary(coef(lasso.lamba, lasso.lamba$lambda.min))$x), "\n",
	sep = "")

# Printing proteins that are used in making prediction in Lasso regression
lasso.usedProteins <- colnames(cortex)[summary(coef(lasso.lamba, lasso.lamba$lambda.min))$i] %>%
	{ .[seq_along(.) - 1] } %>%
	print()

findCorrelation(cor.plot, cutoff = 0.8) %>%
	sort() %>%
	{ which(. == summary(coef(lasso.lamba, lasso.lamba$lambda.min))$i) } %>%
	{ colnames(cortex)[.] }

findCorrelation(cor.plot, cutoff = 0.8) %>%
	{ print(length(.)) }

# ----------------------------------------------------------------------------------------------------------------------
# 1 — Since in Lasso regression we are also performing variable selection in addition to shrinkige, and after running
# the model we can see that it has only 13 variables in it of which one is highly correlated (correlation > 0.8). The
# whole dataset include 27 highly correlated variables which are also included in Ridge regression. This can result in
# multicollinearity which, in turn, affects model coefficients, hence, prediction power.
# ----------------------------------------------------------------------------------------------------------------------

par(mfrow = c(1, 2))
plot(ridge.lamba)
plot(ridge.lamba$glmnet.fit, "lambda", label = FALSE)

par(mfrow = c(1, 2))
plot(lasso.lamba)
plot(lasso.lamba$glmnet.fit, "lambda", label = FALSE)

# ----------------------------------------------------------------------------------------------------------------------
# 2 — On the 1st plot (Ridgre regression) we can see how the deviance increases with increase in lambda. The more
# lambda grows, more penalised model becames. We are choosing the lowest point to determine proper value of lambda.
# On ther 2nd plot (Lasso regression) we can see how value of lambda changes with the decrease in predictors. Deviance
# grows in the directon of increased amount of predictors. As before, we are choosing lowest value of lambda.
# ----------------------------------------------------------------------------------------------------------------------

# Printing models accuracy
cat("Ridge regression accuracy: ", ridge.accuracy, "\n",
    "Lasso regression accuracy: ", lasso.accuracy, "\n",
	"Difference in accuracy: ", Mod(ridge.accuracy - lasso.accuracy), "\n",
	sep = "")

# ----------------------------------------------------------------------------------------------------------------------
# 3 — As we previously saw, model in which Ridge regression is used uses all available predictors (plus intercept).
# Lasso regression, on the other hand, uses only 13 of them (plus intercept). Reduced set of predictors resulted in
# improving accuracy by 0.095. Hence, we can conclude that Lasso works better. However, since lambda used in Lasso
# regression approaches zero, we can say that it acts like OLS regression and penalty doesn't effect selection and
# shrinkage.
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Build a boosting model and compare it to the ridge and lasso models. Are the same variables important for the
# predictions? Do you see evidence for non-linear effects or interactions between the most important predictor
# variables?
# ----------------------------------------------------------------------------------------------------------------------

# Building boost model using bernoulli distribution as family parameter since the output is binary. First, running
# cross-validation to determine the amount of trees that we'll use in model making
boost.bestNumber <- gbm.step(train.data, gbm.x = 1:70, gbm.y = 71, family = "bernoulli", tree.complexity = 5,
					 learning.rate = 0.005, bag.fraction = 0.5, n.folds = 10) %>%
	gbm.perf()
boost.model <- gbm(train.data$Behavior ~., train.data, distribution = "bernoulli", n.trees = boost.bestNumber,
				   shrinkage = 0.01, interaction.depth = 4, cv.folds = 10)
boost.prediction <- predict(boost.model, newdata = test.data, n.trees = boost.bestNumber,
							type = "response")
boost.accuracy <- ifelse(boost.prediction > 0.5, 1, 0) %>%
	{ mean(. == test.data$Behavior) }

cat("Accuracy of the gradient boosted decision tree model: ", boost.accuracy, "\n",
	sep = "")

# ----------------------------------------------------------------------------------------------------------------------
# We see that the prediction power of Gradient Boost Decesion Tree is the same as in the Lasso regression (1).
# ----------------------------------------------------------------------------------------------------------------------

# Vector of influential proteins used in Boost model
boost.usedProteins <- summary(boost.model)$rel.inf %>%
	{ which(. > 1) } %>%
	{ row.names(summary(boost.model))[.]}
print(length(boost.usedProteins))

# Vector of predictors that don't contribute to the model
boost.unusedProtein <- summary(boost.model)$rel.inf %>%
	{ which(. == 0) } %>%
	{ row.names(summary(boost.model))[.]}
print(length(boost.unusedProtein))

# Find common proteins used in both models
intersect(boost.usedProteins, lasso.usedProteins)

# ----------------------------------------------------------------------------------------------------------------------
# In the Boost model there are 12 higly influential parameters (> 1) and 27 that don't contribute to the model (0). Of
# 12 proteins in the Boost model only 3 (pERK_N, DYRK1A_N, nNOS_N) are used in the Lasso regression model. This can be
# an evidance towards non-linear relationships or possible interactions between predictors in the dataset.
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Does Memantine injection influence protein values when controlling for genotype and treatment?
# ----------------------------------------------------------------------------------------------------------------------
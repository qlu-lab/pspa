<h1 align="center">
<p> pspa
</h1>

This repository hosts the R package that implements the `post-prediction adaptive infernece` (`pspa`) method described in the paper: [Assumption-Lean and Data-Adaptive Post-Prediction Inference](https://arxiv.org/abs/2311.14220). 

`pspa` provides valid and powerful inference based on ML predictions for parameters defined through estimation equations.


## Installation         
```
# install.packages("devtools")
devtools::install_github("qlu-lab/pspa")
```

## Useful examples
### ML-predicted labels
Here are examples of pspa with ML-predicted labels for M-estimation tasks including: mean estimation, linear regression, logistic regression, and Poisson regrssion. The main function is `pspa_y()`, where the argument `method` indicates which task to do.


```
# Load the package
library(pspa)

# Load the simulated data
set.seed(999)
data <- sim_data_y()
X_lab = data$X_lab ## Covariates in the labeled data
X_unlab = data$X_unlab ## Covariates in the unlabeled data
Y_lab = data$Y_lab ## Observed outcome in the labeled data
Yhat_lab = data$Yhat_lab ## Predicted outcome in the labeled data
Yhat_unlab = data$Yhat_unlab ## Predicted outcome in the unlabeled data
``````

#### Mean estimation
```
# Run pspa mean estimation
fit_mean <- pspa_y(Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
                  alpha = 0.05, method = "mean")

print(fit_mean)

#   Estimate  Std.Error Lower.CI Upper.CI       P.value    Weight
# 1 1.623601 0.05514429  1.51552 1.731682 1.557956e-190 0.9226747
```

#### Linear regression
```
# Run pspa linear regression
fit_ols <- pspa_y(X_lab = X_lab, X_unlab = X_unlab,
           Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
           alpha = 0.05, method = "ols")

print(fit_ols)

#     Estimate  Std.Error Lower.CI Upper.CI       P.value    Weight
#    1.6227593 0.05385971 1.517196 1.728322 1.998620e-199 0.8741176
# X1 0.8698039 0.07193445 0.728815 1.010793  1.169598e-33 1.0000000
```

#### Logistic regression
```
# Load the simulated data
set.seed(999)
data <- sim_data_y(binary = T)
X_lab = data$X_lab
X_unlab = data$X_unlab
Y_lab = data$Y_lab
Yhat_lab = data$Yhat_lab
Yhat_unlab = data$Yhat_unlab

# Run pspa logistic regression
fit_logistic <- pspa_y(X_lab = X_lab, X_unlab = X_unlab,
                      Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
                      alpha = 0.05, method = "logistic")

print(fit_logistic)

#      Estimate  Std.Error   Lower.CI   Upper.CI      P.value    Weight
#    -0.1289001 0.08347881 -0.2925156 0.03471532 1.225626e-01 0.4290559
# X1  0.5749601 0.08653142  0.4053617 0.74455861 3.041970e-11 0.5337078
```

#### Poisson regression
```
# Load the simulated data
set.seed(999)
data <- sim_data_y()
X_lab = data$X_lab
X_unlab = data$X_unlab
Y_lab = round(data$Y_lab - min(data$Y_lab))
Yhat_lab = round(data$Yhat_lab - min(data$Yhat_lab))
Yhat_unlab = round(data$Yhat_unlab - min(Yhat_unlab))

# Run pspa Poisson regression
fit_poisson <- pspa_y(X_lab = X_lab, X_unlab = X_unlab,
                     Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
                     alpha = 0.05, method = "poisson")

print(fit_poisson)

#     Estimate  Std.Error  Lower.CI  Upper.CI      P.value    Weight
#    0.9732937 0.02261537 0.9289684 1.0176191 0.000000e+00 0.8392517
# X1 0.3188511 0.03125507 0.2575923 0.3801099 1.950752e-24 0.8303991
```
### ML-predicted covariates
Here are examples of pspa with ML-predicted covariates for M-estimation tasks including: mean estimation, linear regression, logistic regression, and Poisson regrssion. The main function is `pspa_x()`, where the argument `method` indicates which task to do.


```
# Load the package
library(pspa)

# Load the simulated data
set.seed(999)
data <- sim_data_x(r = 0.5)
X_lab = data$X_lab ## Covariates in the labeled data
Xhat_lab = data$Xhat_lab ## Predicted covariates in the labeled data
Xhat_unlab = data$Xhat_unlab ## Predicted covariates in the unlabeled data
Y_lab = data$Y_lab ## Observed outcome in the labeled data
Y_unlab = data$Y_unlab ## Observed outcome in the unlabeled data
``````

#### Linear regression
```
# Run pspa linear regression
fit_ols <- pspa_x(X_lab = X_lab, Xhat_lab = Xhat_lab,
           Xhat_unlab = Xhat_unlab, Y_lab = Y_lab, Y_unlab = Y_unlab,
           alpha = 0.05, method = "ols")

print(fit_ols)

#        Estimate   Std.Error     Lower.CI   Upper.CI   P.value    Weight
#     0.503592340 0.009579156  0.484817539 0.52236714 0.0000000 0.7553199
# X1 -0.003066491 0.031646383 -0.065092261 0.05895928 0.9228069 0.3674327
# X2  0.011997846 0.009443294 -0.006510669 0.03050636 0.2039013 0.7553199
```

#### Logistic regression
```
# Load the simulated data
set.seed(999)
data <- sim_data_x(r = 0.5, binary = T)
X_lab = data$X_lab ## Covariates in the labeled data
Xhat_lab = data$Xhat_lab ## Predicted covariates in the labeled data
Xhat_unlab = data$Xhat_unlab ## Predicted covariates in the unlabeled data
Y_lab = data$Y_lab ## Observed outcome in the labeled data
Y_unlab = data$Y_unlab ## Observed outcome in the unlabeled data

# Run pspa logistic regression
fit_logistic <-  pspa_x(X_lab = X_lab, Xhat_lab = Xhat_lab,
           Xhat_unlab = Xhat_unlab, Y_lab = Y_lab, Y_unlab = Y_unlab,
           alpha = 0.05, method = "logistic")

print(fit_logistic)


#       Estimate  Std.Error    Lower.CI   Upper.CI   P.value    Weight
#     0.01293817 0.03826099 -0.06205199 0.08792834 0.7352459 0.7562264
# X1 -0.01242760 0.12668094 -0.26071767 0.23586247 0.9218516 0.3673090
# X2  0.04379359 0.03775431 -0.03020350 0.11779067 0.2460640 0.7562264
```

## Analysis script
We provide the script for analysis in the `pspa` paper [here](https://github.com/jmiao24/POP-Inf_analysis).

## Contact 

Please submit an issue or contact Jiacheng (jiacheng.miao@wisc.edu) or Xinran (xinran.miao@wisc.edu) for questions.

## Reference
[Assumption-lean and Data-adaptive Post-Prediction Inference](https://arxiv.org/abs/2311.14220)

[Valid inference for machine learning-assisted GWAS](https://www.medrxiv.org/content/10.1101/2024.01.03.24300779v1)

## "POP" familial links
* [POP-TOOLS](https://github.com/qlu-lab/POP-TOOLS) (**PO**st-**P**rediction **TOOLS**) is a toolkit for conducting valid and powerful machine learning (ML)-assisted genetic association studies. It currently implements
  * `POP-GWAS`, where statistical and computational methods are optimized for GWAS applications.

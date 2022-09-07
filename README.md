# LRBmat

LRBmat (Binary matrix based on Logistic Regression) is proposed when there is complex interaction between covariates and heterogeneity between covariates and response variables. 

## [model_functions.R](https://github.com/tsnm1/LRBmat/blob/main/model_functions.R "model_functions.R")

This file mainly contains the subsequent required functions, including data generation function, variable selection function, heterogeneity test function, various classification function and so on.

## [model_simulation.R](https://github.com/tsnm1/LRBmat/blob/main/model_simulation.R "model_simulation.R")

This file is mainly the code of the simulation part, which mainly includes three parts: the first is the generation of data, the second is the comparison of several logistic models, and the third is the effect of binary matrix on other classification methods.

## [model_real_data_analysis.R](https://github.com/tsnm1/LRBmat/blob/main/model_real_data_analysis.R "model_real_data_analysis.R")

This part is the analysis of real world data. The main content is similar to simulation, but the difference is that the real data model includes the test of heterogeneity and the analysis of association rules.

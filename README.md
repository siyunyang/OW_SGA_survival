# OW_SGA_survival

The R code for implementing the simulation study as shown in "Propensity Score Methods for Causal Subgroup Analysis with Time-to-Event Outcomes" by Yang et al.

List of Files:

simu_main.R = R program for implementing the main simulation study.

simu_superpopulation.R = R program for calculating the true subgroup causal effect parameters used in the simulation, which takes a long time to run.

_true_val_censor_100_2021-09-13_.Rdata = R data that saves the true subgroup causal effect parameters returned from simu_superpopulation.R to evaluate bias

datagen_520.R = a helper R function to generate the simulated data

helper_overall.R, RMST_sub_overall.R = helper R functions to generates point estimates, standard errors and confidence intervals for the desired subgroup survival causal contrasts of interest. 

PSWeight.SGA = a suite of R functions to perform propensity score weighting analysis for causal subgroup analyses. For more details and illustrative examples, please refer to https://github.com/siyunyang/PSweight.SGA

# cmipb_3rd_cyz
CMI_PB 3rd Challenge

This is R codes for 3rd CMI_PB challenge, the prediction of pertussis vaccination response
We imputed missing data by softImpute, and MOFA algorithm to reduce dimensions, then Lasso regression to predict the response

###############
Procedure:
Features included: Antibody titres, cell frequencies, RNA sequences
softImpute imputation
MOFA: Reduced to 10 principle components
Lasso regression:
  Task 1.1 - IgG-PT Day 14: used 10 MOFA components + IgG-PT on Day 0, in total 11 factors
  Task 1.2 - IgG-PT Fold Change: calculated by results of task 1.1
  Task 2.1 - Monocytes Day 1: used 10 MOFA components + Monocytes on Day 0, in total 11 factors
  Task 2.2 - Monocytes Fold Change: calculated by results of task 2.1
  Task 3.1 - CCL3 Day 3: used 10 MOFA components, 10 factors
  Task 3.2 - CCL3 Fold Change: calculated by results of task 3.1

###############
R packages needed:
tidyverse: tables, data processing
softImpute: missing data imputation
mogsa: MOFA dimension reduction
glmnet: LASSO regression

###############
R scripts:
cmipb_data_clean.R - Data load and clean. Exclude features and subjects with high proportion of missing in training dataset. Imputation.
cmipb_challenge_model3.R - MOFA dimension reduction. Lasso regression prediction.

###############
RData files
data_cleaned.RData - Training and challenge data, with all subjects and features
demo_data.RData - Demographic data
train_challenge_cleaned.RData - Training and challenge data after exclusion of highly missed features and subjects
imputed_softImpute.RData - Imputed data

###############
3rdChallengeSubmissionTemplate_10032024(3).tsv - Submission template
cmipb_challenge_model3.tsv - Results to be submitted

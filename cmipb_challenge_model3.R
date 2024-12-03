library(tidyverse)
library(mogsa)
library(glmnet)

load("imputed_softImpute.RData")
load("data_cleaned.RData")

train_set <- train_softimp %>% select(-dataset)
challenge_set <- challenge_softimp %>% select(-dataset)

################################
#Task 1.1 & 1.2: IgG titre
igg_raw <- challenge_set %>% select(subject_id) %>%
  left_join(challenge_tbl %>% filter(timepoint == 0) %>% select(subject_id, IgG_PT_raw),
            by = "subject_id")

train_set1 <- train_set %>% left_join(train_y %>% select(subject_id, IgG_PT_raw),
                                      by = "subject_id")
train_set1 <- train_set1[!is.na(train_set1$IgG_PT_raw),]

mbpca_rst1 <- mbpca(list(t(as.matrix(train_set1[,5:35])), 
                         t(as.matrix(train_set1[,36:74])),
                         t(as.matrix(train_set1[,75:6677]))), ncomp = 10,
                    method = "globalScore", moa = FALSE)

set.seed(1922)
cv_model <- cv.glmnet(cbind(train_set1$IgG_PT,mbpca_rst1$t), 
                      train_set1$IgG_PT_raw, alpha = 1,
                      family='gaussian')

lasso_model <- glmnet(cbind(train_set1$IgG_PT,mbpca_rst1$t), 
                      train_set1$IgG_PT_raw, alpha = 1, family='gaussian', 
                      lambda = cv_model$lambda.min)

coef(lasso_model)

x_test1 <- as.matrix(challenge_set[,5:35]) %*% mbpca_rst1$pb$data_1
x_test2 <- as.matrix(challenge_set[,36:74]) %*% mbpca_rst1$pb$data_2
x_test3 <- as.matrix(challenge_set[,75:6677]) %*% mbpca_rst1$pb$data_3
x_test_t1 <- matrix(0, nrow = nrow(challenge_set), ncol = 10)
for(p in c(1:10)){
  x_test1[,p] <- x_test1[,p] - mean(x_test1[,p])
  x_test2[,p] <- x_test2[,p] - mean(x_test2[,p])
  x_test3[,p] <- x_test3[,p] - mean(x_test3[,p])
  x_test_t1[,p] <- x_test1[,p]*mbpca_rst1$w[1,p] +
    x_test2[,p]*mbpca_rst1$w[2,p] +
    x_test3[,p]*mbpca_rst1$w[3,p]
}

lasso_pred1 <- predict(lasso_model, s = cv_model$lambda.min, 
                       newx = cbind(challenge_set$IgG_PT, x_test_t1))

result_11 <- challenge_set %>% select(subject_id) %>%
  mutate(pred = lasso_pred1,
         rank = rank(-lasso_pred1))
result_11 <- result_11 %>% arrange(subject_id)

result_12 <- challenge_set %>% select(subject_id) %>%
  mutate(pred = lasso_pred1/igg_raw$IgG_PT_raw,
         rank = rank(-lasso_pred1/igg_raw$IgG_PT_raw))
result_12 <- result_12 %>% arrange(subject_id)

################################
#Task 2.1 & 2.2: Monocytes

train_set2 <- train_set %>% 
  left_join(train_y %>% select(subject_id, Monocytes) %>%
              rename(Monocytes_y = Monocytes),
            by = "subject_id")
train_set2 <- train_set2[!is.na(train_set2$Monocytes_y),]

mbpca_rst2 <- mbpca(list(t(as.matrix(train_set2[,5:35])), 
                         t(as.matrix(train_set2[,36:74])),
                         t(as.matrix(train_set2[,75:6677]))), ncomp = 10,
                    method = "globalScore", moa = FALSE)

set.seed(1922)
cv_model <- cv.glmnet(cbind(train_set2$cell_freq_1, mbpca_rst2$t), 
                      train_set2$Monocytes_y, alpha = 1,
                      family='gaussian')

lasso_model <- glmnet(cbind(train_set2$cell_freq_1, mbpca_rst2$t), 
                      train_set2$Monocytes_y, alpha = 1, family='gaussian', 
                      lambda = cv_model$lambda.min)

coef(lasso_model)

x_test1 <- as.matrix(challenge_set[,5:35]) %*% mbpca_rst2$pb$data_1
x_test2 <- as.matrix(challenge_set[,36:74]) %*% mbpca_rst2$pb$data_2
x_test3 <- as.matrix(challenge_set[,75:6677]) %*% mbpca_rst2$pb$data_3
x_test_t2 <- matrix(0, nrow = nrow(challenge_set), ncol = 10)
for(p in c(1:10)){
  x_test1[,p] <- x_test1[,p] - mean(x_test1[,p])
  x_test2[,p] <- x_test2[,p] - mean(x_test2[,p])
  x_test3[,p] <- x_test3[,p] - mean(x_test3[,p])
  x_test_t2[,p] <- x_test1[,p]*mbpca_rst2$w[1,p] +
    x_test2[,p]*mbpca_rst2$w[2,p] +
    x_test3[,p]*mbpca_rst2$w[3,p]
}

lasso_pred2 <- predict(lasso_model, s = cv_model$lambda.min, 
                       newx = cbind(challenge_set$cell_freq_1, x_test_t2))

result_21 <- challenge_set %>% select(subject_id) %>%
  mutate(pred = lasso_pred2,
         rank = rank(-lasso_pred2))
result_21 <- result_21 %>% arrange(subject_id)

result_22 <- challenge_set %>% select(subject_id) %>%
  mutate(pred = lasso_pred2/challenge_set$cell_freq_1,
         rank = rank(-lasso_pred2/challenge_set$cell_freq_1))
result_22 <- result_22 %>% arrange(subject_id)

################################
#Task 3.1 & 3.2: CCL

train_set3 <- train_set %>% 
  left_join(train_y %>% select(subject_id, tpm_ENSG00000277632.1) %>%
              rename(ENSG_y = tpm_ENSG00000277632.1),
            by = "subject_id")
train_set3 <- train_set3[!is.na(train_set3$ENSG_y),]
train_set3 <- train_set3 %>% mutate(ENSG_y = log(ENSG_y))


mbpca_rst3 <- mbpca(list(t(as.matrix(train_set3[,5:35])), 
                         t(as.matrix(train_set3[,36:74])),
                         t(as.matrix(train_set3[,75:6677]))), ncomp = 10,
                    method = "globalScore", moa = FALSE)

set.seed(1922)
cv_model <- cv.glmnet(mbpca_rst3$t, 
                      train_set3$ENSG_y, alpha = 1,
                      family='gaussian')

lasso_model <- glmnet(mbpca_rst3$t, 
                      train_set3$ENSG_y, alpha = 1, family='gaussian', 
                      lambda = cv_model$lambda.min)

coef(lasso_model)

x_test1 <- as.matrix(challenge_set[,5:35]) %*% mbpca_rst3$pb$data_1
x_test2 <- as.matrix(challenge_set[,36:74]) %*% mbpca_rst3$pb$data_2
x_test3 <- as.matrix(challenge_set[,75:6677]) %*% mbpca_rst3$pb$data_3
x_test_t3 <- matrix(0, nrow = nrow(challenge_set), ncol = 10)
for(p in c(1:10)){
  x_test1[,p] <- x_test1[,p] - mean(x_test1[,p])
  x_test2[,p] <- x_test2[,p] - mean(x_test2[,p])
  x_test3[,p] <- x_test3[,p] - mean(x_test3[,p])
  x_test_t3[,p] <- x_test1[,p]*mbpca_rst3$w[1,p] +
    x_test2[,p]*mbpca_rst3$w[2,p] +
    x_test3[,p]*mbpca_rst3$w[3,p]
}

lasso_pred3 <- predict(lasso_model, s = cv_model$lambda.min, newx = x_test_t3)

result_31 <- challenge_set %>% select(subject_id) %>%
  mutate(pred = lasso_pred3,
         rank = rank(-lasso_pred3))
result_31 <- result_31 %>% arrange(subject_id)

result_32 <- challenge_set %>% select(subject_id) %>%
  mutate(pred = lasso_pred3 - log(challenge_set$tpm_ENSG00000277632.1),
         rank = rank(-lasso_pred3 - log(challenge_set$tpm_ENSG00000277632.1)))
result_32 <- result_32 %>% arrange(subject_id)

################################
#Fill the result table

result_table <- read_tsv("3rdChallengeSubmissionTemplate_10032024 (3).tsv")

result_table$`1.1) IgG-PT-D14-titer-Rank` <- result_11$rank
result_table$`1.2) IgG-PT-D14-FC-Rank` <- result_12$rank
result_table$`2.1) Monocytes-D1-Rank` <- result_21$rank
result_table$`2.2) Monocytes-D1-FC-Rank` <- result_22$rank
result_table$`3.1) CCL3-D3-Rank` <- result_31$rank
result_table$`3.2) CCL3-D3-FC-Rank` <- result_32$rank

result_table %>% View()

write_tsv(result_table, file = "cmipb_challenge_model3.tsv")


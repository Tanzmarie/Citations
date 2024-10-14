library(tidyverse)
library(data.table)
library(mboost)

# load sjr data set
data = fread("data/cit_sjr.csv", sep = ",")


# transform variables accordingly
data[,c("References", "MeSH", "Title", "Length", "Age")] = lapply(data[,c("References", "MeSH", "Title", "Length", "Age")], asinh)
data[,c("Year", "Triangle")] = lapply(data[,c("Year", "Triangle")], as.factor)
data[,c("Language")] = lapply(data[,c("Language")], as.factor)

# Create Train and test data (10204457)
set.seed(123)

split_year = function(df) {
  smp_size = floor(0.8 * nrow(df))
  train_ind = sample(seq_len(nrow(df)), size = smp_size)
  train = df[train_ind, ]
  test = df[-train_ind, ]
  list(train = train, test = test)
}

split_data = data %>%
  group_by(Year) %>%
  group_split() %>%
  map(~ split_year(.x))

train = bind_rows(map(split_data, "train"))
test = bind_rows(map(split_data, "test"))

rm(data, split_year, split_data)

# Boost the models
set.seed(123)
B50 = glmboost(asinh(SJR) ~ ., 
                control = boost_control(mstop = 50, nu = 0.1, trace = TRUE),
                data = train)


B100 = glmboost(asinh(SJR) ~ ., 
                control = boost_control(mstop = 100, nu = 0.1, trace = TRUE),
                data = train)

B250 = glmboost(asinh(SJR) ~ ., 
                 control = boost_control(mstop = 250, nu = 0.1, trace = TRUE),
                 data = train)

# Optimal Stopping criterion
# n.train = round(0.9 * nrow(train))
# n.val = round(0.1 * nrow(train))
# 
# weight.mstop = c(rep(1, times = n.train),rep(0, times = n.val))
# weight.mstop = sample(weight.mstop)
# 
# 
# BOPT = glmboost(asinh(SJR) ~ .,
#                control =  boost_control(trace = TRUE, mstop = 5000, risk = "oobag", nu = 0.1),
#                weights = weight.mstop,
#                data = train)
# 
# StopIT = which.min(risk(BOPT,merge = T))
# 
# train_opt = train[weight.mstop == 1, ]
# 
# BOPT = glmboost(asinh(SJR) ~ .,
#                control =  boost_control(trace = TRUE, mstop = StopIT, nu = 0.1),
#                data = train_opt)

# calculate performance
options(scipen = 999)

fm = list("B50" = B50, "B100 " = B100 , "B250" = B250)

perf_train = rbind(NLL = sapply(fm, function(x) round(-x$logLik(), digits = 0)),
                   R2 = sapply(fm, function(x) {
                     y_pred = x$fitted()
                     y_true = asinh(train$SJR)
                     
                   R2 = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)}),
                   MSEP = sapply(fm, function(x) mean((asinh(train$SJR) - fitted(x, type = "response"))^2)),
                   MAE = sapply(fm, function(x) mean(abs(asinh(train$SJR) - fitted(x, type = "response")))))



perf_test = rbind(NLL = sapply(fm, function(x) {
                  sum((asinh(test$SJR) - predict(x, newdata = test, type = "response"))^2)}),
                  MSEP = sapply(fm, function(x) mean((asinh(test$SJR) - predict(x, newdata = test, type = "response"))^2)),
                  MAE = sapply(fm, function(x) mean(abs(asinh(test$SJR) - predict(x, newdata = test, type = "response")))))

# Get summary of coef and selection frequency
coefs = sapply(fm, function(x) unlist(x$coef()))
sels = sapply(fm, function(x) summary(x)$selprob)




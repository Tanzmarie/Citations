library(tidyverse)
library(data.table)
library(mboost)

# load data set 
chunk = fread("data/cit.csv", sep = ",")

# transform variables accordingly
chunk[,c("References", "MeSH", "Title", "Length", "Age")] = lapply(chunk[,c("References", "MeSH", "Title", "Length", "Age")], asinh)
chunk[,c("Year", "Triangle")] = lapply(chunk[,c("Year", "Triangle")], as.factor)
chunk[,c("Language")] = lapply(chunk[,c("Language")], as.factor)

# Create Train and test data
set.seed(123)

split_year = function(df) {
  smp_size = floor(0.8 * nrow(df))
  train_ind = sample(seq_len(nrow(df)), size = smp_size)
  train = df[train_ind, ]
  test = df[-train_ind, ]
  list(train = train, test = test)
}

split_data = chunk %>%
  group_by(Year) %>%
  group_split() %>%
  map(~ split_year(.x))

train = bind_rows(map(split_data, "train"))
test = bind_rows(map(split_data, "test"))

rm(chunk, split_year, split_data)

# Gradient Boosting (Negative Binomial)
set.seed(123)
B50 = glmboost(Citations ~ ., family = NBinomial(),
                   control = boost_control(mstop = 50, nu = 0.1, trace = TRUE),
                   data = train)

B100 = glmboost(Citations ~ ., family = NBinomial(),
                control = boost_control(mstop = 100, nu = 0.1, trace = TRUE),
                data = train)

B250 = glmboost(Citations ~ ., family = NBinomial(),
                 control = boost_control(mstop = 250, nu = 0.1, trace = TRUE),
                 data = train)

# Optimal Stopping criterion
# n.train = round(0.9 * nrow(train))
# n.val = round(0.1 * nrow(train))
# 
# weight.mstop = c(rep(1, times = n.train),rep(0, times = n.val))
# weight.mstop = sample(weight.mstop)
# 
# BOPT = glmboost(Citations ~ ., family = NBinomial(),
#                    control =  boost_control(trace = TRUE, mstop = 5000, risk = "oobag", nu = 0.1),
#                    weights = weight.mstop,
#                    data = train)
# 
# StopIT = which.min(risk(BOPT,merge = T))
# 
# train_opt = train[weight.mstop == 1, ]
# 
# BOPT = glmboost(asinh(SJR) ~ .,
#                 control =  boost_control(trace = TRUE, mstop = StopIT, nu = 0.1),
#                 data = train_opt)

# calculate performance
options(scipen = 999)

fm = list( "B50" = B50, "B100 " = B100 , "B250" = B250)

perf_train = rbind(NLL = sapply(fm, function(x) round(-x$logLik(), digits = 0)),
                   MSEP = sapply(fm, function(x) mean((train$Citations - fitted(x, type = "response"))^2)),
                   MAE = sapply(fm, function(x) mean(abs(train$Citations - fitted(x, type = "response")))))

perf_test = rbind(NLL = sapply(fm, function(x) {
                  nu = predict(x, newdata = test, type = "response")
                  psi = x$nuisance()[[50]]
                  cit = test$Citations
  
                  logLik_nb = sum(-(lgamma(cit + psi) - lgamma(psi) - lgamma(cit + 1) +
                      psi * log(psi) - psi*log(nu + psi) + cit * log(nu) -
                      cit * log(nu + psi)))
  
                  round(logLik_nb, digits = 0)}),
                  MSEP = sapply(fm, function(x) mean((test$Citations - predict(x, newdata = test, type = "response"))^2)),
                  MAE = sapply(fm, function(x) mean(abs(test$Citations - predict(x, newdata = test, type = "response"))))
)

# Get summary of coef and selection frequency
coefs = sapply(fm, function(x) unlist(x$coef()))
sels = sapply(fm, function(x) summary(x)$selprob)





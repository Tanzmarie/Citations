library(tidyverse)
library(data.table)
library(mboost)

# load sjr data set
data = fread("./cit_sjr.csv", sep = ",")

# transform variables accordingly
data[,c("References", "MeSH", "Title", "Length", "Age")] = lapply(data[,c("References", "MeSH", "Title", "Length", "Age")], asinh)
data[,c("Year", "Triangle")] = lapply(data[,c("Year", "Triangle")], as.factor)
data[,c("Language")] = lapply(data[,c("Language")], as.factor)

# Create Train and test data (10204457)

split_year = function(df, p) {
  smp_size = floor(p * nrow(df))
  train_ind = sample(seq_len(nrow(df)), size = smp_size)
  train = df[train_ind, ]
  test = df[-train_ind, ]
  list(train = train, test = test)
}

set.seed(123)

split_data = data %>%
  group_by(Year) %>%
  group_split() %>%
  map(~ split_year(.x, p = 0.8))

train = bind_rows(map(split_data, "train"))
test = bind_rows(map(split_data, "test"))

# Create train_val and validation set
set.seed(123)

split_data2 = train %>%
  group_by(Year) %>%
  group_split() %>%
  map(~ split_year(.x, p = 0.99))

train_val = bind_rows(map(split_data2, "train"))
val = bind_rows(map(split_data2, "test"))

rm(data, split_data, split_data2)
gc()

# Function to repeat subsampling and glmboost
repeat_glmboost = function(iterations, train_val, val) {
  stopit_vals = vector("numeric", length = iterations)
  
  for (i in 1:iterations) {
    # Subsample
    set.seed(i) 
    split_data3 = train_val %>%
      group_by(Year) %>%
      group_split() %>%
      map(~ split_year(.x, p = 0.01))
    
    train_sub = bind_rows(map(split_data3, "train"))
    opt_val = bind_rows(train_sub, val)
    
    # Create weights
    weight.mstop = c(rep(1, times = nrow(train_sub)), rep(0, times = nrow(val)))
    
    # Train the model with glmboost
    BOPT = glmboost(
      asinh(SJR) ~ .,
      control = boost_control(trace = TRUE, mstop = 3000, risk = "oobag", nu = 0.1),
      weights = weight.mstop,
      data = opt_val
    )
    
    # Get optimal stopping iteration
    stopit_vals[i] = which.min(risk(BOPT, merge = TRUE))
    
    # Clean up
    rm(split_data3, train_sub, opt_val, weight.mstop, BOPT)
    gc(verbose = TRUE)
  }
  
  return(stopit_vals)
}

# Repeat for 10 times
set.seed(123)
stopit_vals = repeat_glmboost(iterations = 10, train_val = train_val, val = val)

# Print the StopIT values
print(stopit_vals)

StopIT = mean(stopit_vals)

rm(stopit_vals, train_val, val, split_year, repeat_glmboost)

set.seed(123)
FOPT = glmboost(asinh(SJR) ~ .,
                control =  boost_control(trace = TRUE, mstop = StopIT, nu = 0.1),
                data = train)

# calculate performance
options(scipen = 999)

fm = list("FOPT" = FOPT)

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

# Print results

print(StopIT)
print(perf_train)
print(perf_test)
print(coefs)
print(sels)


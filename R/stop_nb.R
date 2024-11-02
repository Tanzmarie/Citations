library(tidyverse)
library(data.table)
library(mboost)

# load sjr chunk set
chunk = fread("./cit.csv", sep = ",")

# transform variables accordingly
chunk[,c("References", "MeSH", "Title", "Length", "Age")] = lapply(chunk[,c("References", "MeSH", "Title", "Length", "Age")], asinh)
chunk[,c("Year", "Triangle")] = lapply(chunk[,c("Year", "Triangle")], as.factor)
chunk[,c("Language")] = lapply(chunk[,c("Language")], as.factor)

# Create Train and test chunk (10204457)

split_year = function(df, p) {
  smp_size = floor(p * nrow(df))
  train_ind = sample(seq_len(nrow(df)), size = smp_size)
  train = df[train_ind, ]
  test = df[-train_ind, ]
  list(train = train, test = test)
}

set.seed(123)

split_chunk = chunk %>%
  group_by(Year) %>%
  group_split() %>%
  map(~ split_year(.x, p = 0.8))

train = bind_rows(map(split_chunk, "train"))
test = bind_rows(map(split_chunk, "test"))

# Create train_val and validation set
set.seed(123)

split_chunk2 = train %>%
  group_by(Year) %>%
  group_split() %>%
  map(~ split_year(.x, p = 0.995))

train_val = bind_rows(map(split_chunk2, "train"))
val = bind_rows(map(split_chunk2, "test"))

rm(chunk, split_chunk, split_chunk2)
gc()

# Function to repeat subsampling and glmboost
repeat_glmboost = function(iterations, train_val, val) {
  stopit_vals = vector("numeric", length = iterations)
  
  for (i in 1:iterations) {
    # Subsample
    set.seed(i) 
    split_chunk3 = train_val %>%
      group_by(Year) %>%
      group_split() %>%
      map(~ split_year(.x, p = 0.02))
    
    train_sub = bind_rows(map(split_chunk3, "train"))
    opt_val = bind_rows(train_sub, val)
    
    # Create weights
    weight.mstop = c(rep(1, times = nrow(train_sub)), rep(0, times = nrow(val)))
    
    # Train the model with glmboost
    BOPT = glmboost(Citations ~ ., family = NBinomial(nuirange = c(0,100)),
                    control =  boost_control(trace = TRUE, mstop = 1000, risk = "oobag", nu = 0.1),
                    weights = weight.mstop,
                    data = opt_val)
    
    # Get optimal stopping iteration
    stopit_vals[i] = which.min(risk(BOPT, merge = T))
    
    # Clean up
    rm(split_chunk3, train_sub, opt_val, weight.mstop, BOPT)
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

rm(stopit_vals, val, train_val, split_year)
gc()

set.seed(123)
FOPT = glmboost(Citations ~ ., family = NBinomial(nuirange = c(0,100)),
                control =  boost_control(trace = TRUE, mstop = StopIT, nu = 0.1),
                data = train)

# calculate performance
options(scipen = 999)

fm = list("FOPT" = FOPT)

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

# Print results

print(StopIT)
print(perf_train)
print(perf_test)
print(coefs)
print(sels)

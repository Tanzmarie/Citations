library(tidyverse)
library(data.table)

# load sjr data set
data = fread("data/cit_sjr.csv", sep = ",")

# transform variables accordingly
data[,c("References", "MeSH", "Title", "Length", "Age")] = lapply(data[,c("References", "MeSH", "Title", "Length", "Age")], asinh)
data[,c("Year", "Triangle")] = lapply(data[,c("Year", "Triangle")], as.factor)
data[,c("Language")] = lapply(data[,c("Language")], as.factor)

# Create Train and test data
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


# Linear Model with continuous sjr values
mod_lmc = lm(asinh(SJR) ~ ., data = train)

summary(mod_lmc)

# Only meta_data variables
mod_lmi = lm(asinh(SJR) ~ ., data = train[,c(1:15)])

summary(mod_lmi)

# Only cont variables
mod_lmr = lm(asinh(SJR) ~ ., data = train[,c("SJR","H.Score","A.Score", "C.Score","Title","References","Age","MeSH","Length")])

summary(mod_lmr)

# Performance metrics of test and train data
options(scipen = 999)

fm = list("LMC" = mod_lmc,"LMI" = mod_lmi, "LMR" = mod_lmr)

# Train set
perf_train = rbind(NLL = sapply(fm, function(x) round(- logLik(x), digits = 0)),
      Df = sapply(fm, function(x) attr(logLik(x), "df")),
      AIC = sapply(fm, function(x) round(AIC(x))),
      BIC = sapply(fm, function(x) round(BIC(x))),
      MSEP = sapply(fm, function(x) mean((asinh(train$SJR) - fitted(x))^2)),
      MAE = sapply(fm, function(x) mean(abs(asinh(train$SJR) - fitted(x)))),
      R2 = sapply(fm, function(x) summary(x)$r.squared),
      Adj.R2 = sapply(fm, function(x) summary(x)$adj.r.squared))

# Test set
perf_test = rbind(NLL = sapply(fm, function(x) {
  n = length(test$SJR)
  sigma = summary(x)$sigma
  sigma.ML = sigma*sqrt((n-dim(model.matrix(x))[2])/n) 
  
  - sum(log(dnorm(x = asinh(test$SJR), mean = predict(x, newdata = test), sd = sigma.ML)))
}),
MSEP = sapply(fm, function(x) mean((asinh(test$SJR) - predict(x, newdata = test))^2)),
MAE = sapply(fm, function(x) mean(abs(asinh(test$SJR) - predict(x, newdata = test))))
)



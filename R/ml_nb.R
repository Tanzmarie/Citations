library(tidyverse)
library(data.table)
library(pscl)
library(MASS)

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

# Count data regression  for citation count (adjust model accordingly)
mod_glmc = glm.nb(Citations ~ ., trace = TRUE, data = train)

summary(mod_glmc)

mod_glmi = glm.nb(Citations ~ ., trace = TRUE, data = train[,1:15])

summary(mod_glmi)


mod_glmr = glm.nb(Citations ~ ., trace = TRUE, data = train[,c("Citations","H.Score","A.Score", "C.Score","Title","References","Age","MeSH","Length")])
summary(mod_glmr)

stargazer(mod_glmc, mod_glmi, mod_glmr)

# Performance metrics of test and train data
options(scipen = 999)

fm = list("GLMC" = mod_glmc,"GLMI" = mod_glmi, "GLMR" = mod_glmr)

# Train set
perf_train = rbind(NLL = sapply(fm, function(x) round(-logLik(x), digits = 0)),
      Df = sapply(fm, function(x) attr(logLik(x), "df")),
      AIC = sapply(fm, function(x) round(AIC(x))),
      BIC = sapply(fm, function(x) round(BIC(x))),
      MSEP = sapply(fm, function(x) mean((train$Citations - fitted(x, type = "response"))^2)),
      MAE = sapply(fm, function(x) mean(abs(train$Citations - fitted(x, type = "response")))))

# Test set
perf_test = rbind(NLL = sapply(fm, function(x) {
  nu = predict(x, newdata = test, type = "response")
  psi = x$theta
  cit = test$Citations
  
  logLik_nb = sum(-(lgamma(cit + psi) - lgamma(psi) - lgamma(cit + 1) +
                         psi * log(psi) - psi*log(nu + psi) + cit * log(nu) -
                         cit * log(nu + psi)))
  
  round(logLik_nb, digits = 0)
}),
MSEP = sapply(fm, function(x) mean((test$Citations - predict(x, newdata = test, type = "response"))^2)),
MAE = sapply(fm, function(x) mean(abs(test$Citations - predict(x, newdata = test, type = "response"))))
)



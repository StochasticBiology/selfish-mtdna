# library for some logistic regression summary stats
library(pscl)
library(arm)

# read in data
x<-read.table("kang-competition.txt", header=T)

# construct different (Bayesian) logistic regression models 
model4 <- bayesglm(Did_A_win ~ l_300_diff + l_16345_diff + l_16287_diff + l_16180_diff, family=binomial, data=x)

model3a <- bayesglm(Did_A_win ~ l_300_diff + l_16345_diff + l_16180_diff, family=binomial, data=x)
model3b <- bayesglm(Did_A_win ~ l_300_diff + l_16287_diff + l_16180_diff, family=binomial, data=x)
model3c <- bayesglm(Did_A_win ~ l_300_diff + l_16287_diff + l_16345_diff, family=binomial, data=x)

model2a <- bayesglm(Did_A_win ~ l_300_diff + l_16345_diff, family=binomial, data=x)
model2b <- bayesglm(Did_A_win ~ l_300_diff + l_16180_diff, family=binomial, data=x)
model2c <- bayesglm(Did_A_win ~ l_300_diff + l_16287_diff, family=binomial, data=x)
model1 <- bayesglm(Did_A_win ~ l_300_diff, family=binomial, data=x)

# examine different models
c(BIC(model4), BIC(model3a), BIC(model3b), BIC(model3c), BIC(model2a), BIC(model2b), BIC(model2c), BIC(model1)) 

# now use traditional logistic regression function for comparing chosen and null models
model2a <- glm(Did_A_win ~ l_300_diff + l_16345_diff, family=binomial, data=x)
model1 <- glm(Did_A_win ~ l_300_diff, family=binomial, data=x)

pR2(model2a)
pR2(model1)
summary(model2a)

# select best model
model = model2a

# overall performance
fitted.results <- predict(model, newdata=x, type="response")
fitted.results <- ifelse(fitted.results > 0.5, 1, 0)
fitted.results

# repeatedly sample half of reverting and half of non-reverting observations to train/test model
winners = c(1:4)
losers = c(5:26)
errors = c()
for(i in 1:1000)
{
  # sample random half of each observation type
  winc = sample(winners, length(winners)/2)
  losec = sample(losers, length(losers)/2)

  # construct training and test sets
  trainx = x[c(winc,losec),]
  testx = x[c(setdiff(winners, winc), setdiff(losers, losec)),]

  # fit on training set and predict test set
  fitted.results <- predict(model, newdata=testx, type="response")
  fitted.results <- ifelse(fitted.results > 0.5, 1, 0)
  error <- mean(fitted.results != testx$Did_A_win)
  errors = c(errors, error)
}

# report mean performance
1-mean(errors)
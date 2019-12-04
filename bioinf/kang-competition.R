# library for some logistic regression summary stats
library(pscl)

# read in data
x<-read.table("kang-competition.txt", header=T)

# construct different logistic regression models
model2a <- glm(Did_A_win ~ l_300_diff + l_16345_diff, family=binomial(link='logit'), data=x)
model2b <- glm(Did_A_win ~ l_300_diff + l_16180_diff, family=binomial(link='logit'), data=x)
model1 <- glm(Did_A_win ~ l_300_diff, family=binomial(link='logit'), data=x)

# examine different models
c(AIC(model2a), AIC(model2b), AIC(model1))
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
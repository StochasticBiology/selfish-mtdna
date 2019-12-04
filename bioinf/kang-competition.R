x<-read.table("kang-competition.txt", header=T)
model2a <- glm(Did_A_win ~ l_300_diff + l_16345_diff, family=binomial(link='logit'), data=x)
model2b <- glm(Did_A_win ~ l_300_diff + l_16180_diff, family=binomial(link='logit'), data=x)
model1 <- glm(Did_A_win ~ l_300_diff, family=binomial(link='logit'), data=x)
library(pscl)
c(AIC(model2a), AIC(model2b), AIC(model1))

pR2(model2a)
pR2(model1)
model = model2a
fitted.results <- predict(model, newdata=x, type="response")
fitted.results <- ifelse(fitted.results > 0.5, 1, 0)
fitted.results
winners = c(1:4)
losers = c(5:26)
errors = c()
for(i in 1:1000)
{
  winc = sample(winners, length(winners)/2)
  losec = sample(losers, length(losers)/2)
  trainx = x[c(winc,losec),]
  testx = x[c(setdiff(winners, winc), setdiff(losers, losec)),]
  fitted.results <- predict(model, newdata=testx, type="response")
  fitted.results <- ifelse(fitted.results > 0.5, 1, 0)
  error <- mean(fitted.results != testx$Did_A_win)
  errors = c(errors, error)
}
mean(errors)
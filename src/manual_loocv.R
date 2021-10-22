#First initialize the output vector as an empty object outside the loop.
dt <- mets_all2[-c(6,7)]

fitted_value <- NULL
for(i in 1:29){
  #you did this part right
  validation<-dt[i,]
  training<-dt[-29,]
  model1<-lm(form.vtigl, data = training)
  #when you fit the model, use the newdata argument to predict on a new row
  #also, fitted_value needs the index [i], so the each loop doesn't overwrite the previous
  fitted_value[i] <- predict(model1, newdata = validation)}
rm(list = ls()) #cleaning the workspace before starting analaysis

#formatting data
data<-read.csv('mcf7.csv',header = TRUE) 
rownames(data)<-data[,1]
data<-data[,-1]

testdata<-t(data) 
td<-as.data.frame(testdata)

library(janitor)
td<-remove_constant(td)
controls<-read.csv('mcf7_controls.csv',header = TRUE) 
td$label<-factor(controls$label,labels=c("Yes","No"))
sum(td$label=="Yes")

nr<-50 #noise rate (corresponds to max-min)

library(cwhmisc)
library(freqparcoord)

td_noise<-matrix(nrow = dim(td)[1], ncol = dim(td)[2]-1)


for(i in 1:(dim(td)[2]-1)){
  td_noise[,i]<-jitter(td[,i],factor=nr, amount = 0)
}  
td_noise<-as.data.frame(td_noise)
td_noise$label<-td$label
colnames(td_noise)<-colnames(td)
rownames(td_noise)<-rownames(td)

library(caTools)
set.seed(789)
split = sample.split(td$label, SplitRatio = 0.75)
td_new = subset(td, split == TRUE)
td_validation = subset(td, split == FALSE)
td_new_noise = subset(td_noise, split == TRUE)
td_validation_noise = subset(td_noise, split == FALSE)

td_svm<-as.data.frame(td_new)
td_svm2 <- as.data.frame(td_new_noise)

td_pred<-as.data.frame(td_validation) #noise
td_pred_noise<-as.data.frame(td_validation_noise) # no noise

library(caret)
library(tidyverse) 

set.seed(234)
trctrl <- trainControl(method = "cv", number = 10, savePredictions=TRUE, classProbs=TRUE,
                       summaryFunction=twoClassSummary)
svm_fit <- train(label ~., data = td_svm, method = "svmLinear", trControl=trctrl,
                 tuneLength = 0)
svm_fit

vr<-varImp(svm_fit, scale = FALSE)
#ggplot(vr)
vr2<-vr$importance
write.csv(vr2,'mcf7_no_noise.csv')

pred <- svm_fit$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)
eachfold <- pred %>%                                        
  group_by(Resample) %>%                         
  summarise_at(vars(equal),                     
               list(Accuracy = mean))              
eachfold
mean(eachfold$Accuracy)

ggplot(data=eachfold, aes(x=Resample, y=Accuracy, group=1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()

#validation
pred_data <- td_pred
td_last<-ncol(td)
y_pred = predict(svm_fit, newdata = pred_data[,-td_last])
cm = table(pred_data[, td_last], y_pred)
accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])



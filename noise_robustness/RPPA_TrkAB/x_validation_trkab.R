rm(list = ls()) #cleaning the workspace before starting analysis

#formatting

data<-read.csv('rppa_controls.csv',header = TRUE) 
#data<-read.csv('rppa_controls_red.csv',header = TRUE) 

rownames(data)<-data[,1]
data<-data[,-1]

library(varhandle)
tmp<-data[-nrow(data),]
testdata<-t(unfactor(tmp)) 
td<-as.data.frame(testdata)

library(janitor)
td<-remove_constant(td)

td$label<-factor(t(data[nrow(data),]),labels=c("Yes","No"))
table(td$label)

nr<-0.95 #noise rate 

library(cwhmisc)
library(freqparcoord)

td_new<-td

td_noise<-matrix(nrow = dim(td)[1], ncol = dim(td)[2]-1)

set.seed(123)

for(i in 1:(dim(td)[2]-1)){
  for(j in 1:(dim(td)[1])){
    if(abs(td[j,i])<1e-16){td_noise[j,i] = runif(1,0,nr*td[j,i])} 
    else if(abs(td[j,i]-1)>1e-16 & is.na(td[j,i])==FALSE){td_noise[j,i]=td[j,i]+runif(1,-nr*td[j,i],nr*td[j,i])}
    else{td_noise[j,i]=td_new[j,i]} 
  }
}


td_svm<-as.data.frame(td_new) #clean
td_svm2 <- as.data.frame(td_noise) #noise
td_svm2$label<-td_svm$label
colnames(td_svm2)<-colnames(td_svm)
rownames(td_svm2)<-rownames(td_svm)



table(td_svm$label)
table(td_svm2$label)


library(caret)
library(tidyverse) 

set.seed(678)
trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, savePredictions=TRUE, classProbs=TRUE,
                       summaryFunction=twoClassSummary)
svm_fit <- train(label ~., data = td_svm2, method = "svmLinear", trControl=trctrl,
                 tuneLength = 0)
svm_fit

vr<-varImp(svm_fit, scale = FALSE)
#ggplot(vr)
vr2<-vr$importance
write.csv(vr2,'trkab_noise.csv')


pred <- svm_fit$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)
eachfold <- pred %>%                                        
  group_by(Resample) %>%                         
  summarise_at(vars(equal),                     
               list(Accuracy = mean))              
eachfold
mean(eachfold$Accuracy)

#plot
ggplot(data=eachfold, aes(x=Resample, y=Accuracy, group=1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()




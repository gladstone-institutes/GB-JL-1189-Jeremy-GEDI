#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/
#https://atmos.eoas.fsu.edu/~parker/tr504.pdf
library(mclust)
library(ggplot2)
library(dplyr)

setwd("~/Dropbox (Gladstone)/GB-JL-1189-Jeremy-GEDI")


MN1<-read.csv("TimelapseDataset1_MNs.csv")
MN1$Timepoint=as.factor(MN1$Timepoint)


MN1$Well=as.factor(MN1$Well)


MN2<-read.csv("TimelapseDataset2_MNs.csv")
MN2$Timepoint=as.factor(MN2$Timepoint)


MN2$Well=as.factor(MN2$Well)
table(MN2$Timepoint)
View(subset(MN2,Timepoint==1))

faithfulBIC2 <- mclustBIC(subset(MN2,Timepoint==1)$GEDIratio)
faithfulBIC2

#Top 3 models based on the BIC criterion: 
#  V,4      V,5      V,3 
#23128.92 23119.94 23094.00 

model2=Mclust(subset(MN2,Timepoint==1)$GEDIratio, x = faithfulBIC2)
summary(model2)     
plot(model2)

ICL <- mclustICL(subset(MN2,Timepoint==1)$GEDIratio)
ICL


drmod <- MclustDR(model2, lambda = 1)
plot(drmod, what = 'contour')



model_two=Mclust(subset(MN2,Timepoint==1)$GEDIratio,G=2)
summary(model_two)     
plot(model_two)


#data1 two cluster cases
i=1 #4 8 10 11 12 13
model_two=Mclust(subset(MN1,Timepoint==i)$GEDIratio,G=2,initialization = list(hcPairs = hcRandomPairs(subset(MN1,Timepoint==i)$GEDIratio, seed = 123)))
summary(model_two)     
plot(model_two)


#data2 two cluster cases
i=5#6 8 9 10 13
model_two=Mclust(subset(MN2,Timepoint==i)$GEDIratio,G=2,initialization = list(hcPairs = hcRandomPairs(subset(MN2,Timepoint==i)$GEDIratio, seed = 123)))
summary(model_two)     
plot(model_two)


dens <- densityMclust(subset(MN2,Timepoint==1)$GEDIratio)
br <- seq(min(subset(MN2,Timepoint==1)$GEDIratio), max(subset(MN2,Timepoint==1)$GEDIratio), length = 21)
plot(dens, what = "density", data = subset(MN2,Timepoint==1)$GEDIratio, breaks = br)


###1. using annotated data test MclustDA
###if it works fine, just use it 
annotated<-read.csv("/Users/mingyoungshin/Dropbox (Gladstone)/GB-JL-1189-Jeremy-GEDI/Curation_and_GEDI.csv")
dim(annotated)


set.seed(123)
train <- sample(1:nrow(annotated), size = round(nrow(annotated)*2/3), replace = FALSE)
annotated.train <- annotated[train,'GEDIratio']
Class.train <- annotated[train,c("Curator1_Accurate","Curator2_Accurate","Curator3_Accurate","Curator4_Accurate")]
Class.train[is.na(Class.train)]=0.5
Class.train=round(rowMeans(Class.train))
table(Class.train)


annotated.test <- annotated[-train,'GEDIratio']
Class.test <- annotated[-train,c("Curator1_Accurate","Curator2_Accurate","Curator3_Accurate","Curator4_Accurate")]
Class.test=round(rowMeans(Class.test))
table(Class.test)
 
mod1 <- MclustDA(annotated.train, as.vector((Class.train)), modelType = "EDDA")
 
cv <- cvMclustDA(mod1)
unlist(cv[c('error', 'se')])
#0.267097153 0.007610695 

model_annotated=Mclust(annotated.train,G=2,initialization = list(hcPairs = hcRandomPairs(annotated.train, seed = 123)))
summary(model_annotated)


qqnorm(annotated.train)
qqline(annotated.train)
hist(annotated.train,breaks=1000)
plot.ecdf(annotated.train)


#################Hampel filter
library(EnvStats)
library(DescTools)


upper_bound <- median(sqrt(MN2$GEDIratio)) + 3 * mad(sqrt(MN2$GEDIratio))
upper_bound

hist(sqrt(MN2$GEDIratio),breaks=1000,main=NULL)
abline(v=upper_bound, col="blue")

outlierC=sum(sqrt(MN2$GEDIratio)>upper_bound)
outlierC


test <- rosnerTest(sqrt(MN2$GEDIratio),
                   k = outlierC)
cutoff=min(test$all.stats$Value[test$all.stats$Outlier]) 

abline(v=cutoff,col="red")
mtext(paste0("MN2 combined :cutoff at ",cutoff),cex=1.2)

for(t in 1:max(as.numeric(MN2$Timepoint))){
  par(new=FALSE)
  
  upper_bound <- median(sqrt(subset(MN2,Timepoint==t)$GEDIratio)) + 3 * mad(sqrt(subset(MN2,Timepoint==t)$GEDIratio))
  upper_bound
  
  hist(sqrt(subset(MN2,Timepoint==t)$GEDIratio),breaks=1000,main=NULL)
  abline(v=upper_bound, col="blue")
  
  outlierC=sum(sqrt(subset(MN2,Timepoint==t)$GEDIratio)>upper_bound)
  outlierC
  

  test <- rosnerTest(sqrt(subset(MN2,Timepoint==t)$GEDIratio),
                     k = outlierC)
  cutoff=min(test$all.stats$Value[test$all.stats$Outlier]) 
  
  abline(v=cutoff,col="red")
  mtext(paste0("MN2 Time ",t," :cutoff at ",cutoff),cex=1.2)
  
}



upper_bound <- median(sqrt(MN1$GEDIratio)) + 3 * mad(sqrt(MN1$GEDIratio))
upper_bound

hist(sqrt(MN1$GEDIratio),breaks=1000, main=NULL)
abline(v=upper_bound, col="blue")

outlierC=sum(sqrt(MN1$GEDIratio)>upper_bound)
outlierC


test <- rosnerTest(sqrt(MN1$GEDIratio),
                   k = outlierC)
cutoff=min(test$all.stats$Value[test$all.stats$Outlier]) 

abline(v=cutoff,col="red")

mtext(paste0("MN1 combined :cutoff at ",cutoff),cex=1.2)



for(t in 1:max(as.numeric(MN1$Timepoint))){
  par(new=FALSE)
  
  upper_bound <- median(sqrt(subset(MN1,Timepoint==t)$GEDIratio)) + 3 * mad(sqrt(subset(MN1,Timepoint==t)$GEDIratio))
  upper_bound
  
  hist(sqrt(subset(MN1,Timepoint==t)$GEDIratio),breaks=1000,main=NULL)
  abline(v=upper_bound, col="blue")
  
  outlierC=sum(sqrt(subset(MN1,Timepoint==t)$GEDIratio)>upper_bound)
  outlierC
  
  test <- rosnerTest(sqrt(subset(MN1,Timepoint==t)$GEDIratio),
                     k = outlierC)
  outliers=test$all.stats$Value[test$all.stats$Outlier]
  cutoff=min( outliers[which(outliers>median(sqrt(subset(MN1,Timepoint==t)$GEDIratio)))] )
  
  
  abline(v=cutoff,col="red")
  mtext(paste0("MN1 Time ",t," :cutoff at ",cutoff),cex=1.2)
  
}


################by distribution
library(DescTools)

hist(MN2$GEDIratio,breaks=1000)

Mode(MN2$GEDIratio)
0.013124935
left=subset(MN2,GEDIratio<=0.013124935)

1) unknown σ 
68% of observations falls within the first standard deviation (µ ± σ)
=> 0.013124935 + σ = 34% of the left side data = quantile(left$GEDIratio, probs = seq(0.16)) = 0.01312466
σ=0.01312466-0.013124935

2) known σ
estimate sd of the left side
sd(left$GEDIratio)
0.001449059
=>estimate σ
left side + right side distribution
sqrt(var(2X))=2sqrt(var(X))
σ=0.001449059*2=0.002898118
=> 95% location
qnorm(0.95,mean=0.013124935,sd=0.002898118)  
=> if we decide to ignore 5% extreme values, use this as the cutoff
  
25% quantile median
qnorm(0.25,mean,sd)
0.67*sd=median
get sd
get whole
flexmixed


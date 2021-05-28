
setwd("~/Dropbox (Gladstone)/GB-JL-1189-Jeremy-GEDI")

library(ggplot2)
MN1<-read.csv("TimelapseDataset1_MNs.csv")
MN1$Timepoint=as.factor(MN1$Timepoint)

bp <- ggplot(MN1, aes(x=Timepoint, y=GEDIratio, fill=Timepoint)) + 
  geom_violin(trim=FALSE) + scale_fill_brewer(palette="Blues") + theme_classic()

jpeg("boxplot_time_GEDI_MN1.jpeg")
print(bp)
dev.off()

MN1$Well=as.factor(MN1$Well)

unique_wells=as.character(unique(MN1$Well))
unique_wells=c(unique_wells,rep(NA,5))
for(i in 1:10){
  subdata=MN1[which(as.character(MN1$Well)%in%unique_wells[10*(i-1)+1:10*i]),]
  
  bp <- ggplot(subdata, aes(x=Well, y=GEDIratio, fill=Well)) + 
    geom_violin(trim=FALSE) + theme_classic()
  
  jpeg(paste0("boxplot_well_GEDI_MN1_subset",i,".jpeg"))
  print(bp)
  dev.off()
}


bp <- ggplot(MN1, aes(x=Well, y=GEDIratio, fill=Well)) + 
  geom_boxplot()+ theme_classic()

jpeg("boxplot_well_GEDI_MN1.jpeg")
print(bp)
dev.off()





MN2<-read.csv("TimelapseDataset2_MNs.csv")
MN2$Timepoint=as.factor(MN2$Timepoint)


unique_wells=as.character(unique(MN2$Well))
unique_wells=c(unique_wells,rep(NA,4))
for(i in 1:10){
  subdata=MN2[which(as.character(MN2$Well)%in%unique_wells[10*(i-1)+1:10*i]),]
  
  bp <- ggplot(subdata, aes(x=Well, y=GEDIratio, fill=Well)) + 
    geom_violin(trim=FALSE) + theme_classic()
  
  jpeg(paste0("boxplot_well_GEDI_MN2_subset",i,".jpeg"))
  print(bp)
  dev.off()
}



bp <- ggplot(MN2, aes(x=Timepoint, y=GEDIratio, fill=Timepoint)) + 
  geom_violin(trim=FALSE)+ scale_fill_brewer(palette="Blues") + theme_classic()

jpeg("boxplot_time_GEDI_MN2.jpeg")
print(bp)
dev.off()



MN2$Well=as.factor(MN2$Well)
bp <- ggplot(MN2, aes(x=Well, y=GEDIratio, fill=Well)) + 
  geom_boxplot()+  theme_classic()

jpeg("boxplot_well_GEDI_MN2.jpeg")
print(bp)
dev.off()
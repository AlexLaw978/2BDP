rm(list = ls())
setwd("")#update


if(F){
  install.packages("./TwoBDP_0.5.0.tar.gz",repos = NULL)
  library(TwoBDP)
}

df=read.csv("./TwoBDP/TestData/df.csv")
md=read.csv("./TwoBDP/TestData/md.csv")

#setup for 2bdp
rownames(df)=df[,1];df=df[,-1]
colnames(df)=paste("V",1:ncol(df),sep="") #to make sure the names are usable
rownames(md)=md[,1]

bdp=createTwoBDPObject(data = df,metadata = md,
                   decisionPointColumn="poi",pointOfInterest="TRUE",threads = 10,saveDirectory="./"
                   )

bdp=panelGeneration(batchSize = 50,bdp = bdp)
bdp=frequencyMap(bdp = bdp)
bdp=PanelExpansion(bdp = bdp)
bdp=panelValidation(bdp = bdp)

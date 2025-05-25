rm(list = ls())
setwd("...")

library(openxlsx)
source("2bdp.R")


df=read.xlsx("exampleData/df.xlsx",rowNames = T)
#rowNames = T becauase sample names are first column and are required to be set as rownames
#if the sample names are not the first column then manaully set the rownames of df to that column and then remove
#ex: rownames(df)=df[,num]; df=df[,-num]. where num is the column index of samples  
md=read.xlsx("exampleData/md.xlsx")

bdpObject=create2BDPClass(data = df,metadata = md,sampleIDS = "Sample",bdp = "poi",poi = "F",
                          rfDataSize = ncol(df))

#default rfDataSize uses 50% of all factors due to the expectation that lots of factors will be used
#however this example dataset is small enough where using all factors is acceptable so the additional parameter is changed from default.


bdpObject=run2BDP(bdpObject)
####Setup####
#rm(list=ls())

####Libraries####
#data manipulation
require(stringr)
require(utils)
require(openxlsx)
#machine learning
require(caret)
require(randomForest)
require(pROC)
require(ResourceSelection)
#parallel
require(parallel)
#graphing
library(ggplot2)
####Functions####

##support##
checkClass=function(v){
  if(is.null(v) || !is(v,"2BDP")){
    message("object not valid\nEnter Valid 2BDP class")  
    message("Can be created through create2BDPClass function")
    return(F)
  }
  return(T)
}

get_mean_diff = function(df,poi){#used in getTopProbs
  idx_pos = which(df$cutoff == poi)
  idx_neg = which(df$cutoff != poi)
  
  mean_pos = apply(df[idx_pos,-1],2,mean)
  mean_neg = apply(df[idx_neg,-1],2,mean)
  diff_mean = abs(mean_pos - mean_neg)
  head(sort(diff_mean,decreasing = TRUE))
  diff_mean = sort(diff_mean,decreasing = TRUE)
  return(diff_mean)
}

getModelInformation=function(fit,data,kfold,nullmod,observedValues,numOfGenes,seed,totalFeatures){
  set.seed(seed)
  p = hoslem.test(as.numeric(data$cutoff),fitted(fit),g=kfold)$p.value
  r2 = with(summary(fit), 1 - deviance/null.deviance)
  r2.adjusted = 1-(1-r2)*(observedValues-1)/(observedValues-numOfGenes-1)
  weights = summary(fit)$coefficients[,1]
  weights=c(weights,rep(NA,totalFeatures-length(weights)+1))
  wnm=c("Intercept",paste("Feature",1:totalFeatures,sep = ""))
  result.model = c(p,r2,r2.adjusted,weights)#,error
  names(result.model) = c("P.Value","McFadden’s.R.Squared","Adjusted.R.Square",wnm)#,"Validation.Error"
  return(result.model)
}

validationMetrics=function(pred,roc,test){
  predicted_values=factor(unname(ifelse(pred>0.5,1,0)),levels = c(0,1))
  actual_values=factor(as.numeric(test$cutoff) - 1, levels =c(0,1))
  conf_matrix=table(predicted_values,actual_values)
  sens = sensitivity(conf_matrix)
  spec = specificity(conf_matrix)
  auc=roc$auc
  return(c(AUC=auc,Specificity=spec,Sensitivity=sens))
}

kfoldCrossValidation=function(seed,fm,train,test,trControl,pointOfInterest,
                              kfold,nullmod,observedValues,numOfGenes,totalFeatures){
  set.seed(seed)
  model = train(fm,train,method = "glm",family = "binomial",trControl = trControl)
  pred=predict.train(model,newdata = test,type = "prob")
  pred=pred[,pointOfInterest]
  resultROC=roc(as.numeric(test$cutoff),pred,quiet = T)
  metrics=validationMetrics(pred,resultROC,test)
  mi=getModelInformation(model,train,kfold,nullmod,observedValues,numOfGenes,seed,totalFeatures)
  metrics=c(metrics,mi)
  return(list(keyMetrics=metrics,ROC=resultROC))
}

rsbmr=function(seed,fm,train,test,
               kfold,nullmod,observedValues,numOfGenes,totalFeatures){
  set.seed(seed)
  model = glm(fm,family = "binomial",data = train)
  pred = predict(model,test,type="response")
  resultROC = roc(as.numeric(test$cutoff),pred,quiet = T)
  metrics=validationMetrics(pred,resultROC,test)
  mi=getModelInformation(model,train,kfold,nullmod,observedValues,numOfGenes,seed,totalFeatures)
  metrics=c(metrics,mi)
  return(list(keyMetrics=metrics,ROC=resultROC))
}

getROC=function(bdp,features){
  features=paste(features,collapse = "+")
  fm=formula(paste("cutoff",features,sep="~"))
  nog=length(str_split(features,"\\+")[[1]])
  #get kfold
  set.seed(bdp$seed)
  tr = createDataPartition(y=bdp$data$cutoff,p=bdp$trainSplitRatio,list = F)
  train = bdp$data[tr,]
  test = bdp$data[-tr,]
  nullmod = glm("cutoff~1", data = train, family="binomial")
  
  kfcvResult=kfoldCrossValidation(seed=bdp$seed,fm=fm,
                                  train=train,test=test,
                                  trControl=bdp$trainControl,
                                  pointOfInterest=bdp$pointOfInterest,
                                  kfold=bdp$kfold,nullmod,observedValues=nrow(train),
                                  numOfGene=nog,totalFeatures=max(bdp$amountOfFeatures))[["ROC"]]
  
  #rsbmr
  rsbmrResult=c() 
  seed=bdp$seed
  count=0
  fdf=data.frame()
  repeat{
    set.seed(seed)
    tr = createDataPartition(y=bdp$data$cutoff,p=bdp$trainSplitRatio,list = F)
    train = bdp$data[tr,]
    test = bdp$data[-tr,]
    nullmod = glm("cutoff~1", data = train, family="binomial")
    count=count+1
    results_rsbmr=rsbmr(seed=seed,fm=fm,train = train,test = test,
                        kfold=bdp$kfold,nullmod,observedValues=nrow(train),
                        numOfGene=nog,totalFeatures=max(bdp$amountOfFeatures))[["ROC"]]
    rsbmrResult=append(rsbmrResult,list(results_rsbmr))
    if(count>=bdp$kfold){
      break
    }
    seed=bdp$seed+count
  }
  
  return(list(rsbmr=rsbmrResult,kfold=kfcvResult))
}

plotROC=function(r,title){
  g=ggroc(r,legacy.axes = T,colour="black")+geom_abline()+ggtitle(title)+
    theme(legend.position="none",panel.background = element_blank(),
          panel.border = element_rect(colour = "black",fill = NA),
          panel.grid.major = element_line(colour = "light grey"),
          panel.grid.minor = element_line(colour = "light grey"))
  return(g)
}

validation=function(i,dex,bdp){
  require(caret);require(base);require(stats);require(stringr)
  require(pROC);require(ResourceSelection)#
  result=list()
  dex=dex+i
  if(dex>length(bdp$biomarkerCombinations)){return(NULL)}
  genes=bdp$biomarkerCombinations[dex]
  fm=formula(paste("cutoff",genes,sep="~"))
  nog=length(str_split(genes,"\\+")[[1]])
  
  
  #validation information
  if(bdp$validationMethod=="all" || validationMethod=="kfcv"){
    set.seed(bdp$seed)
    tr = createDataPartition(y=bdp$data$cutoff,p=bdp$trainSplitRatio,list = F)
    train = bdp$data[tr,]
    test = bdp$data[-tr,]
    nullmod = glm("cutoff~1", data = train, family="binomial")
    
    kfcvResult=kfoldCrossValidation(seed=bdp$seed,fm=fm,
                                    train=train,test=test,
                                    trControl=bdp$trainControl,
                                    pointOfInterest=bdp$pointOfInterest,
                                    kfold=bdp$kfold,nullmod,observedValues=nrow(train),
                                    numOfGene=nog,totalFeatures=max(bdp$amountOfFeatures))[["keyMetrics"]]
  }else{
    kfcvResult=NULL
  }
  
  if(bdp$validationMethod=="all" || validationMethod=="rsbmr"){
    rsbmrInfo=c()
    seed=bdp$seed
    count=0
    fdf=data.frame()
    repeat{
      set.seed(seed)
      tr = createDataPartition(y=bdp$data$cutoff,p=bdp$trainSplitRatio,list = F)
      train = bdp$data[tr,]
      test = bdp$data[-tr,]
      nullmod = glm("cutoff~1", data = train, family="binomial")
      count=count+1
      results_rsbmr=rsbmr(seed=seed,fm=fm,train = train,test = test,
                          kfold=bdp$kfold,nullmod,observedValues=nrow(train),
                          numOfGene=nog,totalFeatures=max(bdp$amountOfFeatures))[["keyMetrics"]]
      fdf=as.data.frame(rbind(fdf,results_rsbmr))
      if(count>=bdp$kfold){
        colnames(fdf)=names(results_rsbmr)
        avg=as.numeric(apply(fdf,2,mean,na.rm=T))
        names(avg)=paste("Average",colnames(fdf),sep = ".")
        stdev=as.numeric(apply(fdf,2,sd,na.rm=T)) 
        names(stdev)=paste("StandardDev",colnames(fdf),sep = ".")
        rsbmrInfo=c()
        for(i in 1:length(avg)){
          rsbmrInfo=c(rsbmrInfo,avg[i],stdev[i])
        }
        break
      }
      seed=bdp$seed+count
    }
    rsbmrResult=rsbmrInfo
  }else{
    rsbmrResult=NULL
  }
  return(list(index=dex,genes=genes,size=nog,kfcv=kfcvResult,rsbmr=rsbmrResult))
}

generatePanel=function(x,bdp,seed){
  require(caret);require(base);require(randomForest)
  sd=seed+x
  message("Generating panel from seed: ",sd," (process: ",x,")")
  set.seed(sd)
  tr = createDataPartition(y=bdp$data$cutoff,p = bdp$trainSplitRatio,list = FALSE)
  diff_mean = get_mean_diff(bdp$data[tr,],bdp$pointOfInterest)
  if(bdp$rfDataSize<=1){
    keep = c("cutoff",names(diff_mean[1:as.integer(ncol(bdp$data)*bdp$rfDataSize)]))
  }else{
    keep = c("cutoff",names(diff_mean[1:as.integer(bdp$rfDataSize)]))
  }
  #keep = c("cutoff",names(diff_mean[1:(ncol(bdp$data)*bdp$rfDataSize-1)]))
  rf = train(cutoff ~., data = bdp$data[,keep],
             method="rf",metric=bdp$metric,
             tuneGrid=bdp$tuneGrid, trControl=bdp$trainControl,importance=TRUE)
  imp = data.frame(importance(rf$finalModel))
  imp = imp[order(imp$MeanDecreaseGini,decreasing = TRUE),]
  topProbes = rownames(imp)[1:max(bdp$amountOfFeatures)]
  return(list(seed=sd,topProbes=topProbes))
}

##main##

#sampleIDS is the column name in metadata to link samples between md and data
#subGrouping=c("all","sw")
#sw = slidingwindow
#for tuneGride and trainControl, you can do either the whole obj or the parms
#amountOfFeatures must be a list of numbers
#validationMethod can only be c("kfcv","rsbmr","all")
#kfcv = kfold cross-validation
create2BDPClass=function(data=NULL,
                         metadata=NULL,bdp=NULL,poi=NULL,sampleIDS=NULL,fileBasename=NULL,
                         featureNameMap=NULL,ifFeatureMap=T,validationMethod="all",
                         seed=100,kfold=10,threads=round(detectCores()*.75),
                         totalPanels=2000,topPanels=200,trainSplitRatio=.7,rfDataSize=.5,
                         amountOfFeatures=2:10,subGrouping=c("all","sw")[1],metric="Accuracy",
                         tuneGrid=NULL,tuneGridParms=2:20,
                         trainControl=NULL,
                         trainControlParms=c(method="repeatedcv",
                                             repeats=5,savePredictions = "all"),
                         AUC=.8,pValue=.05){
  
  #check if data is present
  inputCheck=function(){
    if(is.null(data)){
      message("missing data")
      return(FALSE)
    }
    if(is.null(metadata)){
      message("missing metadata")
      return(FALSE)
    }
    if(any(metadata[,sampleIDS]!=rownames(data))){
      message("samples between Metadata and data are not matching, pls match rows")
      return(FALSE)
    }
    if(is.null(bdp) || !(bdp %in% colnames(metadata))){
      message("missing valid decision point var")
      return(FALSE)
    }
    if(length(unique(metadata[,bdp]))!=2){
      message("Selected binomail decision point does not contain 2 unique options")
      return(FALSE)
    }
    if(kfold<1){
      message("Kfold is to low")
      return(FALSE)
    }
    if(totalPanels<0){
      message("totalPanels can not be negitive")
      return(FALSE)
    }
    if(topPanels<1){
      message("totalPanels can less than 1")
      return(FALSE)
    }
    if(!any(c("all","sw")==subGrouping)){
      message("invalid grouping")
      return(FALSE)
    }
    if(!any(c("Accuracy")==metric)){
      message("invalid metric")
      return(FALSE)
    }
    if(is.null(poi)){
      message("Missing point of interest")
      return(FALSE)
    }
    if(!(poi %in% metadata[,bdp])){
      message("Point of interest is not found in decision point")
      message("Make poi = a var from the bdp vector")
      return(FALSE)
    }
    if(1<rfDataSize && rfDataSize<tail(amountOfFeatures[order(amountOfFeatures)],1)){
      message("rfDataSize is larger than 1 yet smaller than the largest feature list")
      message("Please make rfDataSize a percent or a larger number than amount of features")
      return(FALSE)
    }
    if(rfDataSize<=1 && as.integer(ncol(data)*rfDataSize)<tail(amountOfFeatures[order(amountOfFeatures)],1)){
      message(paste(as.integer(ncol(data)*rfDataSize),"is less than largest amount of features"))
      message("Please increase percentage of rfDataSize")
      return(FALSE)
    }
    if(rfDataSize>ncol(data)){
      message("rfDataSize is larger than max amount of features")
      message(paste("Please set rfDataSize between 0 &",ncol(data)+1))
      return(FALSE)
    }
    return(TRUE)
  }
  if(!inputCheck()){return(NULL)}
  
  if(any(!grepl("^[A-Za-z0-9]+$",colnames(data)))){#check for invalid chars in feature names
    message("invalid symbols found in feature names")
    if(!ifFeatureMap){
      message("Feature Mapping disabled")
      message("enable ifFeatureMap or rename features to numeric and alphabetic letters only")
      return(NULL)
    }
    message("Mapping names to valid name scheme")
    feats=colnames(data)
    ids=paste("id",1:length(feats),sep = "")
    featureNameMap=as.data.frame(cbind(ids=ids,features=feats))
    rm(feats,ids)
    colnames(data)=featureNameMap$ids
  }
  
  if(is.null(tuneGrid)){
    tuneGrid=expand.grid(.mtry=tuneGridParms)
  }
  if(is.null(trainControl)){
    trainControl=trainControl(method=trainControlParms["method"], number=kfold, 
                              repeats=trainControlParms["repeats"],
                              savePredictions=trainControlParms["savePredictions"])
  }
  if(!is.null(fileBasename)){
    fileBasename=paste(fileBasename,".xlsx",sep = "")
    if(!file.exists(fileBasename)){
      message("files doesnt exist, making new workbook")
      wb=createWorkbook()
    }else{
      message("files exists, loading workbook")
      wb=loadWorkbook(fileBasename)
    }
    
    if(!"featureNameMap"%in%wb$sheet_names){
      addWorksheet(wb,"featureNameMap")
    }
    
    if(!is.null(featureNameMap)){
      writeData(wb,"featureNameMap",featureNameMap)
    }else{
      writeData(wb,"featureNameMap",data.frame())
    }
    saveWorkbook(wb,fileBasename,overwrite = T)
  }
  
  data=as.data.frame(cbind(metadata[,bdp],data))
  colnames(data)[1]="cutoff"
  data[,1]=factor(data[,1])
  if(length(levels(data[,1]))!=2){
    message("Require a 2 factor data set");return(NULL)
  }
  bl=levels(data[,1])[levels(data[,1])!=poi]
  data[,1]=relevel(data[,1],bl)
  
  
  #make class
  BDP=list(data=data,metadata=metadata,dbp=bdp,sampleIDS=sampleIDS,pointOfInterest=poi,
           validationMethod=validationMethod,featureNameMap=featureNameMap,seed=seed,
           kfold=kfold,threads=threads,rfDataSize=rfDataSize,totalPanels=totalPanels,
           topPanels=topPanels,trainSplitRatio=trainSplitRatio,fileName=fileBasename,
           amountOfFeatures=amountOfFeatures,subGrouping=subGrouping,metric=metric,
           tuneGrid=tuneGrid,trainControl=trainControl,AUC=AUC,pValue=AUC)
  class(BDP)="2BDP"
  return(BDP)
}


ParallelPanelGeneration=function(bdp=NULL,batchSize=1){
  message("Beginning Panel Generation")
  if(!checkClass(bdp)){
    return(NULL)
  }
  if(!is.null(bdp$panels)){
    message("Panels already created for this 2bdp")
  }
  if(is.null(batchSize) || batchSize<1){
    message("Invalid batchSize - Setting batchSize to 1")
    batchSize==bdp$totalPanels
  }
  
  if(!is.null(bdp$fileName) && file.exists(bdp$fileName)){
    wb=loadWorkbook(bdp$fileName)
    sn="panels"
    if(sn %in% wb$sheet_names){
      message("found file with panels, loading data")
      finalPanels=read.xlsx(wb,sn,rowNames = T)
    }else{
      message("found file with no panel data\ncreating new dataset")
      addWorksheet(wb,sn)
      finalPanels=data.frame()
    }
    if(nrow(finalPanels)>=bdp$totalPanels){
      message(bdp$fileName," has requested amount of panels\nAttaching panels to 2BDP Object")
      bdp$panels=finalPanels
      return(bdp)
    }
  }else{
    finalPanels=data.frame()
  }
  
  message("creating clusters")
  cl=makeCluster(bdp$threads)
  clusterSetRNGStream(cl,bdp$seed)
  clusterExport(cl,c("get_mean_diff"))#,"bdp"
  message("Generating Panels in batches, will take time some time")
  st=Sys.time()
  bc=0
  thrs=(bdp$threads * batchSize)
  while(nrow(finalPanels)<bdp$totalPanels){
    seed=bdp$seed+nrow(finalPanels)
    if(bdp$totalPanels-nrow(finalPanels)<thrs){#prevent getting more panels than requested
      thrs=bdp$totalPanels-nrow(finalPanels)
    }
    bc=bc+1
    bst=Sys.time()
    
    result=parLapply(cl,1:thrs,generatePanel,bdp,seed)
    
    brt=Sys.time()-bst
    
    sds=as.integer(lapply(result,function(v){
      return(v[[1]])
    }))
    panels=as.data.frame(lapply(result,function(v){
      return(v[[2]])
    }))
    panels=as.data.frame(t(panels))
    rownames(panels)=sds
    finalPanels=rbind(finalPanels,panels)
    if(!is.null(bdp$fileName)){#update file with completed batch
      writeData(wb,sn,finalPanels,rowNames = T)
      saveWorkbook(wb,bdp$fileName,overwrite = T)
    }
    message("Finished processing batch ",bc," with:");print(round(brt,2))
    message(nrow(finalPanels),"/",bdp$totalPanels," completed")
  }
  rt=Sys.time()-st
  stopCluster(cl)
  message("Finished processing",bdp$totalPanels,"panels with:")
  print(round(rt,2))
  bdp$panels=finalPanels
  message("Finished Panel Generation completed")
  return(bdp)
}

ParallelFrequencyMapping=function(bdp=NULL){
  message("Beginning Frequency Calculations and Mapping")
  if(!checkClass(bdp)){
    message("Require 2bdp object")  
    return(NULL)
  }
  if(is.null(bdp$panels)){
    message("Missing panels to compute frequencies with - add panels or run *PanelGeneration")
    return(bdp)
  }
  if(!is.null(bdp$frequencyMap)){
    message("Frequency map already exists")
    return(bdp)
  }
  
  if(bdp$threads>ncol(bdp$panels)){#create minimal required amount of threads
    thr=ncol(bdp$panels)
  }else{thr=bdp$threads}
  message("Creating clusters")
  cl=makeCluster(thr)
  clusterSetRNGStream(cl,bdp$seed)
  message("calulating frequencies")
  freqs=parLapply(cl,1:ncol(bdp$panels),function(i,bdp){
    c=bdp$panels[,i]
    uid=unique(c)
    freq=c()
    for(id in uid){
      freq=c(freq,sum(c==id))
    }
    names(freq)=uid
    return(freq)
  },bdp)
  stopCluster(cl)
  
  message("Mapping frequencies")
  freqMap=bdp$panels
  for(i in 1:ncol(freqMap)){
    for(j in 1:length(freqs[[i]])){
      freqMap[freqMap[,i]==names(freqs[[i]])[j],i]=freqs[[i]][j]
    }
  }
  sums=apply(freqMap,1,function(v){sum(as.numeric(v))})
  odr=order(sums,decreasing = T)
  freqMap=freqMap[odr,]
  bdp$panels=bdp$panels[odr,]
  bdp$frequencyMap=freqMap
  message("Finished Frequency Calculations and Mapping")
  return(bdp)
}

ParallelPanelExpansion=function(bdp=NULL){
  message("Beginning Panel Expansion")
  if(!checkClass(bdp)){
    message("Require 2bdp object")  
    return(NULL)
  }
  if(is.null(bdp$panels)){
    message("Missing panels to compute feature combinations - add panels or run *PanelGeneration")
    return(bdp)
  }
  if(!bdp$subGrouping %in% c("all","sw")){
    message("invalid subGrouping of panels")
    return(bdp)
  }
  if(!is.null(bdp$biomarkerCombinations)){
    message("Biomarker Combinations already created")
  }
  if(bdp$topPanels>nrow(bdp$panels) || bdp$topPanels<1){
    message("top panels exceed total panels, must be less than all panels and greater than 0")
    return(bdp)
  }
  
  message("Creating Clusters")
  cl=makeCluster(bdp$threads)
  clusterSetRNGStream(cl,bdp$seed)
  message("Generating feature combinations")
  featureCombinations=parLapply(cl,1:bdp$topPanels,function(i,bdp){
    require(base);require(utils)
    rw=bdp$panels[i,]
    comboList=c()
    for(j in bdp$amountOfFeatures){
      if(j>length(rw)){break}
      if(bdp$subGrouping=="all"){
        cdf=combn(rw,j)
        grps=as.character(apply(cdf,2,function(cl){
          cl=as.character(cl)
          cl=cl[order(cl)]
          return(paste(cl,collapse = "+"))
        }))
        comboList=c(comboList,grps)
      }else if(bdp$subGrouping=="sw"){
        grp=rw[1:j]
        grp=grp[order(grp)]
        grp=paste(grp,collapse = "+")
        comboList=c(comboList,grp)
      }
    }
    return(comboList)
  },bdp)
  stopCluster(cl)
  
  message("cleaning combination results")
  finalCombos=c()
  for(i in featureCombinations){
    finalCombos=c(finalCombos,i)
    finalCombos=unique(finalCombos)#this is inside loop to conserve as much memory as possible
  }
  
  bdp$biomarkerCombinations=finalCombos
  
  message("Finished Panel Expansion")
  return(bdp)
}

ParallelBiomarkerValidation=function(bdp=NULL,batchSize=1,saveROCs=T){
  message("Beginning Validation")
  #data checking
  if(!checkClass(bdp)){
    message("Require 2bdp object")  
    return(NULL)
  }
  if(is.null(bdp$biomarkerCombinations)){
    message("Missing combinations in bdp object")
    message("Use *PanelExpansion to generate combinations")
    return(bdp)
  }
  if(is.null(batchSize) || batchSize<1){
    message("Invalid batchSize - Setting batchSize to 1")
    batchSize==bdp$totalPanels
  }
  if(!is.null(bdp$fileName)){
    message("Found file, loading data")
    wb=loadWorkbook(bdp$fileName)
    rI="rsbmrInfo"
    kI="kfcvInfo"
    getBDPTables=function(workBook,sheetName){
      if(sheetName %in% workBook$sheet_names){
        df=read.xlsx(workBook,sheetName)
      }else{
        addWorksheet(workBook,sheetName)
        df=data.frame()
      }
      return(list(workBook,df))
    }
    tmp=getBDPTables(wb,rI)
    wb=tmp[[1]]
    rsbmrTable=tmp[[2]]
    tmp=getBDPTables(wb,kI)
    wb=tmp[[1]]
    kfcvTable=tmp[[2]]
  }else{
    rsbmrTable=data.frame()
    kfcvTable=data.frame()
  }
  
  
  
  message("Creating clusters")
  cl=makeCluster(bdp$threads)
  clusterSetRNGStream(cl,bdp$seed)
  clusterExport(cl,c("getModelInformation","kfoldCrossValidation","rsbmr","validationMetrics"))
  st=Sys.time()
  bc=0
  thrs=(bdp$threads * batchSize)
  #batch processor
  message("Validating Panels in batches, will take time some time to process")
  repeat{
    if(switch(
      bdp$validationMethod,
      "all"={nrow(rsbmrTable)>=length(bdp$biomarkerCombinations) && nrow(kfcvTable)>=length(bdp$biomarkerCombinations)},
      "rsbmr"={nrow(rsbmrTable)>=length(bdp$biomarkerCombinations)},
      "kfcv"={nrow(kfcvTable)>=length(bdp$biomarkerCombinations)}
    )){
      break
    }
    bc=bc+1
    bst=Sys.time() 
    panelValidations=parLapply(cl,1:thrs,validation,
    min(nrow(rsbmrTable),nrow(kfcvTable)), bdp)
    # #TESTING single Thread
    # panelValidations=lapply(1:thrs,validation,nrow(modelTable),bdp)
    brt=Sys.time()-bst
    message("Finished processing batch ",bc," with:");print(round(brt,2))
    message("Collaping results")
    for(i in panelValidations){
      if(is.null(i)){next}
      #pull data
      rowNum=i[["index"]]
      genes=i[["genes"]]
      numOfGenes=i[["size"]]
      kfcvRes=i[["kfcv"]]
      rsbmrRes=i[["rsbmr"]]
      
      sInfo=c(Index=rowNum,Features=genes,NumberOfGenes=numOfGenes)
      
      rsbmrRes=c(sInfo,rsbmrRes)
      rsbmrTable=as.data.frame(rbind(rsbmrTable,rsbmrRes))
      colnames(rsbmrTable)=names(rsbmrRes)
      kfcvRes=c(sInfo,kfcvRes)
      kfcvTable=as.data.frame(rbind(kfcvTable,kfcvRes))
      colnames(kfcvTable)=names(kfcvRes)
    }
    
    if(!is.null(bdp$fileName)){
      message("Saving progress to ",bdp$fileName)
      writeData(wb,rI,rsbmrTable)
      writeData(wb,kI,kfcvTable)
      saveWorkbook(wb,bdp$fileName,overwrite = T)
    }
    message(min(nrow(rsbmrTable),nrow(kfcvTable)),
            "/",length(bdp$biomarkerCombinations)," left")
    
  }
  rt=Sys.time()-st
  message("Finished processing  ",bc," with:");print(round(rt,2))
  stopCluster(cl)
  #rsbmrTable=rsbmrTable[order(rsbmrTable$Average.AUC,decreasing = T),]
  bdp$rsbmrValidation=rsbmrTable
  #kfcvTable=kfcvTable[order(kfcvTable$AUC,decreasing = T),]
  bdp$kfcvValidation=kfcvTable
  
  message("Finished Validating Feature List")
  return(bdp)
}

trimByThreshold=function(bdp=NULL,pth=.05,ath=.8){
  message("Beginning trim & merge of tables for best results by cutoff values")
  if(!checkClass(bdp)){
    return(NULL)
  }
  
  ri=bdp$rsbmrValidation
  ki=bdp$kfcvValidation
  ri=ri[ri$Average.P.Value <= pth & ri$Average.AUC >= ath,]
  ki=ki[ki$P.Value <= pth & ki$AUC >= ath,]
  
  ki=ki[order(ki$AUC,decreasing = T),]
  ri=ri[order(ri$Average.AUC,decreasing = T),]
  
  if(!is.null(bdp$fileName)){
    message("Saving data to ",bdp$fileName)
    wb=loadWorkbook(bdp$fileName)
    sn1="bestKFOLD"
    sn2="bestRSMBR"
    if(!sn1 %in% wb$sheet_names){addWorksheet(wb,sn1)}
    if(!sn2 %in% wb$sheet_names){addWorksheet(wb,sn2)}
    writeData(wb,sn1,ki);writeData(wb,sn2,ri)
    saveWorkbook(wb,bdp$fileName,overwrite = T)
  }
  message("finished trimming")
  bdp$bestKFOLD=ki;bdp$bestRSMBR=ri
  return(bdp)
}

convertNames=function(bdp=NULL,updateTables=T){
  message("converting id names to feature names")
  if(!checkClass(bdp)){
    return(NULL)
  }
  if(is.null(bdp$featureNameMap)){
    #message("names already converted")
    return(bdp)
  }
  #the three tables should be equal to each other
  if(!is.null(bdp$rsbmrValidation)){
    fts=bdp$rsbmrValidation$Features
  }else if(!is.null(bdp$kfcvValidation)){
    fts=bdp$kfcvValidation$Features
  }
  for(i in 1:9){
    con = any(grepl(bdp$featureNameMap$ids[i],fts))
    if(con){break}
  }
  if(!con){
    message("no ids found, assuming conversion has been done already")
    return(bdp)
  }
  
  fts=str_split(fts,"\\+")
  #making assumption that all feature ids exist in this feature map
  
  fts=as.character(lapply(fts,function(ft,mp=bdp$featureNameMap){
    for(i in 1:length(ft)){
      ft[i]=mp$features[mp$ids==ft[i]]
    }
    rs=paste(ft,collapse = ", ")
    return(rs)
  }))
  tbs=c("rsbmrValidation","kfcvValidation")
  for(i in tbs){
    bdp[[i]]$Features=fts
  }
  
  
  if(updateTables & !is.null(bdp$fileName)){
    message("Saving data to ",bdp$fileName)
    fns=c("rsbmrInfo","kfcvInfo")
    wb=loadWorkbook(bdp$fileName)
    for(i in 1:length(fns)){
      writeData(wb,fns[i],bdp[[tbs[i]]])
    }
    saveWorkbook(wb,bdp$fileName,overwrite = T)
  }
  
  message("finished name conversion")
  return(bdp)
}

run2BDP=function(bdp,batchSize=1){
  if(!checkClass(bdp)){
    message("invalid 2bdp object")
    message("check that object creation was correctly done")
    return(NULL)
  }
  bdp=ParallelPanelGeneration(bdp=bdp,batchSize = batchSize)
  gc()
  bdp=ParallelFrequencyMapping(bdp = bdp)
  bdp=ParallelPanelExpansion(bdp = bdp)
  bdp=ParallelBiomarkerValidation(bdp = bdp,batchSize = batchSize)
  gc()
  bdp=convertNames(bdp=bdp)
  bdp=trimByThreshold(bdp=bdp)
  gc()
  return(bdp)
}
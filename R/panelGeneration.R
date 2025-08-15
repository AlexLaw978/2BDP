#' Generate a panel for TwoBDP
#'
#' @param i An integer to offset seed value in modelParameters.
#' @param modelParameters A list of parameters that are used to generate a panel.
#' @param maxFeatures The max size of features required. Used as cutoff to select features from random forest.
#' @param getModel A logical to know if the model should be returned
#'
#' @return TwoBDP class with processed data appended.
#'
#' @description
#'  Primary function that is used to generate a panel using random forest
#'
#' @import randomForest
#' @importFrom caret train
#'
#' @rdname generatePanel
#' @export
generatePanel=function(i,modelParameters,maxFeatures,getModel=F){
  #library(randomForest)
  #library(caret,include.only = c("train","createDataPartition"))

  rngkind=modelParameters$seedParms$RNGkind
  seed_base=modelParameters$seedParms$seed
  seed_data=seed_base
  seed_model=seed_base

  if(modelParameters$generationType["data"]){
    seed_data=seed_data+i
  }
  if(modelParameters$generationType["model"]){
    seed_model=seed_model+i
  }

  #we arent here to validate the random forest, so there isnt any need to subset the data.
  #will include option though

  if(!is.null(modelParameters$indexes) && length(modelParameters$indexes)==1){
    data=modelParameters$data[as.numeric(modelParameters$indexes[["0"]][-c(1,2)]),]
    seed_data=seed_base
  }else if(!is.null(modelParameters$indexes) && length(modelParameters$indexes)!=1){
    data=modelParameters$data[as.numeric(modelParameters$indexes[[as.character(i)]][-c(1,2)]),]
  }else{
    seed_data=NA #since no subset was performed, seed_data was not used or needed to regenerate data
    data=modelParameters$data
  }

  data=cbind.data.frame(POI=modelParameters$pointsOfInterest,data)
  #data$POI=factor(data$POI)

  setRNG(seed = seed_model,rngKind = rngkind)
  rf = train(POI ~., data = data,
             method="rf",metric=modelParameters$metric,
             tuneGrid=modelParameters$tuneGrid, trControl=modelParameters$trainControl,importance=TRUE)
  imp = data.frame(importance(rf$finalModel))
  imp = imp[order(imp$MeanDecreaseGini,decreasing = TRUE),]
  TopFeatures = rownames(imp)[1:maxFeatures]
  result=c(seed_data=seed_data,seed_model=seed_model,TopFeatures=TopFeatures)
  if(getModel){
    result=list()
    result$randomForest=rf
    result$allProbes=imp
    result$seed_data=seed_data
    result$seed_model=seed_model
    result$TopFeatures=TopFeatures
  }
  return(result)
}

#' Two BDP Automation
#'
#' @param bdp The TwoBDP object to process
#' @param batchSize How many datasets to process before saving. Only meaningful when progress is saved.
#'
#' @return TwoBDP class with processed data appended.
#'
#' @description
#' These functions are used to automate, parallize, and process the data in batches.
#' Once the function has been completed, the TwoBDP object will then have the required data
#' needed for the next step.
#' Order of functions are as follows:
#' bdp=createTwoBDPObject(bdp)
#' bdp=panelGeneration(bdp)
#' bdp=frequencyMap(bdp);bdp=PanelExpansion(bdp)
#' bdp=panelValidation(bdp)
#'
#'
#' @examples
#' df=read.csv("../TwoBDP/TestData/df.csv")
#' md=read.csv("../TwoBDP/TestData/md.csv")
#' rownames(df)=df[,1];df=df[,-1]
#' colnames(df)=paste("V",1:ncol(df),sep="")
#' rownames(md)=md[,1]
#'
#' bdp=createTwoBDPObject(data = df,metadata = md,
#'   decisionPointColumn="poi",pointOfInterest="TRUE",threads = 10,saveDirectory = ".",
#'   panelGenerationParameters=TwoBDP_makePanelParameterList(totalPanels=100)
#'
#' bdp=panelGeneration(batchSize = 50,bdp = bdp)
#' bdp=frequencyMap(bdp = bdp)
#' bdp=PanelExpansion(bdp = bdp)
#' bdp=panelValidation(bdp = bdp)
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom stringr str_count
#'
#' @rdname TwoBDPAutomation
#' @export
panelGeneration=function(batchSize=5,bdp){

  if(is.null(bdp$panelGenerationParameters)){
    return(message_missingValues_function("bdp$panelGenerationParameters","panelGeneration"))
  }

  panelCount=0
  totalPanels=bdp$panelGenerationParameters$totalPanels
  batchSize=batchSize*bdp$threads
  batchesFinished=0
  fp="RF"
  modelParameters=generateModelParametersFromTwoBDPClass(bdp,"rf")
  bdp$panelGenerationParameters$modelParameters=modelParameters #saving it in case user needs it later
  finalFile=paste(bdp$saveDirectory,"GeneratedPanels.csv",sep = "/")
  panels=data.frame()

  #--------------------------------------------------------------------------------------------------------------------------handle directory load
  if(file.exists(finalFile)){
    message(paste(finalFile,"exists."))
    message("Loading file and appending to TwoBDP object")
    bdp$panels=read.csv(finalFile,row.names = 1)
    message("Will skip generation")
    return(bdp)
  }
  #library(stringr)
  if(!is.null(bdp$saveDirectory)){
    fls=list.files(bdp$workDir,pattern = fp)
    batchesFinished=length(fls)
    panelCount=batchesFinished*batchSize
  }

  #--------------------------------------------------------------------------------------------------------------------------iterator for panel generation
  if(bdp$threads>1){
    cl=makeCluster(bdp$threads)
    clusterExport(cl,c("setRNG"))
    clusterEvalQ(cl,expr = {library(caret,include.only = c("train","createDataPartition"));library(randomForest)})
  }
  message("Starting panel generation")
  while(panelCount<totalPanels){
    startBatch=panelCount
    endBatch=startBatch+batchSize-1 #need -1 here so its exclusive not inclusive
    if(endBatch>totalPanels){endBatch=totalPanels-1}#starting dex is 0 so need -1

    if(bdp$threads>1){
      p_res=parLapply(cl,startBatch:endBatch,generatePanel,
                      modelParameters,max(bdp$frequencyMapParameters$amountOfFeatures))
    }else{
      p_res=lapply(startBatch:endBatch,generatePanel,
                   modelParameters,max(bdp$frequencyMapParameters$amountOfFeatures))
    }

    p_res=convertFlipList(p_res)
    rownames(p_res)=startBatch:endBatch
    panelCount=endBatch+1
    msg_c(c("finished",panelCount,"of",bdp$panelGenerationParameters$totalPanels,"panels"))
    batchesFinished=batchesFinished+1

    #--------------------------------------------------------------------------------------------------------------------------Save panels
    if(is.null(bdp$saveDirectory)){
      panels=rbind.data.frame(panels,p_res)
    }else{#save in batches below
      write.csv(p_res,file = paste(bdp$workDir,
                             createTempFileName(pre=fp,num=batchesFinished),sep = "/")
                )

    }
  }
  if(bdp$threads>1){stopCluster(cl)}

  if(!is.null(bdp$saveDirectory)){
    fls=list.files(bdp$workDir,pattern = fp,full.names = T)
    for(f in fls){
      subpanels=read.csv(f,row.names = 1)
      panels=rbind.data.frame(panels,subpanels)
    }
    write.csv(panels,file = finalFile)

    #remove all temp files as final file has been saved
    for(f in fls){
      file.remove(f)
    }
  }


  message("Finished panel generation")
  bdp$panels=panels
  return(bdp)
}

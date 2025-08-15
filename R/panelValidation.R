#' Generate validation model information
#'
#' @description
#' The key function in validating a panel. It was designed so it can be ran in parallel,
#' so i acts as an seed offset when the function is called in parLapply.
#'
#' @param i An integer to offset the seed by.
#' @param modelParameters A named list that has all arguments required to validate a panel.
#' @param method A string that identify what validation method to use.
#' @param getModel A logical to determine if you want a list of models for external use.
#'
#' @returns A character vector of the model information if getModel is FALSE,
#' otherwise returns a list containing the model information in addition to the actual models.
#'
#' @importFrom stringr str_count
#' @importFrom stringr str_split
#'
#'
#' @rdname generateModelAndInfo
#' @export
generateModelAndInfo=function(i,modelParameters,method,getModel=F){

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

  data=cbind.data.frame(POI=modelParameters$pointsOfInterest,modelParameters$data)

  if(is.na(seed_data)){#this prevents any validation
    return(NULL)
  }
  pdt=modelParameters$featureCombinations[i+1]
  if(method=="kfold"){#------------------------------------------------------------------------kfold
    setRNG(seed = seed_data,rngKind = rngkind)
    indexes=getIndex(data$POI,ratio = modelParameters$ratio,
             seed = seed_data,rngKind = rngkind,
             count = 1)[[1]]
    indexes=as.numeric(indexes[-c(1,2)])
    train=data[indexes,]
    test=data[-indexes,]

    modelData=getModelAndPred(seed = seed_model,rngkind = rngkind,
                              fm=as.formula(paste("POI",pdt,sep = "~")),
                              data = train,test=test,
                              trainControl = modelParameters$trainControl)

    validationInfo=getValidationInformation(modelData$roc)
    modelInfo=getModelInformation(fit=modelData$model,
                                  observedValues=nrow(train),
                                  numOfFeatures=str_count(pdt,"\\+")+1 )

    if(getModel){
      return(list(modelData=modelData,validationInfo=validationInfo,
                  modelInfo=modelInfo,panel=pdt))
    }

    info=c(seed_data=seed_data,seed_model=seed_model,RNGkind=rngkind,
           panel=pdt,validationInfo,modelInfo)

  }else{#-----------------------------------------------------------------------rsbmr

    indexes=getIndex(data$POI,ratio = modelParameters$ratio,
                     seed = seed_data,rngKind = rngkind,
                     count = modelParameters$repeats)
    infos=list();count=1;
    for(dex in indexes){
      dex=as.numeric(dex[-c(1,2)])
      train=data[dex,]
      test=data[-dex,]
      modelData=getModelAndPred(seed = seed_model,rngkind = rngkind,
                                fm=as.formula(paste("POI",pdt,sep = "~")),
                                data = train,test=test,
                                trainControl = modelParameters$trainControl)

      validationInfo=getValidationInformation(modelData$roc)
      modelInfo=getModelInformation(fit=modelData$model,
                                    observedValues=nrow(train),
                                    numOfFeatures=str_count(pdt,"\\+")+1 )
      info=c(subseed=names(indexes)[count],validationInfo,modelInfo)
      infos[[count]]=info;count=count+1
    }
    #collapse data
    #eed_model=seed_model,RNGkind=rngkind,panel=pdt
    rseed=paste(as.character(seed_data),"-",sep = "")
    valInfo=list();statsInfo=list();count=1
    for(j in infos){
      rseed=paste(rseed,as.character(j["subseed"]),sep = "|") #take care of seed
      valInfo[[count]]=j[2:7]#table stats
      statsInfo[[count]]=j[8:length(j)]
      count=count+1
    }
    valInfo=convertFlipList(valInfo)
    statsInfo=convertFlipList(statsInfo)

    valInfo=apply(valInfo,2,function(cl){mean(as.numeric(cl))})
    statsInfo=apply(statsInfo,2,function(cl){
      cldf=convertFlipList(str_split(cl,"\\|\\|\\|"))
      mcl=c(cldf[1,1],mean(as.numeric(cldf[,2])),
            mean(as.numeric(cldf[,3])),mean(as.numeric(cldf[,4])))
      mcl=paste(mcl,collapse = "|||")
      return(mcl)
    })


    info=c(seed_data=rseed,seed_model=seed_model,RNGkind=rngkind,
           panel=pdt,valInfo,statsInfo)
  }


  return(info)
}

#' panelValidation
#'
#' @importFrom stringr str_count
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterExport
#' @rdname TwoBDPAutomation
#' @export
panelValidation=function(batchSize=10,bdp){

  if(is.null(bdp$featureCombinations)){
    message("Missing combinations in bdp object")
    message("Use *PanelExpansion to generate combinations")
    return(bdp)
  }
  if(is.null(batchSize) || batchSize<1){
    message("Invalid batchSize - Setting batchSize to 1")
    batchSize==bdp$totalPanels
  }

  #-----------------------------------------------------------------------------kfold checks
  if("kfold"%in%bdp$validationParameters$validationMethods){
    if(is.null(bdp$validationParameters$trainControls$kfold)){
      message("Missing train control for kfold")
      return(bdp)
    }
  }
  #-----------------------------------------------------------------------------rsbmr checks
  if("rsbmr"%in%bdp$validationParameters$validationMethods){
    if(is.null(bdp$validationParameters$trainControls$kfold)){
      message("Missing train control for rsbmr, returning NULL")
      return(bdp)
    }
  }
  #-----------------------------------------------------------------------------process data

  res_baseCol=10
  res_featSize=max(str_count(bdp$featureCombinations,"\\+"))+2 #' 1 offset + 1 for intercept
  res_colSize=res_featSize+res_baseCol
  totalmodels=length(bdp$featureCombinations)
  batchSize=batchSize*bdp$threads
  maxCount=length(bdp$featureCombinations)
  bdp$validatedFeatures=list()

  for(i in bdp$validationParameters$validationMethods){
    modelCount=0
    batchesFinished=0
    modelParameters=generateModelParametersFromTwoBDPClass(bdp=bdp,method=i)
    finalFile=paste(bdp$saveDirectory,"/",i,"_ValidatedFeatures.csv",sep = "")
    if(i=="kfold"){
      fp="KF"
    }else{fp="RS"}

    if(file.exists(finalFile)){
      message(paste(finalFile,"exists."))
      message("Loading file and appending to TwoBDP object")
      bdp$validatedFeatures[[i]]=read.csv(finalFile,row.names = 1)
      message("Will skip generation")
      next
    }

    if(!is.null(bdp$saveDirectory)){
      fls=list.files(bdp$workDir,pattern = fp)
      batchesFinished=length(fls)
      modelCount=batchesFinished*batchSize
    }

    message(paste("Begining",i,"validation"))

    if(bdp$threads>1){
      cl=makeCluster(bdp$threads)
      clusterExport(cl,c("setRNG","getModelAndPred",
                         "getValidationInformation",
                         "getModelInformation",
                         "getIndex","convertFlipList"))
      clusterEvalQ(cl,expr = {
        library(caret,include.only = c("train","createDataPartition","predict.train"))
        library(randomForest);library(pROC);library(stringr)
        })
    }
    while(modelCount<totalmodels){
      startBatch=modelCount
      endBatch=startBatch+batchSize-1
      if(endBatch>=totalmodels){endBatch=totalmodels-1}

      if(bdp$threads>1){#parallel
        modelInfos=parLapply(cl,startBatch:endBatch,generateModelAndInfo,modelParameters,i)
      }else{#serial
        modelInfos=lapply(startBatch:endBatch,generateModelAndInfo,modelParameters,i)
      }


      #due to the differing sizes from the input formula
      #Ill need to manually join and add buffers so the tables are actually uniformed
      #max size should be 9 + maxfeatures
      p_res=data.frame()
      szs=as.numeric(lapply(modelInfos,length))
      for(j in unique(szs)){
        sdf=convertFlipList(modelInfos[j==szs])
        if(res_colSize-ncol(sdf)!=0){
          ndf=as.data.frame(matrix(NA,nrow = nrow(sdf),ncol = res_colSize-ncol(sdf)))
          colnames(ndf)=paste("FeatureStats",(ncol(sdf)-res_baseCol+1):res_featSize,sep = "")
          sdf=cbind.data.frame(sdf,ndf)
        }
        p_res=rbind.data.frame(p_res,sdf)
      }

      rownames(p_res)=startBatch:endBatch
      modelCount=endBatch+1
      msg_c(c("finished",i,modelCount,"of",maxCount,"feature combinations"))
      batchesFinished=batchesFinished+1

      if(!is.null(bdp$saveDirectory)){
        write.csv(p_res,file = paste(bdp$workDir,
                                     createTempFileName(pre=fp,num=batchesFinished),sep = "/")
        )
      }else{
        bdp$validatedFeatures[[i]]=rbind.data.frame(bdp$validatedFeatures[[i]],p_res)
      }

    }
    if(bdp$threads>1){stopCluster(cl)}

    if(!is.null(bdp$saveDirectory)){
      fls=list.files(bdp$workDir,pattern = fp,full.names = T)
      validatedFeatures=data.frame()
      for(f in fls){
        subfeat=read.csv(f,row.names = 1)
        validatedFeatures=rbind.data.frame(validatedFeatures,subfeat)
      }
      write.csv(validatedFeatures,file = finalFile)
      bdp$validatedFeatures[[i]]=validatedFeatures
      message("Saved final file")
      #remove all temp files as final file has been saved
      for(f in fls){
        file.remove(f)
      }
      message("Cleaned tempory files")
    }

  }

  return(bdp)
}






#' An internal constant thats repeated
#' user never need to call this
#' @keywords internal
constant_panelSubgroupMethods=c("all","slidingWindow")

#' An internal constant thats repeated
#' user never need to call this
#' @keywords internal
default_generationTypes=c(data=T,model=T)

#' An internal constant thats repeated
#' user never need to call this
#' @keywords internal
default_seedParams=list(seed=100,RNGkind="L'Ecuyer-CMRG")

#' createTwoBDPObject
#'
#' @param data Data frame of features "columns" and samples "rows" to analyze. Rownames should be sample identifiers which are all in the rownames in metadata.
#' @param metadata Data frame which has a 2 factor column which will be used as the point of interest.
#' @param saveDirectory A string path where results will be saved.
#' @param decisionPointColumn A string for the column name in metadata, which is the 2 factor for analysis.
#' @param pointOfInterest A string for the non-control value in decisionPointColumn.
#' @param seedParms A named list that contains the "seed" and "RNGkind". see ?RNGkind for more details and options.
#' @param threads An int for the amount of cores use in parallel.
#' @param panelGenerationParameters A list of required parameters to run the random forest. See TwoBDP_makePanelParameterList function for more details.
#' @param frequencyMapParameters A list of required parameters to run the frequency map. See TwoBDP_makeFrequencyParameterList function for more details.
#' @param validationParameters A list of required parameters to run the validation. See TwoBDP_makeValidationParameterList function for more details.
#'
#' @return TwoBDP class that will be uses for all TwoBDP related functions.
#'
#' @examples
#' bdp=createTwoBDPObject(data=df,metadata=metaDF,
#' decisionPointColumn="metadataColumn",pointOfInterest="Case")
#'
#' @rdname TwoBDPClass
#' @export
createTwoBDPObject=function(data=NULL,metadata=NULL,saveDirectory=NULL,
                            decisionPointColumn=NULL,pointOfInterest=NULL,
                            seedParms=default_seedParams,
                            threads=1,
                            panelGenerationParameters=TwoBDP_makePanelParameterList(),
                            frequencyMapParameters=TwoBDP_makeFrequencyParameterList(),
                            validationParameters=TwoBDP_makeValidationParameterList()){

  message("checking 2bdp creation")
  #--------------------------------------------------------------------------------------------------------------------------nulls&types #TODO types
  #check primary nulls
  if(is.null(data)){return(message_missingValues_function("data","TwoBDP"))}
  if(is.null(metadata)){return(message_missingValues_function("metadata","TwoBDP"))}
  if(is.null(decisionPointColumn)){return(message_missingValues_function("decisionPointColumn","TwoBDP"))}
  if(is.null(pointOfInterest)){return(message_missingValues_function("pointOfInterest","TwoBDP"))}

  ##--------------------------------------------------------------------------------------------------------------------------data/metadata
  #check rownames between data & metadata match/exists
  if(nrow(data)!=nrow(metadata)){return(message_notEqualWithin("data","metadata"))}
  else if(!all(rownames(data)%in%rownames(metadata))){return(message_missMatch_dataframe("rownames","data","metadata"))}
  else{#as everything is in place, we will align the rows
    #this will make appending data for models easier
    data=data[rownames(metadata),]
  }

  #--------------------------------------------------------------------------------------------------------------------------decision/point #TODO
  #make sure decisionPointColumn & pointOfInterest is correct
  if(length(decisionPointColumn)!=1 || length(pointOfInterest)!=1){
    message_incorrectSize("decisionPointColumn",length(decisionPointColumn),1,rnull=F)
    return(message_incorrectSize("pointOfInterest",length(pointOfInterest),1))
  }
  if(!decisionPointColumn%in%colnames(metadata)){
    return(message_notEqualWithin_NULL("DecisionPointColumn","metadata",type=2))
  }else if(!pointOfInterest%in%metadata[[decisionPointColumn]]){
    return(message_notEqualWithin("pointOfInterest","DecisionPointColumn",type=2))
  }


  if(!is.factor(metadata[[decisionPointColumn]])){
    dpc=metadata[[decisionPointColumn]]
    if(is.logical(dpc)){
      message("Decision Point vector was found to be a logical, converting it to a usable factor")
      dpc[dpc==FALSE]="control"
      dpc[dpc==TRUE]="target"
    }else if(is.numeric(dpc) || is.integer(dpc)){
      message("Decision Point vector was found to be a numeric, converting it to a usable factor")
      message("Lowest value will be set as control")
      dpc[dpc==min(dpc)]="control"
      dpc[dpc==max(dpc)]="target"
    }else{
      message("Unable to convert DecisionPointColumn in metadata into a valid factor.")
      return(NULL)
    }
    dpc=factor(dpc)
    metadata[[decisionPointColumn]]=dpc
  }

  #--------------------------------------------------------------------------------------------------------------------------panel Params #TODO
  if(!is.null(panelGenerationParameters)){
    if(panelGenerationParameters$rfInputDataSize!=1 && panelGenerationParameters$rfInputDataSize!=ncol(df)){#trim table to top probes
      msg_c(
        c("rfInputDataSize was not found to be the full data set.\nwill trim based on decreasing diff mean calculated from pointOfInterest")
      )
      #data is aligned by this point so we can append pointOfInterest to df to trim data
      namedMeanDiffVector=get_mean_diff(df,metadata[[decisionPointColumn]]==pointOfInterest)
      if(panelGenerationParameters$rfInputDataSize<1){#subset vector by percent or an actual limit
        namedMeanDiffVector=namedMeanDiffVector[
                                1:floor(panelGenerationParameters$rfInputDataSize*length(namedMeanDiffVector))
                                ]
      }else{#actual
        namedMeanDiffVector=namedMeanDiffVector[1:panelGenerationParameters$rfInputDataSize]
      }
      data=data[,names(namedMeanDiffVector)]
    }
  }

  #--------------------------------------------------------------------------------------------------------------------------frequency Parms #TODO
  if(!is.null(frequencyMapParameters)){
    if(!is.null(frequencyMapParameters$append)){
      if(!(all(frequencyMapParameters$append%in%colnames(data)))){
        goodDex=frequencyMapParameters$append%in%colnames(data)
        message("Warning:")
        message(paste("Removing",frequencyMapParameters$append[!goodDex],"from append as they were not found in columns of data"))
        frequencyMapParameters$append=frequencyMapParameters$append[goodDex]
      }
    }
    if(!is.null(frequencyMapParameters$sub)){
      if(!(all(frequencyMapParameters$sub%in%colnames(data)))){
        goodDex=frequencyMapParameters$sub%in%colnames(data)
        message("Warning:")
        message(paste("Removing",frequencyMapParameters$sub[!goodDex],"from sub as they were not found in columns of data"))
        frequencyMapParameters$sub=frequencyMapParameters$sub[goodDex]
      }
      if(length(frequencyMapParameters$sub)>max(frequencyMapParameters$amountOfFeatures)){
        message("sub length exceeds length of max combinations")
        message("sub is being set to NULL")
        frequencyMapParameters$sub=NULL
      }
    }
  }

  #--------------------------------------------------------------------------------------------------------------------------validation Parms #TODO
  if(!is.null(validationParameters)){

  }

  #--------------------------------------------------------------------------------------------------------------------------build list
  bdp=list(data=data,metadata=metadata,saveDirectory=saveDirectory,
           decisionPointColumn=decisionPointColumn,pointOfInterest=pointOfInterest,
           seedParms=seedParms,threads=threads,
           panelGenerationParameters=panelGenerationParameters,
           frequencyMapParameters=frequencyMapParameters,
           validationParameters=validationParameters)

  #--------------------------------------------------------------------------------------------------------------------------directory
  if(!is.null(saveDirectory)){
    workDir=paste(saveDirectory,"temp",sep = "/")
    if(!dir.exists(workDir)){dir.create(workDir, recursive = TRUE)}
    bdp$workDir=workDir
  }
  #--------------------------------------------------------------------------------------------------------------------------finished
  class(bdp)="TwoBDP"
  msg_c("TwoBDP object successfully created")
  return(bdp)
}

#' TwoBDP_makePanelParameterList
#'
#' @param totalPanels Data frame of features "columns" and samples "rows" to analyze. Rownames should be sample identifiers which are all in the rownames in metadata.
#' @param topPanels Data frame which has a 2 factor column which will be used as the point of interest.
#' @param generationType A string path where results will be saved.
#' @param rfInputDataSize A numeric value to determine if the data should be subset before running 2bdp.
#' @param maxPanelSize A string for the column name in metadata, which is the 2 factor for analysis.
#' @param trainSplitRatio A string for the non-control value in decisionPointColumn.
#' @param ntree A named list that contains the "seed" and "RNGkind".
#' @param nodeSize An int for the amount of cores use in parallel.
#' @param metric A list of required parameters to run the random forest. See TwoBDP_makePanelParameterList function for more details.
#' @param tuneGrid A list of required parameters to run the frequency map. See TwoBDP_makeFrequencyParameterList function for more details.
#' @param trainControl A list of required parameters to run the validation. See TwoBDP_makeValidationParameterList function for more details.
#'
#' @return List of named parameters that will be called in panelGeneration function.
#'
#' @importFrom caret trainControl
#'
#' @rdname TwoBDP_makePanelParameterList
#' @export
TwoBDP_makePanelParameterList=function(totalPanels=2000,topPanels=200,
                                       generationType=default_generationTypes,
                                       rfInputDataSize=1,maxPanelSize=10,
                                       trainSplitRatio=1,
                                       ntree=500,nodeSize=1,metric="Accuracy",
                                       tuneGrid=NULL,trainControl=NULL){
  #create control and grid if left null
  if(is.null(tuneGrid)){
    message("--------------------------------------------------------------------")
    message("tuneGrid is null, creating default tuneGrid")
    tuneGrid=expand.grid(.mtry=1:10)
  }
  if(is.null(trainControl)){
    message("--------------------------------------------------------------------")
    message("trainControl is null, creating default trainControl")
    trainControl=trainControl(method = "cv", number = 5,
                              returnData=TRUE, savePredictions=TRUE, classProbs=TRUE)
  }

  #check for edge cases of every variable to make sure data is valid enough
  #hard to accurately check tune Grid & train control so just hoping users did it correctly if they want to use something custom
  if(totalPanels<1){
    message("total panels are less than one so no panels will be generated. Returning NULL")
    message("--------------------------------------------------------------------")
    return(NULL)
    }else if(totalPanels<500){
    message("--------------------------------------------------------------------")
    message("total panels are less than 500, 1000-2000 is recommend but will continue")
    }else if(totalPanels>2000){
      message("--------------------------------------------------------------------")
      message("total panels are over 2000, 1000-2000 is recommend but will continue")
    }

  if(topPanels<1){
    message("total panels are less than one so no panels will be generated. Returning NULL")
    message("--------------------------------------------------------------------")
    return(NULL)
  }else if(topPanels>totalPanels){
    message("--------------------------------------------------------------------")
    message("topPanels exceeds totalPanels, topPanels will default to totalPanels")
    topPanels=totalPanels
  }

  if(length(generationType)!=2){
    message("--------------------------------------------------------------------")
    message("generationType != 2, Need to know if seed and data is dynamic/static")
    message("setting generationType to default")
    generationType=default_generationTypes
  }else if(all(generationType==T)){ # need at least something thats dynamic
    message("--------------------------------------------------------------------")
    message("generationType is static for both data and model")
    message("Output would be identical for every iteration, reverting to default")
    message("At least one value must be dynamic")
    generationType=default_generationTypes
  }else if(!all(is.logical(generationType))){
    message("--------------------------------------------------------------------")
    message("invalid generation type found that wasnt dynamic/static")
    message("using default generation type")
    generationType=default_generationTypes
  }else if(trainSplitRatio==1 && generationType["model"]==F){
    message("--------------------------------------------------------------------")
    message("no data split is performed, so model seed is required to be dynamic")
    message("using default generation type")
    generationType=default_generationTypes
  }

  #need to add edge case checks for rfInputDataSize,maxPanelSize,ntree,nodeSize,metric="ROC"

  return(list(totalPanels=totalPanels,topPanels=topPanels,rfInputDataSize=rfInputDataSize,
              generationType=generationType,trainSplitRatio=trainSplitRatio,
              maxPanelSize=maxPanelSize,ntree=ntree,nodeSize=nodeSize,
              metric=metric,tuneGrid=tuneGrid,trainControl=trainControl))
}

#' TwoBDP_makeFrequencyParameterList
#'
#' @param amountOfFeatures An integers vector that are the size of the combinations you want.
#' @param append A character vector of features you want appended to every sub panel of features identified.
#' @param sub A character vector of features you want to use to replace the tail of panels from. Panels that are equal or smaller than sub's length will be ignored.
#' @param subGroupMethod A string that is either "all" or "slidingWindow" to determine how feature combinations are performed.
#'
#' @return List of named parameters that will be called in frequency mapping and feature combination.
#'
#' @rdname TwoBDP_makeFrequencyParameterList
#' @export
TwoBDP_makeFrequencyParameterList=function(amountOfFeatures=1:10,
                                           append=NULL,sub=NULL,
                                           subGroupMethod="all"){
  #remove neg features
  amountOfFeatures=amountOfFeatures[amountOfFeatures>0]

  if(length(amountOfFeatures)==0){message("after removing feature sizes <1, no features remained. Please provide a list of positve features. Returning NULL");return(NULL)}
  if(length(subGroupMethod)!=1){message("length(subGroupMethod)!=1 only one subGroup method allowed. returning NULL");return(NULL)}
  else if(!subGroupMethod%in%constant_panelSubgroupMethods){message("unknown subGroupMethod provided. Returning NULL");return(NULL)}

  #other edge case tests will occur later

  return(list(amountOfFeatures=amountOfFeatures,subGroupMethod=subGroupMethod,
              append=append,sub=sub))
}

#' TwoBDP_makeValidationParameterList
#'
#' @param cutoffs A named numeric vector containing the "AUC" cutoff and "pvalue" cutoff used to trim results by.
#' @param metric A string for the metric used in caret Train function.
#' @param generationType A named logical vector that determines if "data" and "model" seeds are statically or dynamically used.
#' @param trainSplitRatio A numeric value used to split the data to a training data set. The remaining values will be used to validate.
#' @param kfold A numeric value used for how many folds "kfold" method will perform.
#' @param repeats A numeric named vector that determines how many repeats are performed for each validation method. Only effects "rsbmr".
#' @param validationMethods A string vector used to identify all validation methods to run. Currently "kfold", and "rsbmr" are the only options.
#' @param tuneGrid A caret valid tuneGrid used for all validation methods.
#' @param trainControls A list of caret trainControls used for each validation method.
#'
#' @return List of named parameters that will be called in frequency mapping and feature combination.
#'
#' @rdname TwoBDP_makeFrequencyParameterList
#' @export
TwoBDP_makeValidationParameterList=function(cutoffs=c(AUC=0.8,pvalue=0.05),
                                            metric="Accuracy",generationType=default_generationTypes,
                                            trainSplitRatio=.7,
                                            kfold=10,repeats=c(kfold=0,rsbmr=10),
                                            validationMethods=c("kfold","rsbmr"),
                                            tuneGrid=NULL,trainControls=NULL){

  if(is.null(tuneGrid)){
    message("--------------------------------------------------------------------")
    message("tuneGrid is null, creating default tuneGrid")
    tuneGrid=expand.grid(.mtry=1:20)
  }
  if(is.null(trainControls)){
    message("--------------------------------------------------------------------")
    message("trainControl is null, creating default trainControl for each method")
    trainControls=list()
    if("kfold"%in%validationMethods){
      trainControls[["kfold"]]=getDefaultTrainControl("kfold")
    }
    if("rsbmr"%in%validationMethods){
      trainControls[["rsbmr"]]=getDefaultTrainControl("rsbmr")
    }
  }else{#check there is a train control for each method
    #TODO check that all custom traincontrols have: returnData=TRUE, savePredictions=TRUE, classProbs=TRUE
  }

  return(list(cutoffs=cutoffs,metric=metric,generationType=generationType,trainSplitRatio=trainSplitRatio,repeats=repeats,
              validationMethods=validationMethods,tuneGrid=tuneGrid,trainControls=trainControls))
}

#' getDefaultTrainControl
#' @description
#' Will return the default train control used for kfold or rsbmr.
#'
#' @param m The method that the train control will be built for.
#' @return a default caret trainControl.
#'
#' @importFrom caret trainControl
#'
#' @rdname defaultTrainControl
#' @export
getDefaultTrainControl=function(m){
  if(m=="kfold"){
    return(trainControl(method = "repeatedcv",repeats = 5,number = 10,
                        returnData=TRUE, savePredictions=TRUE, classProbs=TRUE))
  }else if(m=="rsbmr"){
    return(trainControl(method = "repeatedcv",repeats = 10,number = 2,
                        returnData=TRUE, savePredictions=TRUE, classProbs=TRUE))
  }
}


#' An internal function thats used to trim dataset
#' user never need to call this, can just call bdp$data after creating if needed
#' @keywords internal
get_mean_diff = function(df,boolDex){#used in getTopProbs
  idx_pos = boolDex
  idx_neg = !boolDex

  mean_pos = apply(df[idx_pos,],2,mean)
  mean_neg = apply(df[idx_neg,],2,mean)
  diff_mean = abs(mean_pos - mean_neg)
  names(diff_mean)=colnames(df)
  #^^this shouldn't be required but should be safe if apply doesnt properly keep the names col
  diff_mean = sort(diff_mean,decreasing = TRUE) #decrease makes goes 1000->1 order
  #this should be a named vector you will need to subset using names
  return(diff_mean)
}

#' TwoBDP_add functions
#'
#' @description
#' These functions were designed to allow the user to add their own dataframe of panels to skip random forest,
#' or their own vector of feature combinations to only perform validation.
#'
#' @param bdp A TwoBDP class that you want to manually add information to.
#' @param v A dataframe or vector depending if its used to skip the random forest or feature combinations
#'
#' @return TwoBDP object with the data append.
#'
#' @examples
#' bdp=createTwoBDPObject(...)
#' bdp=TwoBDP_addPanels(bdp,v=df)
#' #or
#' features=c("feat1","feat1+feat2",...)
#' bdp=TwoBDP_addFeatureCombinations(bdp,v=features)
#'
#' @rdname TwoBDP_add
#' @export
TwoBDP_addPanels=function(bdp,v){bdp$panels=v;return(bdp)}

#' @rdname TwoBDP_add
#' @export
TwoBDP_addFeatureCombinations=function(bdp,v){dbp$featureCombinations=v;return(bdp)}

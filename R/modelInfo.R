#' Generate model parameters to run key functions
#'
#' @description
#' Extract key arguments from TwoBDP object in order to run specific functions.
#'
#' @param bdp The TwoBDP object that will be parsed to generate a list of parameters to run a section of the TwoBDP algorithm.
#' @param method A String to identify what method the function will run through and compile the required arguments.
#'
#' @return A named list of required arguments needed to run a specific method.
#'
#' @rdname generateModelParametersFromTwoBDPClass
#' @export
generateModelParametersFromTwoBDPClass=function(bdp,method){
  if(is.null(method)){return(message_missingValues_function("method","generateModelParameters"))}

  modelParameters=list()
  modelParameters$data=bdp$data
  modelParameters$seedParms=bdp$seedParms
  modelParameters$pointsOfInterest=bdp$metadata[[bdp$decisionPointColumn]]

  if(method=="rf"){
    modelParameters$generationType=bdp$panelGenerationParameters$generationType
    modelParameters$metric=bdp$panelGenerationParameters$metric
    modelParameters$tuneGrid=bdp$panelGenerationParameters$tuneGrid
    modelParameters$trainControl=bdp$panelGenerationParameters$trainControl
    modelParameters$indexes=NULL

    doSplit=(bdp$panelGenerationParameters$trainSplitRatio!=1 && bdp$panelGenerationParameters$trainSplitRatio != ncol(bdp$data))
    if(doSplit){
      if(bdp$panelGenerationParameters$generationType["data"]){#dynamic
        indexes=getIndex(vec=bdp$metadata[bdp$decisionPointColumn][,1],
                         ratio=bdp$panelGenerationParameters$trainSplitRatio,
                         seed = bdp$seedParms$seed,rngKind = bdp$seedParms$RNGkind,count = bdp$panelGenerationParameters$totalPanels)
      }else{#serial
        indexes=getIndex(vec=bdp$metadata[bdp$decisionPointColumn][,1],
                         ratio=bdp$panelGenerationParameters$trainSplitRatio,
                         seed = bdp$seedParms$seed,rngKind = bdp$seedParms$RNGkind,count = 1)
      }
      modelParameters$indexes=indexes
    }
  }else if(method=="kfold" || method=="rsbmr"){
    modelParameters$cutoffs=bdp$validationParameters$cutoffs
    modelParameters$metric=bdp$validationParameters$metric
    modelParameters$generationType=bdp$validationParameters$generationType
    modelParameters$trainControl=bdp$validationParameters$trainControls[[method]]
    modelParameters$featureCombinations=bdp$featureCombinations
    modelParameters$ratio=bdp$validationParameters$trainSplitRatio
    modelParameters$repeats=bdp$validationParameters$repeats[method]

  }else{#default is rsbmr
    message("unknown TwoBDP method. Returning NULL.")
    return(NULL)
  }
  return(modelParameters)
}

#-------------------------------------------------------------------------------
#' This function handles generating the validation models
#' Shouldnt be directly called, there are maco functions that handle this
#' @importFrom pROC roc
#' @keywords internal
#' @export
getModelAndPred=function(seed,rngkind,
                  fm,data,test,
                  trainControl){
  setRNG(seed = seed,rngKind = rngkind)
  model = train(fm,data,
                 method = "glm",family = "binomial",
                 trControl = trainControl)
  pred=predict.train(object = model,newdata = test,type = "prob")
  roc_curve=roc(test$POI, pred$target)
  return(list(model=model,pred=pred,roc=roc_curve))
}
#-------------------------------------------------------------------------------
#' This function handles getting the roc information
#' Shouldnt be directly called, there are maco functions that use this
#' @keywords internal
#' @export
getValidationInformation=function(roc){
  return(c(AUC=roc$auc,
           sensitivity=mean(roc$sensitivities),
           specificity=mean(roc$specificities)))
}
#-------------------------------------------------------------------------------
#' This function handles getting the model information
#' Shouldnt be directly called, there are maco functions that use this
#' @keywords internal
#' @export
getModelInformation=function(fit,observedValues,numOfFeatures){
  r2 = with(summary(fit), 1 - deviance/null.deviance)
  r2.adjusted = 1-(1-r2)*(observedValues-1)/(observedValues-numOfFeatures-1)
  weights = summary(fit)$coefficients[,c(1,2,4)]
  p=mean(weights[,3])
  stats=apply(weights,1,function(r){paste(r,collapse = "|||")})
  stats=paste(names(stats),stats,sep = "|||")
  return(c(pvalue=p,R2=r2,R2.adj=r2.adjusted,FeatureStats=stats))
}


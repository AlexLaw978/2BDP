#' @rdname TwoBDPAutomation
#' @export
frequencyMap=function(bdp){
  if(is.null(bdp$panels)){
    message("Missing panels to compute frequencies with - add panels or run *PanelGeneration")
    return(bdp)
  }
  if(!is.null(bdp[["frequencyMap"]])){#cant use $ has shorthand equates to longer name
    message("Frequency map already exists")
    return(bdp)
  }

  #filter for panels only
  panels=bdp$panels
  panels=panels[,grepl("Feature",colnames(panels))]
  if(ncol(panels)==0){
    message("grepl panels for Feature and found nothing.")
    message("Make sure colnames has \"Features\" for columns that will be processed.")
    return(bdp)
  }
  message("Calculating frequencies")

  freqs=apply(panels,2,function(cl,panels){
    uid=unique(cl)
    freq=c()
    for(id in uid){
      freq=c(freq,sum(cl==id))
    }
    names(freq)=uid
    return(freq)
  },panels)

  message("Mapping frequencies")
  freqMap=panels
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
  #-----------------------------------------------------------------------------handle sub and appends
  #also add note that unless the sub/append values are excluded from panel generation
  #there is a possibility that panels will have duplicated values from here
  #panel combinations should account for that when creating subsets though
  if(!is.null(bdp$frequencyMapParameters$sub) && length(bdp$frequencyMapParameters$sub)<ncol(bdp$panels)){# handle substitution
    message("Subsituting tail of features")
    v=bdp$frequencyMapParameters$sub
    panels=bdp$panels
    p_pos=ncol(panels)
    for(i in length(v):1){
      replacement=rep(v[i],nrow(panels))
      panels[,p_pos]=replacement
      p_pos=p_pos-1
    }
    bdp$panels=panels
  }
  if(!is.null(bdp$frequencyMapParameters$append)){# handle appending
    message("Appending features")
    v=bdp$frequencyMapParameters$append
    panels=bdp$panels
    sdf=data.frame()
    for(i in v){
      if(nrow(sdf)==0){sdf=data.frame(rep(i,nrow(panels)));next}
      sdf=cbind.data.frame(sdf,rep(i,nrow(panels)))
    }
    colnames(sdf)=paste("AppendFeatures",1:ncol(sdf),sep = "")
    panels=cbind.data.frame(panels,sdf)
    bdp$panels=panels
  }

  return(bdp)

}
#-------------------------------------------------------------------------------
#' @rdname TwoBDPAutomation
#' @export
PanelExpansion=function(bdp){
  if(is.null(bdp$panels)){
    message("Missing panels to compute feature combinations - add panels or run *PanelGeneration")
    return(bdp)
  }
  if(!bdp$frequencyMapParameters$subGroupMethod %in% c("all","sw")){
    message("invalid subGrouping of panels")
    return(bdp)
  }
  if(!is.null(bdp$featureCombinations)){
    message("Biomarker Combinations already created")
    return(bdp)
  }
  if(bdp$panelGenerationParameters$topPanels>nrow(bdp$panels) || bdp$panelGenerationParameters$topPanels<1){
    message("top panels exceed total panels, must be less than or equal to all panels and greater than 0")
    return(bdp)
  }
  if(!is.null(bdp$saveDirectory)){
    fl=paste(bdp$saveDirectory,"FeatureCombinations.csv",sep = "/")
    if(file.exists(fl)){
      message("found feature combination file, loading that and skipping frequency map")
      bdp$featureCombinations=read.csv(fl,row.names = 1)[,1]
      return(bdp)
    }
  }

  panels=bdp$panels
  panels=panels[,grepl("Feature",colnames(panels))]
  if(ncol(panels)==0){
    message("grepl panels for Feature and found nothing.")
    message("Make sure colnames has \"Feature\" in column names that will be processed.")
    return(bdp)
  }
  amountOfFeatures=bdp$frequencyMapParameters$amountOfFeatures
  subGroupMethod=bdp$frequencyMapParameters$subGroupMethod


  getCombinatinos=function(i,panels,amountOfFeatures,subGroupMethod){
    rw=panels[i,]
    comboList=c()
    for(j in amountOfFeatures){
      if(j>length(rw)){break}
      if(subGroupMethod=="all"){
        cdf=combn(rw,j)
        grps=as.character(apply(cdf,2,function(cl){
          cl=as.character(cl)
          cl=cl[order(cl)]
          cl=unique(cl)
          return(paste(cl,collapse = "+"))
        }))
        comboList=c(comboList,grps)
      }else if(subGroupMethod=="slidingWindow"){
        grp=rw[1:j]
        grp=grp[order(grp)]
        grp=unique(grp)
        grp=paste(grp,collapse = "+")
        comboList=c(comboList,grp)
      }
    }
    return(comboList)
  }
  message("Generating combinations")
  if(bdp$threads==1){
    featureCombinations=lapply(1:bdp$panelGenerationParameters$topPanels,getCombinatinos,
                               panels,amountOfFeatures,subGroupMethod)
  }else{
    cl=makeCluster(bdp$threads)
    featureCombinations=parLapply(cl,1:bdp$panelGenerationParameters$topPanels,getCombinatinos,
                                  panels,amountOfFeatures,subGroupMethod)
    stopCluster(cl)
  }

  message("cleaning combination results")
  featureCombinations=unlist(featureCombinations)
  finalCombos=unique(featureCombinations)
  bdp$featureCombinations=finalCombos

  if(!is.null(bdp$saveDirectory)){
    fl=paste(bdp$saveDirectory,"FeatureCombinations.csv",sep = "/")
    tmp=write.csv(bdp$featureCombinations,fl)
  }

  message("Finished Panel Expansion")
  return(bdp)
}

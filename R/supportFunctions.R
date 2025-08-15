#' setRNG managing reproducibility
#'
#' @description
#' this function is called at key parts in functions to manage reproducibility within functions.
#'
#' @param seed used to set seed in set.seed function.
#' @param rngKind used to set kind in RNGkind function.
#'
#' @return NULL
#'
#' @rdname supportFunctions_setRNG
#' @export
setRNG=function(seed=100,rngKind="L'Ecuyer-CMRG"){
  RNGkind(kind = rngKind)
  set.seed(seed)
}

#' getIndex managing index creation for train & test sets
#'
#' @description
#' this function is called at key parts in functions to manage reproducibility within functions.
#' Note that if count!=1 then random seeds will be generated using that number,
#' which are then used to split the data.
#'
#' @param vec A vector of values used to determine how the data will be split.
#' @param ratio A numeric value to know what indexes will be used in training data set.
#' @param seed A numeric value to set set the seed for reproducibility.
#' @param rngKind A string value to set set the RNGkind for reproducibility.
#' @param count An integer used to know how many indexes should be generated.
#'
#' @return A named list of indexes to subset a dataset into a train & test set.
#'
#' @importFrom caret createDataPartition
#'
#' @rdname supportFunctions_getIndex
#' @export
getIndex=function(vec,ratio=0.7,seed=100,rngKind="L'Ecuyer-CMRG",count=1){
  res=list()
  if(count==1){
    setRNG(seed,rngKind)
    res[["0"]]=c(seed=seed,
                 RNGkind=rngKind,
                 subset=createDataPartition(vec,p=ratio,list = FALSE))
  }else{
    setRNG(seed,rngKind)
    newSeeds=runif(count,0,2^25)
    for(i in newSeeds){
      i=floor(i)
      setRNG(i,rngKind)
      res[[as.character(i)]]=c(seed=i,
                     RNGkind=rngKind,
                     subset=createDataPartition(vec,p=ratio,list = FALSE))
    }
  }

  return(res)
}

#' convertFlipList
#'
#' @description
#' this function is used to convert and transform a list of vectors into a dataframe.
#'
#' @param lt A list that will be converted into a dataframe.
#'
#' @return dataframe
#'
#' @rdname supportFunctions_convertFlipList
#' @export
convertFlipList=function(lt){
  return(as.data.frame(t(as.data.frame(lt))))
}

#' createTempFileName for file managment
#'
#' @description
#' this function returns a temp file name that pads 0 before the number so the file output is uniformed and easy to manage.
#'
#' @param pre A string that will be used as the prefix for the file.
#' @param num an integer that will be used to identify the data set.
#'
#' @return String file name
#'
#' @rdname supportFunctions_createTempFileName
#' @export
createTempFileName=function(pre,num){
  #If the user somehow found something that requires over 999,999,999,999 batches of panels or validations they burn in hell
  return(
    paste(pre,sprintf("%.12d",num),".csv",sep = "")
  )
}

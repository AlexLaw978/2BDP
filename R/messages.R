#all functions here are created to attempt creating uniformed output messages
#throughout the package using a functional approach. As such no function will
#be exported as they are not ment to be called by anyone.


#' An internal message function to allow clean code
#' user never need to call this
#' @keywords internal
msg_c=function(v){
  message(cat(v))
}

#' An internal message function to allow clean code
#' user never need to call this
#' @keywords internal
message_NULL=function(){
  msg_c(
    c("Returned Null Value")
  )
  return(NULL)
}

#' An internal message function to allow clean code
#' user never need to call this
#' @keywords internal
message_missingValues_function=function(dataName,functionName,rnull=T){
  msg_c(
    c(dataName,"is required to run",functionName,"function.")
  )
  if(rnull){return(message_NULL)}
}

#' An internal message function to allow clean code
#' user never need to call this
#' @keywords internal
message_notEqualWithin=function(v1,v2,type=1,rnull=T){
  mtype=c("equal to","within")[type]
  msg_c(
    c(v1,"was not",type,v2,".",
      "\nMake sure",v1,"is",type,v2)
  )
  if(rnull){return(message_NULL)}
}

#' An internal message function to allow clean code
#' user never need to call this
#' @keywords internal
message_missMatch_dataframe=function(comp,v1,v2,rnull=T){
  msg_c(
    c(comp,"between",v1,"&",v2,
          "did not match. \nMake sure all",
          comp,"from",v1,"are in the",comp,"of",v2)
  )
  if(rnull){return(message_NULL)}
}

#' An internal message function to allow clean code
#' user never need to call this
#' @keywords internal
message_incorrectSize=function(param,sz,exsz,rnull=T){
  msg_c(
    c(param,"is expected to have a length of",exsz,".",
      "\n",param,"was had a length of",sz)
  )
  if(rnull){return(message_NULL)}
}




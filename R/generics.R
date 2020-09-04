# generics.R - DESC
# /generics.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@jrc.ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

setGeneric("performance", function(x, ...) standardGeneric("performance"))

setGeneric("mcN", function(x, ...) standardGeneric("mcN"))

find.original.name <- function(fun) {

  # 'NULL' function
  if(is.null(formals(fun)))
     if(is.null(do.call(fun, args=list())))
       return("NULL")
  
  ns <- environment(fun)
  objects <- ls(envir = ns)
  
  if(isNamespace(ns))
    name <- getNamespaceName(ns)
  else
    name <- environmentName(ns)

  for (i in objects) {
    if (identical(fun, get(i, envir = environment(fun)))) {
        return(paste(name, i, sep="::"))                   
    }
  }
  return("NULL")
}






### functions for combining FLStock object
### based on FLCore functionality, but keep additional attributes

setGeneric("combine_attr", function(stk1, stk2, ...) {
  standardGeneric("combine_attr")
})

### stk1 = FLStock, stk2 = FLStock, ...
#' @rdname combine_attr
setMethod(f = "combine_attr",
          signature = signature(stk1 = "FLStock", stk2 = "FLStock"),
          definition = function(stk1, stk2, ...) {
            
  ### list with arguments
  stks <- c(list(stk1, stk2), list(...))
  ### combine stocks with FLCore functionality
  res <- do.call(FLCore::combine, stks)
  
  ### check for additional non-standard attributes
  attrs <- setdiff(names(attributes(stk1)), c(slotNames("FLStock"), "class"))
  ### keep only attributes available in both stocks
  attrs <- intersect(attrs, names(attributes(stk2)))
  
  for (attr_i in attrs) {
    
    ### list with attributes
    attr_i_lst <- lapply(stks, attr, attr_i)
    
    if (any(is(attr(stk1, attr_i)) %in% c("FLComp", "FLStock", "FLQuants", 
                                          "FLQuant", "FLModel"))) {
      
      attr(res, attr_i) <- do.call(FLCore::combine, attr_i_lst)
      
    } else if (any(is(attr(stk1, attr_i)) %in% c("FLPar"))) {
      
      attr(res, attr_i) <- do.call(cbind2, attr_i_lst)
      ### dimnames
      itns <- unlist(lapply(lapply(attr_i_lst, dimnames), "[[", "iter"))
      if (length(unique(itns)) != length(itns)) itns <- seq(length(itns))
      dimnames(attr(res, attr_i))$iter <- itns
      
      ### otherwise simply add list
    } else {
      
      attr(res, attr_i) <- attr_i_lst
      
    }
    
  }
  
  return(res)
  
})



### subset iter and keep attributes
iter_attr <- function(object, iters, subset_attributes = TRUE) {
  
  ### subset object to iter
  res <- FLCore::iter(object, iters)
  
  if (isTRUE(subset_attributes)) {
    
    ### get default attributes of object class
    attr_def <- names(attributes(new(Class = class(object))))
    
    ### get additional attributes
    attr_new <- setdiff(names(attributes(object)), attr_def)
    
    ### subset attributes
    for (attr_i in attr_new) {
      attr(res, attr_i) <- FLCore::iter(attr(res, attr_i), iters)
    }
    
  }
  
  return(res)
  
}

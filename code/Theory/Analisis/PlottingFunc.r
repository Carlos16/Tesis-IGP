## Construct parameter set, for an initial list of parameter values
ConstructParamSet <- function(params){
    for(i in seq(2, length(params))){
        params = CombParams(params,i)
    }
    return(params)
}
##
# Assuming the same number of rows for all the list elements before i,
# add a new level of combination for the values in column i. The idea is
# to expand vertically each row until i to encompass all possible values of i
##
CombParams <- function(params, i){
    n <- length(params[[i]])
    params[[i]] <- rep(params[[i]],length(params[[1]]))
    for(j in seq(1,i-1)){
        params[[j]] <- as.vector(sapply(params[[j]], rep, n))
    }
    return(params)        
}


##
# Create a sequence by putting together subsequences between each adjacent pair of breaks for the specified length.
##
createSequence <- function(breaks, lengths){
    v <- vector()
    for(i in seq(1,length(breaks) - 1)){
        range <-  breaks[i + 1] - breaks[i]
        v <- c(v,seq(breaks[i],breaks[i + 1], range / lengths[i])[1: lengths[i] + 1])
        
    }
    return(c(breaks[1],v))
}

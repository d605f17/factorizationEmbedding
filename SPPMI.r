source('matrixFactorization.r')

createHashIJ <- function(ratingMatrix){
  nrItems <- ncol(ratingMatrix)
  hashIJ <- matrix(nrow=nrItems, ncol=nrItems)
  
  for(itemOne in 1:nrItems){
    print(itemOne)
    for(itemTwo in itemOne:nrItems){
      itemsRated <- which(ratingMatrix[, itemOne] > 3 & ratingMatrix[, itemTwo] > 3)
      hashIJ[itemOne, itemTwo] <- length(itemsRated)
    }
  }
  
  return(hashIJ)
}


calculateD <- function(hashIJ){
  nrItems <- ncol(hashIJ)
  D <- 0
  
  for(itemPair in 1:nrItems){
    D <- D + length(which(!is.na(hashIJ[itemPair,]) & hashIJ[itemPair,] != 0))
  }
  
  return(D-nrItems)
}

PMI <- function(hashIJ, D, i, j){
  if(i < j){
    top <- hashIJ[i, j]*D
    bottom <- hashIJ[i, i]*hashIJ[j, j]
  }
  else{
    top <- hashIJ[j, i]*D
    bottom <- hashIJ[i, i]*hashIJ[j, j]
  }
  
  return(log(top/bottom))
}

createSPPMI <- function(hashIJ, D){
  nrItems <- 1682
  SPPMI = matrix(nrow=nrItems, ncol=nrItems)
  k = 10
  
  for(i in 1:nrItems){
    for(j in 1:nrItems){
      SPPMI[i, j] <- max(PMI(hashIJ, D, i, j) - log(k), 0)
    }
  }
  
  write.table(SPPMI, file = paste("SPPMI.csv", sep = ""), sep = ",", row.names = F, col.names = F) #virker ikke ordenligt
}

executeSPPMI <- function(fileName){
  hashIJ <- createHashIJ(makeRatingsMatrix(fileName, 943, 1682))
  D <- calculateD(hashIJ)
  createSPPMI(hashIJ, D)
}

getSPPMI <- function(fileName){
  SPPMI <- read.csv(paste(getwd(), "/", fileName, sep = ""), header=FALSE, sep = ",")
  return(as.matrix(SPPMI))
}


source('matrixFactorization.r')

createHashIJ <- function(ratingMatrix){
  nrItems <- ncol(ratingMatrix)
  hashIJ <- matrix(nrow=nrItems, ncol=nrItems)
  
  for(itemOne in 1:nrItems){
    print(itemOne)
    for(itemTwo in 1:nrItems){
      if(itemOne > itemTwo)
        hashIJ[itemOne, itemTwo] <- hashIJ[itemTwo, itemOne] 
      else {
        itemsRated <- which(ratingMatrix[, itemOne] > 3 & ratingMatrix[, itemTwo] > 3)
        hashIJ[itemOne, itemTwo] <- length(itemsRated)
      }
    }
  }
  
  return(hashIJ)
}


calculateD <- function(hashIJ){
  nrItems <- ncol(hashIJ)
  D <- 0
  
  for(itemPair in 1:nrItems){
    D <- D + sum(hashIJ[itemPair, itemPair:nrItems])
  }
  
  return(D)
}

PMI <- function(hashIJ, D, i, j){
  hashI <- sum(hashIJ[i, ])
  hashJ <- sum(hashIJ[, j])
  top <- hashIJ[i, j]*D
  bottom <- hashI * hashJ
  
  if(bottom == 0)
    return(0)
  
  return(log(top/bottom))
}

createSPPMI <- function(hashIJ, D){
  nrItems <- 1682
  SPPMI = matrix(nrow=nrItems, ncol=nrItems)
  k = 1
  
  for(i in 1:nrItems){
    print(Sys.time())
    print(i)
    for(j in 1:nrItems){
      SPPMI[i, j] <- max(PMI(hashIJ, D, i, j) - log(k), 0)
    }
  }
  
  write.table(SPPMI, file = paste("SPPMIk_1.csv", sep = ""), sep = ",", row.names = F, col.names = F)
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


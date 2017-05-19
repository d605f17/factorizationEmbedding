source('SPPMI.r')
source('matrixFactorization.r')

solveTheta <- function(theta, userId, beta, K, alpha, lambda, ratingMatrix, l){
  firstSum <- matrix(0, nrow = K, ncol = K)
  secondSum <- matrix(0, nrow = K, ncol = 1)
  ratedItems <- which(!is.na(ratingMatrix[userId, ]))
  
  for(item in ratedItems){
    cui <- l
    rating <- ratingMatrix[userId, item]
    if(!is.na(rating)){
    cui <- l * (1 + rating * alpha)
    }
    firstSum <- firstSum + (cui * as.matrix(beta[item, ]) %*% beta[item, ])
    secondSum <- secondSum + (cui * rating * as.matrix(beta[item, ]))
  }
  
  firstHalf <- solve(firstSum + (lambda * diag(K)))
  result <- firstHalf %*% secondSum

  return(t(result))
}

solveBeta <- function(theta, itemId, gamma, K, alpha, SPPMI, lambda,  ratingMatrix, wi, cj, l){
  firstSum <- matrix(0, nrow = K, ncol = K)
  secondSum <- matrix(0, nrow = K, ncol = K)
  thirdSum <- matrix(0, nrow = K, ncol = 1)
  fourthSum <- matrix(0, nrow = K, ncol = 1)
  ratedUsers <- which(!is.na(ratingMatrix[, itemId]))
  nonZeroContexts <- which(SPPMI[itemId, ] != 0)
  
  for(user in ratedUsers){
    rating <- ratingMatrix[user, itemId]
    cui <- l
    
    if(!is.na(rating)){
    cui <- cui + l * rating * alpha
    }
    firstSum <- firstSum + (cui * as.matrix(theta[user, ]) %*% theta[user, ])
    thirdSum <- thirdSum + (cui * rating * as.matrix(theta[user, ]))
  }
  
  for(context in 1:length(nonZeroContexts)){
    if(length(nonZeroContexts) != 0){
      secondSum <- secondSum + as.matrix(gamma[context, ]) %*% gamma[context, ]
      fourthSum <- fourthSum + ((SPPMI[itemId, context] - wi[itemId, ] - cj[context, ]) * as.matrix(gamma[context, ]))
    }else{
      secondSum <- secondSum + 0
      fourthSum <- fourthSum + (-wi[itemId, ])
    }
  }
  
  secondSum <- secondSum + as.numeric(lambda) * diag(K)
  firstHalf <- solve(firstSum + secondSum)
  secondHalf <- thirdSum + fourthSum
  result <- firstHalf %*% secondHalf
  
  return(t(result))
};

solveGamma <- function(itemId, beta, K, SPPMI, wi, cj, lambda){
  firstSum <- matrix(0, nrow = K, ncol = K)
  secondSum <- matrix(0, nrow = K, ncol = 1)
  nonZeroContexts <- which(SPPMI[itemId, ] != 0)
  
  for(context in 1:length(nonZeroContexts)){
    if(length(nonZeroContexts) != 0){
      firstSum <- firstSum + (as.matrix(beta[context, ]) %*% beta[context, ])
      secondSum <- secondSum + (SPPMI[itemId, context] - wi[context, ] - cj[itemId, ] * as.matrix(beta[context, ]))
    }else{
      firstSum <- firstSum
      secondSum <- secondSum + (- cj[itemId, ])
    }
  }
  
  firstSum <- solve(firstSum + lambda * diag(K))
  result <- firstSum %*% secondSum
  
  return(t(result))
}

solveWi <- function(beta, itemId, gamma, SPPMI, cj){
  firstSum <- 0
  nonZeroContexts <- which(SPPMI[itemId, ] != 0)
  
  for(context in 1:length(nonZeroContexts)){
    if(length(nonZeroContexts) != 0){
      firstSum <- firstSum + (SPPMI[itemId, context] - (beta[itemId, ] %*% 
                  as.matrix(gamma[context, ])) - cj[context, ])  
    }else{
      firstSum <- firstSum
    }
  }
  
  if(length(nonZeroContexts) != 0){
    return((1/length(nonZeroContexts)) * as.numeric(firstSum))
  }else{
    return(0)
  }
}

solveCj <- function(beta, itemId, gamma, SPPMI, wi){
  firstSum <- 0
  nonZeroContexts <- which(SPPMI[itemId, ] != 0)
  
  for(context in 1:length(nonZeroContexts)){
    if(length(nonZeroContexts) != 0){
      firstSum <- firstSum + (SPPMI[context] - (beta[itemId, ] %*% as.matrix(gamma[context, ])) - wi[context, ])
    }else{
      firstSum <- firstSum
    }
  }
  
  if(length(nonZeroContexts) != 0){
    return((1/length(nonZeroContexts)) * as.numeric(firstSum))
  }else{
    return(0)
  }
}

solveLco <- function(gamma, theta, beta, SPPMI, wi, cj, lambda, ratingMatrix, alpha, numberOfUsers, numberOfItems){
  firstSum <- 0
  secondSum <- 0
  thirdSum <- 0
  fourthSum <- 0
  fifthSum <- 0
  
  for(user in 1:numberOfUsers){
    userRatedItems <- which(!is.na(ratingMatrix[user, ]))
    for(item in 1:length(userRatedItems)){
      cui <- 1 + ratingMatrix[user, userRatedItems[item]] * alpha
      firstSum <- firstSum + (cui * (ratingMatrix[user, userRatedItems[item]] - theta[user, ] %*% as.matrix(beta[userRatedItems[item], ]))^2)
    }
  }

  for(item in 1:nrow(SPPMI)){
    nonZeroContexts <- which(SPPMI[item, ] != 0)
    for(context in 1:length(nonZeroContexts)){
      if(length(nonZeroContexts) != 0){
        secondSum <- secondSum + (SPPMI[item, nonZeroContexts[context]] - (beta[item, ] %*% as.matrix(gamma[nonZeroContexts[context], ])) - wi - cj)^2  
      }else{
        secondSum <- secondSum + (-wi - cj)^2
      }
    }
  }
  
  for(user in 1:numberOfUsers){
    thirdSum <- thirdSum + (norm(as.matrix(theta[user, ]), type="f"))^2
  }
  thirdSum <- lambda * thirdSum
  
  for(item in 1:numberOfItems){
    fourthSum <- fourthSum + (norm(as.matrix(beta[item, ]), type="f"))^2
    fifthSum <- fifthSum + (norm(as.matrix(gamma[item, ]), type="f"))^2
  }
  fourthSum <- lambda * fourthSum
  fifthSum <- lambda * fifthSum
  
  result <- firstSum + secondSum + thirdSum + fourthSum + fifthSum
  return(result)
}

calculateNDCG <- function(userId, predictions, M, testData, ratingMatrix){
  pi <- order(predictions[userId, ], decreasing = TRUE)
  rel <- order(ratingMatrix[userId, ], decreasing = TRUE)
  consideredItems <- 0
  DCG <- 0
  IDCG <- 0
  knownItems <- as.matrix(testData[which(testData[, 1] == userId), 2])
  
  for(item in pi){
    if(item %in% knownItems){
      consideredItems <- consideredItems + 1
      consumed <- 0
      rating <- testData[which(testData[, 1] == userId & testData[, 2] == item), 3]
        if(rating >= 4){
          consumed <- 1
        }
      DCG <- DCG + ((2^consumed) - 1)/(log(item + 1))
      if(consideredItems >= M){
        break
      }
    }
  }
  
  consideredItems <- 0
  for(item in rel){
    consideredItems <- consideredItems + 1
    rating <- testData[which(testData[, 1] == userId & testData[, 2] == item), 3]
    consumed <- 0
    if(rating >= 4){
      consumed <- 1
    }
    IDCG <- IDCG + ((2^consumed) - 1)/log(item + 1)
    if(consideredItems >= M){
      break
    }
  }
  
  if(is.na(DCG/IDCG)){
    return(0)
  }
  return(DCG/IDCG)
}


train <- function(filename){
  numberOfItems <- 1682
  numberOfUsers <- 943
  K <- 100
  beta <- matrix(runif(numberOfItems * K, 0, 1), nrow = numberOfItems, ncol = K)
  gamma <- matrix(runif(numberOfItems * K, 0, 1), nrow = numberOfItems, ncol = K)
  theta <- matrix(runif(numberOfUsers * K, 0, 1), nrow = numberOfUsers, ncol = K)
  wi <- matrix(0, nrow = numberOfItems, ncol = 1)
  cj <- matrix(0, nrow = numberOfItems, ncol = 1)
  yi <- 0
  l <- 1
  alpha <- 40
  lambda <- l * 0.5
  totalNDDCG <- 0
  
  trainData <- read_delim(paste(getwd(), "/ml-100k/", filename, ".base", sep = ""),
                          "\t", escape_double = FALSE, trim_ws = TRUE, 
                          col_names = c("userId", "movieId", "rating", "timestamp"),
                          col_types = cols(
                            userId = col_integer(),
                            movieId = col_integer(),
                            rating = col_integer(),
                            timestamp = col_integer()
                          )
  );
  trainData <<- as.matrix(trainData)
  SPPMI <- as.matrix(getSPPMI("SPPMI_k1.csv"))
  
  print('Training started')
  startTime <- Sys.time()
  for(step in 1:3){
    Lco <- 0
    
    for(user in 1:numberOfUsers){
      theta[user, ] <- solveTheta(theta, user, beta, K, alpha, lambda, ratingMatrix, l)
    }
    
    print("SKIFT!")
    for(item in 1:numberOfItems){
      beta[item, ] <- solveBeta(theta, item, gamma, K, alpha, SPPMI, lambda, ratingMatrix, wi, cj, l)
      gamma[item, ] <- solveGamma(item, beta, K, SPPMI, wi, cj, lambda)
      wi[item, ] <- solveWi(beta, item, gamma, SPPMI, cj)
      cj[item, ] <- solveCj(beta, item, gamma, SPPMI, wi)
    }
    predictions <- theta %*% t(beta)

    for(user in 1:numberOfUsers){
      totalNDDCG <- totalNDDCG + calculateNDCG(user, predictions, 10, trainData, ratingMatrix)
    }
    print(totalNDDCG/numberOfUsers)
    totalNDDCG <- 0
  }
  print(startTime)
  print(Sys.time())
}
predictions <- train("ua")





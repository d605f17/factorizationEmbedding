source('SPPMI.r')
#http://stackoverflow.com/questions/15170399/changing-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work

sumVectorNorm <- function(id, matrix){
  result <- 0
  elements <- matrix[id, ]
  elements <- norm(matrix, type = "f")
  
  for(element in  elements){
    result <- result + element
  }
  return(result)
}

train <- function(filename){
  numberOfItems <- 1682
  numberOfUsers <- 943
  K <- 20
  beta <- matrix(runif(numberOfItems * K, 0, 1), nrow = numberOfItems, ncol = K)
  gamma <- matrix(runif(numberOfItems * K, 0, 1), nrow = numberOfItems, ncol = K)
  theta <- matrix(runif(numberOfUsers * K, 0, 1), nrow = numberOfUsers, ncol = K)
  yi <- 0
  alpha <- 40
  lambda <- 0.05
  
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
  SPPMI <- as.matrix(getSPPMI("SPPMI.csv"))
  
  print('Training started')
  for(step in 1:5000){
    Lco <- 0
    
    for(row in 1:nrow(trainData)){
      userId <- as.numeric(trainData[row, 1])
      itemId <- as.numeric(trainData[row, 2])
      rating <- as.numeric(trainData[row, 3])
      
      yui <- rating
      cui <- 1 + alpha * yui
      wi <- 0
      cj <- 0
      
      Lco <- as.numeric(Lco + (cui * (yui - t(theta[userId, ]) %*% beta[itemId, ])^2 +
                               (SPPMI[userId, itemId] - t(beta[itemId, ]) %*% gamma[itemId, ] - wi - cj)^2 +
                                 lambda * (theta[userId, ] * sumVectorNorm(userId ,theta)^2) +
                                 lambda * (beta[itemId, ] * sumVectorNorm(itemId, beta)^2) +
                                 lambda * (gamma[itemId, ] * sumVectorNorm(itemId, gamma)^2))
                                 )
      theta[userId, ] <- (cui * beta[itemId, ] %*% t(beta[itemId, ]) + lambda * numberOfItems)^-1 * (cui * yui * beta[itemId, ])
      beta[itemId, ] <- ((cui * theta[userId, ] %*% t(theta[userId, ])) + (gamma[itemId, ] %*% t(gamma[itemId, ]) + lambda * numberOfItems))^-1 *
        (cui * yui * theta[userId, ] + (SPPMI[itemId, itemId] - wi - cj) * gamma[itemId, ])
      gamma[itemId, ] <- (beta[itemId, ] %*% beta[itemId, ] + lambda * numberOfItems)^-1 * ((SPPMI[itemId, itemId] - wi - cj) * beta[itemId, ])
    }
  }
}





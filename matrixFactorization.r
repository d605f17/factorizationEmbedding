install.packages('readr')
install.packages('tidyverse')

makeRatingsMatrix <- function(filename, numberOfUsers, numberOfItems) {
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
  trainData <- as.matrix(trainData)
  
  ratingsMatrix <- matrix(nrow = numberOfUsers, ncol = numberOfItems)
  for(row in 1:nrow(trainData)){
    ratingsMatrix[as.numeric(trainData[row, 1]), as.numeric(trainData[row, 2])] <- as.numeric(trainData[row, 3])
  }
  return(ratingsMatrix)
}

ratingMatrix<-makeRatingsMatrix('ua', 943, 1682)


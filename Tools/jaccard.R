# Jaccar Index
# Your dataset
library(dplyr)
df <- data.frame(t(data.frame(c1=rnorm(100),
                              c2=rnorm(100),
                              c3=rnorm(100),
                              c4=rnorm(100),
                              c5=rnorm(100),
                              c6=rnorm(100))))

df[df > 0] <- 1
df[df <= 0] <- 0
df
# Function returns the Jaccard index and Jaccard distance
# Parameters:
# 1. df, dataframe of interest
# 2. margin, axis in which the apply function is meant to move along
jaccard <- function(df, margin=1) {
  if (margin == 1 | margin == 2) {
    M_00 <- apply(df, margin, sum) == 0
    M_11 <- apply(df, margin, sum) == 2
    if (margin == 1) {
      df <- df[!M_00, ]
      JSim <- sum(M_11) / nrow(df)
    } else {
      df <- df[, !M_00]
      JSim <- sum(M_11) / length(df)
    }
    JDist <- 1 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}
jaccard(df[1:2,], margin=2)
jaccard_per_row <- function(df, margin=1){
  key_pairs <- expand.grid(row.names(df), row.names(df))
  results <- t(apply(key_pairs, 1, function(row) jaccard(df[c(row[1], row[2]),], margin=margin)))
  key_pair <- key_pairs %>% mutate(pair = paste(Var1,"_",Var2,sep=""))
  results <- data.frame(results)
  row.names(results) <- key_pair$pair
  results
}
jaccard_per_row(df, margin=2)
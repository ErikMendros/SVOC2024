library(microbenchmark)
library(ggplot2)


#
#
#  Algorithms to count the number of k-simplices containing the origin
#
#


#
#* Naive approach--------------------------------------------------------------
#


#Naive approach to count K-simplices
#INPUT:  bivariate dataset X (as.matrix); integer K >= 3
#OUTPUT: number of K-simplices containing the origin / number of all K-simplices.     

Naive <- function(X,K){
  n <- dim(X)[1]
  index <- combn((1:n), K)
  L <- dim(index)[2]
  cnt <- 0
  for (i in (1:L)){
    hull <- as.vector(convhulln(rbind(c(0,0),X[index[,i],])))
    cnt <- cnt + !(1 %in% hull)
  }
  return(cnt)
}


#
#* CountSD_K Implementation--------------------------------------------------------------
#


#The new O(n log n) algorithm CountSD_K
#INPUT:  bivariate dataset X (as.matrix); integer K >= 3
#OUTPUT: number of K-simplices containing the origin / number of all K-simplices.     

CountSD_K = function(X,K){
  #Transforming data X to ordered 01-sequence s     <- O(n*log(n))
  X <- cbind(X, as.numeric(X[,1] > 0))
  s <- X[,3][order(X[,2] / X[,1])]
  n <- length(s)

  cts0 <- rep(0,2+(K*(K+1)/2))
  cts0[1] <- 1
  cts1 <- rep(0,2+(K*(K+1)/2))
  cts1[1] <- 1
  G1.green <- c()
  G1.yellow <- c(2)
  G0.green <- c()
  G0.yellow <- c(1)
  var1 <- 2 #temporary variable to help initialize the indices
  
  #Initialization of indices                        <- O(K^2)
  for(k in (1:K)){
    if(k %% 2 != 0){ 
      G1.green <- c(G1.green, var1 + (1:k))
      G0.green <- c(G0.green, var1 + (1:k) - k)
      var1 <- G1.green[length(G1.green)] + k + 1
    } else {
      var2 <- G1.yellow[length(G1.yellow)] + k - 1
      G1.yellow <- c(G1.yellow, var2 + (1:k))
      G0.yellow <- c(G0.yellow, var2 + (1:k) - k)
    }
  }
  #Updating the counts                              <- O(n*K^2)  
  for(i in 1:n){
    if(s[i]==0){
      # if s_i is 0 
      temp <- cts0[G1.green] + cts1[G0.green]
      cts1[G1.yellow] <- cts1[G1.yellow] + cts0[G0.yellow]
      cts0[G1.green] <- temp} 
    else{
      # if s_i is 1  
      temp <- cts0[G1.yellow] + cts1[G0.yellow]
      cts1[G1.green] <- cts1[G1.green] + cts0[G0.green]
      cts0[G1.yellow] <- temp
    }
  }
  #Sumation of all simple sequences of length K    <- O(K) 
  res = sum(cts0[((length(cts0)-K):(length(cts0)-1))]+
              cts1[((length(cts0)-K):(length(cts0)-1))])
  return(choose(n,K) - res)
}

n=15
K=5
X = matrix(rnorm(n*2),ncol=2)
plot(X)
CountSD_K(X,K)/choose(n,K)
Naive(X,K)/choose(n,K)

#
#* C++ implementation of CountSD_K for K=3 --------------------------------------------------------------
#


#The function SD below utilizes the package KDepth which was written in C++ by Stanislav Nagy.
#The principle of the algorithm is the same as CountSD_K. The only difference is that instead 
#of counting triangles not containing the origin, we count those contaning it. That is, we count the 
#occurences of subsequences 101 and 010 within transformed ordered 01-sequence s obtained from X.

#The package KDepth can be downloaded here 
#https://github.com/NagyStanislav/Kdepth
install.packages("KDepth_0.1.5.tar.gz", repos=NULL)
library(KDepth)

SD = function(X){
  x = X[order(X[,2] / X[,1]),1]
  res = sum(KSign(sign(x),3)[3,])/choose(length(x),3)
  return(res)
}


#
#* Comparisons with the best R implementations --------------------------------------------------------------
#


library(mrfDepth)
library(ddalpha)

set.seed(1)
n=1000
X = matrix(rnorm(n*2),ncol=2)
plot(X)

sdepth(X, matrix(c(0,0), nrow = 1))$depthZ
depth.simplicial(c(0,0), X, exact=T)
SD(X)
bm = bench::mark(SD(X),
                 sdepth(X, matrix(c(0,0), nrow = 1))$depthZ,
                 depth.simplicial(c(0,0), X, exact = T), 
                 check=FALSE, iterations = 100, filter_gc = TRUE)
autoplot(bm)



set.seed(1)
n=5000
X = matrix(rnorm(n*2),ncol=2)
plot(X)


sdepth(X, matrix(c(0,0), nrow = 1))$depthZ
depth.simplicial(c(0,0), X, exact=T) #depth.simplicial only works for medium-sized datasets (up to n=2000)
SD(X)
bm = bench::mark(SD(X),
                 sdepth(X, matrix(c(0,0), nrow = 1))$depthZ,
                 check=FALSE, iterations = 100, filter_gc = TRUE)
autoplot(bm)





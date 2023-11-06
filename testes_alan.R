## PACOTES ##

library(dplyr)

## CONSTRUÇÃO DOS DADOS ##

# Dados menores #

set.seed(543357)
dados=as.data.frame(rnorm(100))
dados$x1 <- 1
dados <- dados%>%rename(y=`rnorm(100)`)
dados$x2 <- rpois(100, 2)


# Dados maiores #

set.seed(543357)
dados=as.data.frame(rnorm(100000000))
dados$x1 <- 1
dados <- dados%>%rename(y=`rnorm(1e+08)`)
dados$x2 <- rpois(100000000, 2)
dados$x3 <- rpois(100000000, 5)

## FUNÇÕES DO R BASE PARA COMPARAÇÃO ##

t.test(dados$y)

model1 <- lm(data=dados, formula=y~x2)
summary(model1)

## MINHAS FUNÇÕES ##

regr <- function(x, y, int=TRUE){
  n <- length(y)
  if (int){
    x <- cbind(rep(1, n), x)
  }
  else{
    x <- matrix(x, nrow=n, ncol=1)
  }
  Beta <- solve(t(x)%*%x)%*%t(x)%*%y
  erro <- y-x%*%Beta
  sigma2 <- as.numeric((t(erro)%*%erro)/(n-ncol(x)))
  var_beta <- as.numeric(diag(solve(t(x)%*%x)*sigma2))
  t <- Beta/sqrt(var_beta)
  print(Beta)
  print(t)
}

regr2 <- function(x, y, int=TRUE){
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- nrow(y)
  if (int){
    x <- cbind(rep(1, n), x)
  }
  Beta <- solve(t(x)%*%x)%*%t(x)%*%y
  erro <- y-x%*%Beta
  sigma2 <- as.numeric((t(erro)%*%erro)/(n-ncol(x)))
  var_beta <- as.numeric(diag(solve(t(x)%*%x)*sigma2))
  t <- Beta/sqrt(var_beta)
  print(Beta)
  print(t)
}

regr3 <- function(x=c(dados$x1, dados$x2), y, int=TRUE){
  y <- as.matrix(y)
  n <- nrow(y)
  x <- matrix(x, n)
  if (int){
    x <- cbind(rep(1, n), x)
  }
  Beta <- solve(t(x)%*%x)%*%t(x)%*%y
  erro <- y-x%*%Beta
  sigma2 <- as.numeric((t(erro)%*%erro)/(n-ncol(x)))
  var_beta <- as.numeric(diag(solve(t(x)%*%x)*sigma2))
  t <- Beta/sqrt(var_beta)
  print(Beta)
  print(t)
}

## TESTES ##

# sem intercepto #

startTime1 <- Sys.time()
regr2(x=dados["x1"], y=dados["y"], int=F)
endTime1 <- Sys.time()
endTime1-startTime1

startTime2 <- Sys.time()
regr2(x=dados$x1, y=dados$y, int=F)
endTime2 <- Sys.time()
endTime2-startTime2

startTime3 <- Sys.time()
regr(x=dados$x1, y=dados$y, int=F)
endTime3 <- Sys.time()
endTime3-startTime3

startTime4 <- Sys.time()
regr3(x=dados$x1, y=dados$y, int=F)
endTime4 <- Sys.time()
endTime4-startTime4

startTime_lm1 <- Sys.time()
model1 <- lm(data=dados, formula=y~0+x2)
startTime_lm1 <- Sys.time()
summary(model1)
endTime_lm1-startTime_lm1

# com intercepto #

startTime5 <- Sys.time()
regr2(x=dados["x2"], y=dados["y"], int=T)
endTime5 <- Sys.time()
endTime5-startTime5

startTime6 <- Sys.time()
regr2(x=dados$x2, y=dados$y, int=T)
endTime6 <- Sys.time()
endTime6-startTime6

startTime7 <- Sys.time()
regr(x=dados$x2, y=dados$y, int=T)
endTime7 <- Sys.time()
endTime7-startTime7

startTime8 <- Sys.time()
regr3(x=dados$x2, y=dados$y, int=T)
endTime8 <- Sys.time()
endTime8-startTime8

startTime_lm2 <- Sys.time()
model1 <- lm(data=dados, formula=y~x2)
endTime_lm2 <- Sys.time()
summary(model1)
endTime_lm2-startTime_lm2

## sem intercepto - mais de 1 var ##

startTime1 <- Sys.time()
regr2(x=dados[c("x1", "x2")], y=dados["y"], int=F)
endTime1 <- Sys.time()
endTime1-startTime1

startTime4 <- Sys.time()
regr3(x=c(dados$x1, dados$x2), y=dados$y, int=F)
endTime4 <- Sys.time()
endTime4-startTime4

## FORMULAS ##

teste <- function(df, formula){
  x <- formula
}

teste(dados, y~x1)

#na fórmula, [[1]] é o ~, [[2]] é o y e [[3]] é o x1+x2
#no elemento explicativo ([[3]]), length() funciona para pegar a quantidade de variáveis (+1, por causa do +)
#    e [[3]][[1]] será o +, [[3]][[2]] será o x1, [[3]][[3]], será o x2, etc...

#mas e se tiver vários x?

regr3 <- function(df, x, y, int=TRUE){
  y <- df$y
  x <- sapply(parse(text=paste(deparse(substitute(df)), "$", x, sep="")), eval)
  n <- length(y)
  if (int){
    x <- cbind(rep(1, n), x)
  }
  Beta <- solve(t(x)%*%x)%*%t(x)%*%y
  erro <- y-x%*%Beta
  sigma2 <- as.numeric((t(erro)%*%erro)/(n-ncol(x)))
  var_beta <- as.numeric(diag(solve(t(x)%*%x)*sigma2))
  t <- Beta/sqrt(var_beta)
  print(Beta)
  print(t)
}

regr3 <- function(x, y, int=TRUE){
  n <- length(y)
  if (int){
    x <- cbind(rep(1, n), x)
  }
  else{
    x <- matrix(x, nrow=n)
  }
  Beta <- solve(t(x)%*%x)%*%t(x)%*%y
  erro <- y-x%*%Beta
  sigma2 <- as.numeric((t(erro)%*%erro)/(n-ncol(x)))
  var_beta <- as.numeric(diag(solve(t(x)%*%x)*sigma2))
  t <- Beta/sqrt(var_beta)
  print(Beta)
  print(t)
}



startTime5 <- Sys.time()
regr3(dados, x="x2", y="y", int=T)
endTime5 <- Sys.time()
endTime5-startTime5

startTime6 <- Sys.time()
regr2(x=c(dados$x1, dados$x2), y=dados$y, int=T)
endTime6 <- Sys.time()
endTime6-startTime6

startTime7 <- Sys.time()
regr2(x=dados["x2"], y=dados["y"], int=T)
endTime7 <- Sys.time()
endTime7-startTime7

startTime6 <- Sys.time()
regr3(x=c(dados$x1, dados$x2), y=dados$y, int=T)
endTime6 <- Sys.time()
endTime6-startTime6

library(DescTools)
regr4 <- function(formula, data, int=TRUE){ #formulas
  fo = ParseFormula(formula, data)
  y <- as.matrix(fo$lhs$mf)
  x <- as.matrix(fo$rhs$mf)
  n <- nrow(y)
  if (int){
    x <- cbind(rep(1, n), x)
  }
  Beta <- solve(t(x)%*%x)%*%t(x)%*%y
  erro <- y-x%*%Beta
  sigma2 <- as.numeric((t(erro)%*%erro)/(n-ncol(x)))
  var_beta <- as.numeric(diag(solve(t(x)%*%x)*sigma2))
  t <- Beta/sqrt(var_beta)
  print(Beta)
  print(t)
}

## sem intercepto - mais de 1 var ##

startTime1 <- Sys.time()
regr2(x=dados[c("x1", "x2")], y=dados["y"], int=F)
endTime1 <- Sys.time()
endTime1-startTime1

startTime2 <- Sys.time()
regr3(x=c(dados$x1, dados$x2), y=dados$y, int=F)
endTime2 <- Sys.time()
endTime2-startTime2


startTime3 <- Sys.time()
regr4(formula=y~x1+x2, data=dados, int=F)
endTime3 <- Sys.time()
endTime3-startTime3

bw <- ggwr.sel(PctBach~PctEld+TotPop90,
               data=georgia,
               family=gaussian(),
               coords=matrix(georgia$X, georgia$Y, nrow(georgia)),
               longlat=TRUE)

bw <- ggwr.sel(PctBach~PctEld+TotPop90,
               data=georgia)


ggwr(formula=PctBach~PctEld+TotPop90, data = georgia, bandwidth=136)

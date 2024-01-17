## CARREGANDO DADOS DA GEORGIA ##

if (!require(spgwr)){
  install.packages("spgwr")
  library(spgwr)
}

data(georgia)
georgia <- as.data.frame(gSRDF)

## PADRONIZANDO DADOS ##

georgia_std <- georgia

#TotPop90
georgia_std$TotPop90 <- (georgia_std$TotPop90-
                           mean(georgia_std$TotPop90))/
  sd(georgia_std$TotPop90)
#PctEld
georgia_std$PctEld <- (georgia_std$PctEld-
                         mean(georgia_std$PctEld))/
  sd(georgia_std$PctEld)
#PctFB
georgia_std$PctFB <- (georgia_std$PctFB-
                        mean(georgia_std$PctFB))/
  sd(georgia_std$PctFB)
#PctBlack
georgia_std$PctBlack <- (georgia_std$PctBlack-
                           mean(georgia_std$PctBlack))/
  sd(georgia_std$PctBlack)

## DISCRETIZANDO RESPOSTA ##

georgia_nb_std <- georgia_std

georgia_nb_std$PctBach <- as.integer(georgia_nb_std$PctBach)

## TESTES ##

## gwr

#mgwnbr
startTime1 <- Sys.time()
mgwnbr3(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="fixed_bsq",
        bandwidth="cv", model="gaussian", mgwr=FALSE)
endTime1 <- Sys.time()
endTime1-startTime1

## spgwr
startTime2 <- Sys.time()
gwr.sel(formula=PctBach~PctBlack+PctFB+TotPop90+PctEld, data=georgia_std, coords=as.matrix(georgia_std[, c("X", "Y")]), adapt=FALSE,
        method = "cv")
endTime2 <- Sys.time()
endTime2-startTime2
#, longlat=NULL
#, gweight=gwr.Gauss
#, weights

startTime3 <- Sys.time()
out <- gwr(formula=PctBach~PctBlack+PctFB+TotPop90+PctEld, data=georgia_std, coords=as.matrix(georgia_std[, c("X", "Y")]),
    bandwidth=128769.9)
out
endTime3 <- Sys.time()
endTime3-startTime3
#, fit.points


#mgwnbr demora uns 50 segundos, spgwr (somando tempo das duas funções) demora 1 segundo

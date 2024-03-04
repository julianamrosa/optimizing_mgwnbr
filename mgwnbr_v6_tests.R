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

#versão 5
startTime <- Sys.time()
mgwnbr5(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#7.5 segs
#7.9 segs
#8 segs

#versão 6
startTime <- Sys.time()
mgwnbr6(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#8.2 segs
#9.1 segs
#8.8 segs


## mgwr

#versão 5
startTime <- Sys.time()
mgwnbr5(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian")
endTime <- Sys.time()
endTime-startTime
#

#versão 6
startTime <- Sys.time()
mgwnbr6(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian")
endTime <- Sys.time()
endTime-startTime
#

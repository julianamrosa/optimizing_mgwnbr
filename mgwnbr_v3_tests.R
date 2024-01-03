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

#com as.numeric() (original)
startTime <- Sys.time()
mgwnbr2(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#49.6 segs
#1.12 mins
#57.5 segs

#com as.vector()
startTime <- Sys.time()
mgwnbr3(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#43.4 segs
#57.0 segs
#49.1 segs

## mgwr

#com as.numeric() (original)
startTime <- Sys.time()
mgwnbr2(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian")
endTime <- Sys.time()
endTime-startTime
#30.18 mins
#11.95 mins
#14.34 mins
#12.91

#com as.vector()
startTime <- Sys.time()
mgwnbr3(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
        lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq",
        bandwidth="cv", model="gaussian")
endTime <- Sys.time()
endTime-startTime
#12.61 mins
#13.35 mins
#12.06 mins
#12.96


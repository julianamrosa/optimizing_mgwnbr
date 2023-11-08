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

## EXEMPLO ##

startTime <- Sys.time()
mgwnbr(DATA=georgia_std, YVAR="PctBach",
       XVAR=c("PctBlack", "PctFB", "TotPop90", "PctEld"),
       LAT="Y", LONG="X", GLOBALMIN="NO", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="GAUSSIAN", MGWR="NO")
endTime <- Sys.time()
endTime-startTime
#1.46 mins (original)

startTime <- Sys.time()
mgwnbr2(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld,
       LAT="Y", LONG="X", GLOBALMIN="NO", METHOD="ADAPTIVE_BSQ",
       BANDWIDTH="CV", MODEL="GAUSSIAN", MGWR="NO")
endTime <- Sys.time()
endTime-startTime
#1.45 mins

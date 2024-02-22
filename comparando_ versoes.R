library(mgwnbr)
data(georgia)
for (var in c("TotPop90", "PctRural",
              "PctEld", "PctFB", "PctPov", "PctBlack")){
  georgia[, var] <- as.data.frame(scale(georgia[, var]))
}
georgia$PctBach <- as.integer(georgia$PctBach)


#Teste 1 do artigo - gwnbr

startTime <- Sys.time()
mgwnbr(DATA=georgia, YVAR="PctBach", XVAR=c("TotPop90", "PctRural",
                                                   "PctEld", "PctFB", "PctPov", "PctBlack"), LAT="Y", LONG="X", GLOBALMIN="no",
       METHOD="ADAPTIVE_BSQ", BANDWIDTH="aic", MODEL="negbin", MGWR="no")
endTime <- Sys.time()
endTime-startTime
#1.2 mins

#Versão 2

startTime <- Sys.time()
mgwnbr2(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
       method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#2.5 mins

#Versão 3

startTime <- Sys.time()
mgwnbr3(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#1.7 mins

#Versão 4

startTime <- Sys.time()
mgwnbr4(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#19 segs

#Versão 5 - cran

startTime <- Sys.time()
mgwnbr5(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime

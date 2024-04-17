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
#1.57 mins --> 94 segs
#1.59 mins --> 95 segs

#Versão 2

startTime <- Sys.time()
mgwnbr2(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
       method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#1.56 mins --> 94 segs
#1.58 mins --> 95 segs

#Versão 3

startTime <- Sys.time()
mgwnbr3(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#1.57 mins --> 94 segs
#1.57 mins --> 94 segs

#Versão 4

startTime <- Sys.time()
mgwnbr4(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#25 segs
#25 segs

#Versão 5

startTime <- Sys.time()
mgwnbr5(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#25 segs
#23 segs

#Versão 6

startTime <- Sys.time()
mgwnbr6(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime
#26 segs
#26 segs

#################

#Teste 2 do artigo - mgwnbr

startTime <- Sys.time()
mgwnbr(DATA=georgia, YVAR="PctBach", XVAR=c("TotPop90", "PctRural",
                                            "PctEld", "PctFB", "PctPov", "PctBlack"), LAT="Y", LONG="X", GLOBALMIN="no",
       METHOD="ADAPTIVE_BSQ", BANDWIDTH="aic", MODEL="negbin")
endTime <- Sys.time()
endTime-startTime


#Versão 2

startTime <- Sys.time()
mgwnbr2(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin")
endTime <- Sys.time()
endTime-startTime


#Versão 3

startTime <- Sys.time()
mgwnbr3(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin")
endTime <- Sys.time()
endTime-startTime


#Versão 4

startTime <- Sys.time()
mgwnbr4(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin")
endTime <- Sys.time()
endTime-startTime
#5.96 --> 357 segs

#Versão 5

startTime <- Sys.time()
mgwnbr5(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin")
endTime <- Sys.time()
endTime-startTime
#5.79 mins --> 347 segs

#Versão 6

startTime <- Sys.time()
mgwnbr6(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin")
endTime <- Sys.time()
endTime-startTime
#6 mins --> 360 segs



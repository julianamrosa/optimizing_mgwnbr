library(sp)
library(mgwnbr)

data(georgia)

for (var in c("TotPop90", "PctRural", "PctEld", "PctFB", "PctPov", "PctBlack")){
  georgia[, var] <- as.data.frame(scale(georgia[, var]))
}

georgia$PctBach <- as.integer(georgia$PctBach)

startTime <- Sys.time()
mgwnbr1(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
       method="adaptive_bsq", bandwidth="aic", model="negbin", mgwr=FALSE)
endTime <- Sys.time()
endTime-startTime

startTime <- Sys.time()
mgwnbr1(data=georgia, formula=PctBach~TotPop90+PctRural+PctEld+PctFB+PctPov+PctBlack, lat="Y", long="X", globalmin=FALSE,
        method="adaptive_bsq", bandwidth="aic", model="negbin")
endTime <- Sys.time()
endTime-startTime

library(roxygen2) # In-Line Documentation for R 
library(devtools) # Tools to Make Developing R Packages Easier
library(testthat) # Unit Testing for R
library(usethis)  # Automate Package and Project Setup
library(rhub)

setwd("C:\\Juliana\\CRAN")

document("mgwnbr")
install("mgwnbr")

# Check for CRAN specific requirements using rhub and save it in the results 
# objects
verify_gwbr <- rhub::check_for_cran("mgwbr")
# Get the summary of your results
verify_gwbr$cran_summary()

# Generate your cran-comments.md, then you copy-paste the output from the function above
setwd("C:\\Juliana\\CRAN\\mgwnbr")
usethis::use_cran_comments()

# Submeter o pacote
devtools::release()

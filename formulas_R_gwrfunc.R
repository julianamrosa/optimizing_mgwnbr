formula_r <- function(formula, data = list()) {
  mf <- match.call(expand.dots = FALSE) #salva a chamada da função sem expandir o ...
  m <- match(c("formula", "data"), names(mf), 0) #salva os índices
  #de cada argumento (ordem em que aparecem na chamada da função)
  mf <- mf[c(1, m)] #arruma a chamada da função na sua ordem tradicional
  mf$drop.unused.levels <- TRUE #adiciona drop.unused.levels como parâmetro da função e passa
  #TRUE como argumento na chamada
  mf[[1]] <- as.name("model.frame") #monta a chamada do model.frame com os argumentos desta função
  #model.frame retorna um dataframe com as variáveis necessárias (vindas da fórmula)
  mf <- eval(mf) #executa model.frame() e salva o dataframe em mf
  mt <- attr(mf, "terms") #pega os termos do modelo gerado
  y <- model.extract(mf, "response") #pega a variável resposta
  dp.n <- length(y)
  x <- model.matrix(mt, mf) #monta a matriz de respostas no formato adequado para um modelo de regressão
  #já com intercepto
  print(x)
  print(y)
}

if (!require(spgwr)){
  install.packages("spgwr")
  library(spgwr)
}

data(georgia)
georgia <- as.data.frame(gSRDF)

formula_r(formula=PctBach~PctEld+TotPop90+PctPov, data = georgia)

formula_r <- function(formula, data = list()) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf)
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  dp.n <- length(y)
  x <- model.matrix(mt, mf)
  print(x)
  print(y)
}
formula_r(data=georgia_std, formula=PctBach~PctBlack+PctFB+TotPop90+PctEld)

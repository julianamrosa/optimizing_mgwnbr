formula_r <- function(formula, data = list()) {
  this.call <- match.call() #salva a chamada da função
  mf <- match.call(expand.dots = FALSE) #salva a chamada da função sem expandir o ...
  m <- match(c("formula", "data", "weights"), names(mf), 0) #salva os índices
  #de cada argumento (ordem em que aparecem na chamada da função)
  mf <- mf[c(1, m)] #arruma a chamada da função na sua ordem tradicional
  mf$drop.unused.levels <- TRUE #adiciona drop.unused.levels como parâmetro da função e passa
  #TRUE como argumento na chamada
  mf[[1]] <- as.name("model.frame") #monta a chamada do model.frame com os argumentos desta função
  #model.frame retorna um dataframe com as variáveis necessárias (vindas da fórmula)
  mf <- eval(mf, parent.frame()) #executa model.frame() e salva o dataframe em mf
  #envir=parent.frame é default da eval()
  mt <- attr(mf, "terms") #pega os termos do modelo gerado
  dp.n <- length(model.extract(mf, "response")) #pega length(y)
  y <- model.extract(mf, "response") #pega a variável resposta
  x <- model.matrix(mt, mf) #monta a matriz de respostas no formato adequado para um modelo de regressão
  #já com intercepto
  print(x)
  print(y)
}

formula_r(formula=PctBach~PctEld+TotPop90+PctPov, data = georgia)

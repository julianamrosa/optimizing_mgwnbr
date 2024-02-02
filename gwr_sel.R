gwr.sel <- function(formula, data = list(), coords, adapt=FALSE, 
                    gweight=gwr.Gauss, method="cv", verbose=TRUE, longlat=NULL,
                    RMSE=FALSE, weights, tol=.Machine$double.eps^0.25,
                    show.error.messages=FALSE) {
  if (!is.logical(adapt)) stop("adapt must be logical")
  if (is(data, "Spatial")) {
    if (!missing(coords))
      warning("data is Spatial* object, ignoring coords argument")
    coords <- coordinates(data)
    if (is.null(longlat) || !is.logical(longlat)) {
      if (!is.na(is.projected(data)) && !is.projected(data)) {
        longlat <- TRUE
      } else {
        longlat <- FALSE
      }
    }
    data <- as(data, "data.frame")
  }
  if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
  if (missing(coords))
    stop("Observation coordinates have to be given")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  dp.n <- length(model.extract(mf, "response"))
  #	mt <- terms(formula, data = data)
  #	mf <- lm(formula, data, method="model.frame", na.action=na.fail)
  #	dist2 <- (as.matrix(dist(coords)))^2
  weights <- as.vector(model.extract(mf, "weights"))
  # set up default weights
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (is.null(weights)) weights <- rep(as.numeric(1), dp.n)
  if (any(is.na(weights))) stop("NAs in weights")
  if (any(weights < 0)) stop("negative weights")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  #	if (NROW(x) != NROW(dist2))
  #		stop("Input data and coordinates have different dimensions")
  #### parei aqui (antes disso é só conferência de argumentos e organização incial do modelo)
  if (!adapt) { #FIXED
    bbox <- cbind(range(coords[,1]), range(coords[,2])) #cria uma matriz onde a primeira coluna tem o máximo e o mínimo da longitude e a segunda o mesmo para a latitude
    difmin <- spDistsN1(bbox, bbox[2,], longlat)[1] #calcula a distância euclidiana entre o ponto máximo e o mínimo
    #o spDistsN1 vai retornar um vetor cujo tamanho equivale ao número de linhas da matriz passada
    #(pois para cada ponto - linha - da matriz, é calculada a distância em relação ao ponto passado)
    #nesse caso, retorna duas distâncias, sendo a primeira entre o ponto mínimo e o máximo e a segunda entre o máximo e o máximo (0)
    #por isso pegamos o primeiro valor
    if (any(!is.finite(difmin)))
      difmin[which(!is.finite(difmin))] <- 0
    beta1 <- difmin/1000
    beta2 <- difmin
    if (method == "cv") { #CV
      opt <- optimize(gwr.cv.f, lower=beta1, upper=beta2, #buscar essa função gwr.cv.f no repositório
                      maximum=FALSE, y=y, x=x, coords=coords, 
                      gweight=gweight, verbose=verbose, 
                      longlat=longlat, RMSE=RMSE, weights=weights, 
                      show.error.messages=show.error.messages,
                      tol=tol) #acha o valor ótimo do cv
    }
    else { #AIC
      opt <- optimize(gwr.aic.f, lower=beta1, upper=beta2, 
                      maximum=FALSE, y=y, x=x, coords=coords, 
                      gweight=gweight, verbose=verbose, 
                      longlat=longlat, 
                      show.error.messages=show.error.messages,
                      tol=tol)
    }
    bdwt <- opt$minimum
    res <- bdwt #valor mínimo do cv/aic
  }
  else { #ADAPTIVE
    beta1 <- 0
    beta2 <- 1
    if (method == "cv") {
      opt <- optimize(gwr.cv.adapt.f, lower=beta1, 
                      upper=beta2, maximum=FALSE, y=y, x=x, 
                      coords=coords, gweight=gweight, 
                      verbose=verbose, longlat=longlat, RMSE=RMSE, 
                      weights=weights, 
                      show.error.messages=show.error.messages,
                      tol=tol) #a função cv é diferente e alguns argumentos também --> verificar no repo
    } else {
      opt <- optimize(gwr.aic.adapt.f, lower=beta1, 
                      upper=beta2, maximum=FALSE, y=y, x=x, 
                      coords=coords, gweight=gweight, 
                      verbose=verbose, longlat=longlat, 
                      show.error.messages=show.error.messages,
                      tol=tol)
    }
    q <- opt$minimum
    res <- q
  }
  if (isTRUE(all.equal(beta2, res, tolerance=.Machine$double.eps^(1/4))))
    warning("Bandwidth converged to upper bound:", beta2) #gera aviso se o valor mínimo de cv é igual ao beta2
  res
}

gwr.cv.f <- function(bandwidth, y, x, coords, gweight, verbose=TRUE, 
                     longlat=FALSE, RMSE=FALSE, weights, show.error.messages=TRUE) {
  n <- NROW(x)
  #    m <- NCOL(x)
  cv <- numeric(n) #salva um vetor numérico de tamanho n
  options(show.error.messages = show.error.messages) #
  for (i in 1:n) {
    xx <- x[i, ] #loop que percorre linha por linha de x
    dxs <- spDistsN1(coords, coords[i,], longlat=longlat) #calcula a distância de todos os pontos até o i-ésimo ponto
    #dxs é um vetor de tamanho n com as distâncias ao i-ésimo local
    if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bandwidth)
    #w.i <- exp((-0.5)*((dxs^2)/(bandwidth^2))) --> a linha de cima está fazendo isso
    #ou seja, ela aplica uma 'transformação' nas distâncias quadradas, o resultado ainda é um vetor de tamanho n
    #	w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bandwidth)
    w.i[i] <- 0 #o peso do i-ésimo ponto na calibração dele mesmo deve ser zero (pela formula acima, seria 1)
    w.i <- w.i * weights
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- try(lm.wfit(y = y, x = x, w = w.i)) #tenta executar a regressão linear ponderada
    if(!inherits(lm.i, "try-error")) { #se o lm não deu erro
      b <- coefficients(lm.i)
      cv[i] <- weights[i] * y[i] - (t(b) %*% (weights[i] * xx)) #calcula cv para a i-ésima linha e salva no vetor cv
    }
  }
  score <- sum(t(cv) %*% cv) #depois de calcular todos esses cvs, faz sua soma quadrática para pegar o escore de fato
  if (RMSE) score <- sqrt(score/n)
  #    score <- sqrt(sum(t(cv) %*% cv)/n)
  if (!show.error.messages) options(show.error.messages = TRUE)
  if (verbose) cat("Bandwidth:", bandwidth, "CV score:", score, "\n")
  score
}

long <- georgia[, 'X']
lat <- georgia[, 'Y']
COORD <<- matrix(c(long, lat), ncol=2)
distance <- dist(COORD, "euclidean")
distance <- as.matrix(distance)
sum(distance[, 1]==distance[1, ])

distance[, 1] #fazer isso usando spDistsN1()
spDistsN1(COORD, COORD[1,])
sum(spDistsN1(COORD, COORD[1,])==distance[, 1])

#lembrar de pegar a distância máxima

##implementar isso na função mgwnbr4
##testar a diferença de tempo!!

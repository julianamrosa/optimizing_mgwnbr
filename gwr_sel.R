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
    else { #AIC?
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
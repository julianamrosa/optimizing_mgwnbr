minha_func <- function(formula, data){
  model_desc <- parse.formula(formula)
  print(model_desc)
  y <- as.vector(data[[model_desc$response]])
  print(y)
}

library(mosaicCore)
minha_func(y~x1+x2, dados)
model_desc <- parse.formula(y~x1+x2)
model_desc
eval(model_desc$lhs)

library(DescTools)
form = ParseFormula(fo, data=dados)
form$lhs
form$rhs

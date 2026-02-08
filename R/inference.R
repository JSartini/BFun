#' Produce posterior credible interval for any matrix-structured parameter
#'
#' @details Estimates posterior credible intervals for each matrix entry
#' (point-wise)
#'
#' @param Var, list of matrix variate samples
#' @param Alpha, alpha level (type 1 error) for the intervals
#'
#' @return named list of credible interval bounds in matrix format
#'
Mat_Infer <- function(Var, Alpha){

  Lower_Quant = Alpha/2
  Upper_Quant = 1-Lower_Quant

  Var_Samp = Var %>%
    abind(along = 3)

  LB = Var_Samp %>%
    apply(c(1,2), quantile, probs = c(Lower_Quant))
  UB = Var_Samp %>%
    apply(c(1,2), quantile, probs = c(Upper_Quant))

  return(list(Lower = LB, Upper = UB))
}

#' Calculate posterior credible intervals for latent smooth functions
#'
#' @details Calculates posterior credible intervals using FPCA component samples
#'
#' @param Samples, aligned and formated posterior samples
#' @param B, orthonormal basis matrix
#' @param Data, data used to fit the model
#' @param Type, model type/structure
#' @param Alpha, alpha level (type 1 error) for the intervals
#'
#' @return named list of credible interval bounds in matrix format (# subjects
#' x # observations)
#'
Smooth_Infer <- function(Samples, B, Data, Type, Alpha){
  X = Data$X
  Lower_Quant = Alpha/2
  Upper_Quant = 1-Lower_Quant

  if(Type == "Fixed Effects"){
    Smooths = map(Samples$wFE, function(beta){
      Latent = X %*% t(beta)
      return(Latent %*% t(B))
    }) %>%
      abind(along = 3)
  }
  else if(Type == "Single-level"){
    Smooths = pmap(list(Samples$wFE, Samples$W1, Samples$S1),
                  function(beta, w, s){
                    Latent = X %*% t(beta) + (s %*% t(w))[Data$ID,]
                    return(Latent %*% t(B))
                  }) %>%
      abind(along = 3)
  }
  else{
    Smooths = pmap(list(Samples$wFE, Samples$W1, Samples$S1,
                       Samples$W2, Samples$S2),
                  function(beta, w1, s1, w2, s2){
                    Latent = X %*% t(beta) +
                      (s1 %*% t(w1))[Data$ID,] +
                      (s2 %*% t(w2))
                    return(Latent %*% t(B))
                  }) %>%
      abind(along = 3)
  }

  LB = Smooths %>%
    apply(c(1,2), quantile, probs = c(Lower_Quant))
  UB = Smooths %>%
    apply(c(1,2), quantile, probs = c(Upper_Quant))

  return(list(Lower = LB, Upper = UB))
}

#' Calculate posterior credible intervals for all components of interest
#'
#' @details Calculates credible intervals using quantiles of drawn samples
#'
#' @param Samples, aligned and formated posterior samples
#' @param B, orthonormal basis matrix
#' @param Data, data used to fit the model
#' @param Type, model type/structure
#' @param Alpha, alpha level (type 1 error) for the intervals
#'
#' @return named list of posterior credible interval quantities, mirrors the
#' output of extract and Posterior_Estimates
#'
Posterior_Infer <- function(Samples, B, Data, Type, Alpha = 0.05){

  output = list()

  # Fixed effects and parameters present in all models
  output$FE = Mat_Infer(Samples$FE, Alpha)

  # Random effects by structure
  if(Type %in% c("Single-level", "Multilevel")){
    output$S1 = Mat_Infer(Samples$S1, Alpha)
    output$EV1 = Mat_Infer(Samples$EV1, Alpha)
    output$EF1 = Mat_Infer(Samples$EF1, Alpha)

    if(Type == "Multilevel"){
      output$S2 = Mat_Infer(Samples$S2, Alpha)
      output$EV2 = Mat_Infer(Samples$EV2, Alpha)
      output$EF2 = Mat_Infer(Samples$EF2, Alpha)
    }
  }

  # Latent smooths
  output$YHat = Smooth_Infer(Samples, B, Data, Type, Alpha)

  return(output)
}

#' Produce posterior estimate of the eigenfunctions
#'
#' @details Estimates mean spline coefficients for each observed data function,
#' then leverages that the eigenvectors in the latent spline space correspond
#' to eigenfunctions in the data space.
#'
#' @param W, list of eigenfunction spline weight (Psi) samples
#' @param S, list of eigenfunction score samples
#' @param B, orthonormal basis matrix
#' @param K, number of principal components
#'
#' @return posterior estimate of the eigenfunctions in matrix form
#'
FPC_Est <- function(W, S, B, K){

  W_hat = pmap(list(W, S), function(w, s){
    return(s %*% t(w))
  }) %>%
    abind(along = 3) %>%
    apply(c(1,2), mean)
  FPC_hat = B %*% svd(W_hat, nv = K)$v

  return(FPC_hat)
}

#' Produce posterior estimate of the fixed effects
#'
#' @details Estimates mean spline coefficients for the fixed effects, then
#' applies the basis transformation
#'
#' @param W_FE, list of fixed effects spline weight (Psi) samples
#' @param B, orthonormal basis matrix
#'
#' @return posterior estimate of the fixed effects in matrix form
#'
FE_Est <- function(W_FE, B){

  W_hat = W_FE %>%
    abind(along = 3) %>%
    apply(c(1,2), mean)
  FE_hat = B %*% W_hat

  return(FE_hat)
}

#' Produce posterior estimate of any matrix-structured parameters
#'
#' @details Estimates posterior mean value for each matrix entry
#'
#' @param Var, list of matrix variate samples
#'
#' @return posterior estimate of the matrix parameter
#'
Mat_Est <- function(Var){

  Var_hat = abind(Var, along = 3) %>%
    apply(c(1,2), mean)

  return(Var_hat)
}

#' Calculate posterior estimates for latent smooth functions
#'
#' @details Calculates posterior mean estimates using FPCA component samples
#'
#' @param Samples, aligned and formated posterior samples
#' @param B, orthonormal basis matrix
#' @param Data, data used to fit the model
#' @param Type, model type/structure
#'
#' @return #participants x #observations matrix of posterior mean smooth estimates
#'
Smooth_Est <- function(Samples, B, Data, Type){
  X = Data$X

  if(Type == "Fixed Effects"){
    Latent = map(Samples$wFE, function(beta){
      return(X %*% t(beta))
    }) %>%
      abind(along = 3) %>%
      apply(c(1,2), mean)
  }
  else if(Type == "Single-level"){
    Latent = pmap(list(Samples$wFE, Samples$W1, Samples$S1),
                      function(beta, w, s){
      return(X %*% t(beta) + (s %*% t(w))[Data$ID,])
    }) %>%
      abind(along = 3) %>%
      apply(c(1,2), mean)
  }
  else{
    Latent = pmap(list(Samples$wFE, Samples$W1, Samples$S1,
                       Samples$W2, Samples$S2),
                  function(beta, w1, s1, w2, s2){
                    output = X %*% t(beta) +
                      (s1 %*% t(w1))[Data$ID,] +
                      (s2 %*% t(w2))
                    return(output)
                  }) %>%
      abind(along = 3) %>%
      apply(c(1,2), mean)
  }
  return(Latent %*% t(B))
}

#' Calculate posterior estimates for all components of interest
#'
#' @details Calculates posterior mean estimates from drawn samples
#'
#' @param Samples, aligned and formated posterior samples
#' @param B, orthonormal basis matrix
#' @param Data, data used to fit the model
#' @param Type, model type/structure
#'
#' @return named list of posterior estimates mirroring the output of extract
#'
Posterior_Estimates <- function(Samples, B, Data, Type){

  output = list()

  # Fixed effects and parameters present in all models
  output$FE = FE_Est(Samples$wFE, B)

  # Random effects by structure
  if(Type %in% c("Single-level", "Multilevel")){
    output$S1 = Mat_Est(Samples$S1)
    output$EV1 = Mat_Est(Samples$EV1)
    output$EF1 = FPC_Est(Samples$W1, Samples$S1, B, Data$K1)

    if(Type == "Multilevel"){
      output$S2 = Mat_Est(Samples$S2)
      output$EV2 = Mat_Est(Samples$EV2)
      output$EF2 = FPC_Est(Samples$W2, Samples$S2, B, Data$K2)
    }
  }

  # Latent smooths
  output$YHat = Smooth_Est(Samples, B, Data, Type)

  return(output)

}

#' Calculates spline bases and corresponding penalty matrices
#'
#' @details Able to generate a number of spline bases: B-splines,
#' orthogonalized B-splines, localized orthonormal splines, and vector
#' orthonormal splines.
#'
#' @param basis_name, name of the spline basis to generate objects for
#' @param basis_dim, dimension of the desired spline basis
#' @param args, argument points along which to
#'
#' @return List containing the "Basis" and penalty matrix "Pen"
#'
calc_splines <- function(basis_name, basis_dim, args){
  if(basis_name == "Splinets"){
    B_f = Splinets_basis(basis_dim)
    B = map(B_f$B, function(f){
      return(f(args))
    }) %>%
      abind(along = 2)

    D2 = Splinets_deriv2(basis_dim, bobj$cInt, bobj$cSlo, 10)
    P_mat = matrix(0, ncol = basis_dim, nrow = basis_dim)
    for(i in 1:basis_dim){
      for(j in 1:i){
        f2 <- function(x){ return(D2[[i]](x) * D2[[j]](x)) }
        P_mat[i,j] = stats::integrate(f2, lower = 0, upper = 10,
                                      subdivisions = 10000)$value
        if(i != j){
          P_mat[j,i] = P_mat[i,j]
        }
      }
    }
  }
  else if(basis_name == "B"){
    B_f = SplineBasis(c(rep(0, 3), seq(0, 1, length.out = basis_dim-2), rep(1, 3)))
    B = evaluate(B_f, args)
    P_mat = OuterProdSecondDerivative(B_f)
  }
  else if(basis_name == "OB"){
    B_f = OBasis(c(rep(0, 3), seq(0, 1, length.out = basis_dim-2), rep(1, 3)))
    B = evaluate(B_f, args)
    P_mat = OuterProdSecondDerivative(B_f)
  }
  else if(basis_name == "Vec"){
    B_f = OBasis(c(rep(0, 3), seq(0, 1, length.out = basis_dim-2), rep(1, 3)))
    B = evaluate(B_f, args)

    # Vector orthonormalize
    svd_objs = svd(B)
    A = svd_objs$v %*% diag(1/svd_objs$d) %*% t(svd_objs$v)
    B_ortho = B %*% A

    P_mat = t(A) %*% OuterProdSecondDerivative(B_f) %*% A
  }
  else{
    stop("Basis not supported")
  }
  return(list(Basis = B, Penalty = P_mat))
}

#' Numerical integration to derive penalty matrices
#'
#' @details Uses numerical integration techniques to calculate penalty
#' matrices for arbitrary basis functions
#'
#' @param d2, evaluated second derivatives of the basis
#' @param basis_dim, dimension of the desired spline basis
#' @param upper, upper extent of the evaluation integral
#'
#' @return Outer product of second derivative functions for penalty
#'
pen_numeric <- function(d2, basis_dim, upper = 1){
  P2 = matrix(0, ncol = basis_dim, nrow = basis_dim)
  for(i in 1:Q){
    for(j in 1:i){
      f2 <- function(x){ return(d2[[i]](x) * d2[[j]](x)) }
      P2[i,j] = stats::integrate(f2, lower = 0, upper = upper,
                                 subdivisions = 10000)$value
      if(i != j){
        P2[j,i] = P2[i,j]
      }
    }
  }
  return(P2)
}

#' Generate Splinet basis augmented with intercept and slope
#'
#' @details Generates orthonormalized B-spline basis which maintains locality
#' of the standard B-splines, assisting in numerical stability and efficiency
#'
#' @param basis_dim, dimension of the desired spline basis
#'
#' @return List containing the basis (B - list of functions) and
#' inner products between intercept/slope and splinets (cInt, cSlo)
#'
Splinets_basis <- function(basis_dim){
  B_f = splinet(knots = seq(0, 1, length.out = basis_dim+2), norm = T)
  basis = map(1:(basis_dim-2), function(q){
    func <- function(x){
      return(evspline(B_f$os, sID = q, x = x)[,-1])
    }
    return(func)
  })

  cs = map(1:(basis_dim-2), function(q){
    return(integral(function(x){return(x*basis[[q]](x))},
                    xmin = 0, xmax = 1))
  }) %>% unlist()
  slopeF = function(x){
    output = x - evspline(B_f$os, x = x)[,-1] %*% cs
    mag = 1/3 - sum(cs^2)
    return(output/sqrt(mag))
  }
  basis[[basis_dim-1]] = slopeF

  ci = map(1:(basis_dim-1), function(q){
    return(integral(basis[[q]], xmin = 0, xmax = 1))
  }) %>% unlist()
  intF = function(x){
    output = rep(1, length(x)) -
      evspline(B_f$os, x = x)[,-1] %*% ci[1:(basis_dim-2)] -
      slopeF(x)*ci[basis_dim-1]
    mag = 1 - sum(ci^2)
    return(output/sqrt(mag))
  }
  basis[[basis_dim]] = intF
  return(list(B = basis, cInt = ci, cSlo = cs))
}

#' Generate second derivatives of Splinet basis augmented with intercept
#' and slope
#'
#' @details Generates second derviatives of orthonormalized B-spline basis
#' which maintains locality of the standard B-splines
#'
#' @param basis_dim, dimension of the desired spline basis
#' @param cInt, inner products between splinets and intercept
#' @param cSlo, inner products between splinets and slope
#' @param Upper, upper bound of splines, can be increased for
#' numerical stability
#'
#' @return List of functions containing the second derivatives of basis
#'
Splinets_deriv2 <- function(basis_dim, cInt, cSlo, Upper = 1){

  B_f = splinet(knots = seq(0, 1, length.out = basis_dim+2), norm = T)
  second_deriv = deriva(deriva(B_f$os))

  derivs = map(1:(basis_dim-2), function(q){
    func <- function(x){
      return(evspline(second_deriv, sID = q, x = x)[,-1])
    }
    return(func)
  })

  slopeD2 = function(x){
    output = - evspline(second_deriv, x = x)[,-1] %*% cSlo
    mag = 1/3 - sum(cSlo^2)
    return(output/sqrt(mag))
  }
  derivs[[basis_dim-1]] = slopeD2

  intD2 = function(x){
    output = - evspline(second_deriv, x = x)[,-1] %*% cInt[1:(basis_dim-2)] -
      slopeD2(x)*cInt[basis_dim-1]
    mag = 1 - sum(cInt^2)
    return(output/sqrt(mag))
  }
  derivs[[basis_dim]] = intD2
  return(derivs)
}

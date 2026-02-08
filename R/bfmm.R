#' Fits a Bayesian Functional (Mixed) Model with scalar predictors
#'
#' @details Fits a flexible Bayesian functional mixed effects model using the
#' FAST approach to functional PCA for representing random effects
#'
#' @param form, a formula indicating the functional response (in matrix form)
#' and the predictor covariates (scalar)
#' @param data, dataframe containing the data to be fit. Functional data
#' should be in matrix format.
#' @param alpha, alpha level (type 1 error rate) for pposterior inferences
#' @param id, column name containing subject-level identifiers (if available)
#' @param K1, number of FPCs used at the first level, default chosen empirically
#' to account for 95% variability
#' @param visit, column name containing visit-level identifiers (if available)
#' @param K2, number of FPCs used at the second level, default chosen empirically
#' to account for 95% variability
#' @param spline_basis, type of spline basis used. Several are supported:
#' * **B**: Simple degree 3 B-splines, the default
#' * **Splinets**: A localized, orthonormal basis
#' * **OB**: Orthonormalized B-splines
#' * **Vec**: Vector orthonormal basis created through SVD
#' @param spline_dim, dimension of the spline basis, defaults to 25
#' @param args, argument points at which the function is evaluated, defaults to
#' equally-spaced points between 0 and 1
#' @param out_args, output argument points at which functions are evaluated, defaults
#' to args
#' @param method, Bayesian computation algorithm used, options are:
#' * **MCMC**: Markov chain Monte Carlo simulation
#' * **MF**: Mean field variational approximation (EXPERIMENTAL)
#' @param n_iter, number of total sampling iterations for MCMC, defaults to 2000
#' @param n_burnin, number of burn-in samples for MCMC, defaults to n_iter/2
#' @param n_chains, number of chains for MCMC, defaults to 1
#' @param n_cores, number of machine cores to use when running MCMC, defaults to n_chains
#' @param treedepth, maximum tree depth for HMC sampler
#' @param adapt_delta, adaptation target acceptance for HMC sampler
#'
#' @return BFM object containing posterior samples suitable for analysis
#'
#' @references Sartini, J., Zhou, X., Selvin, L., Zeger, S., & Crainiceanu, C. (2025).
#'   Fast Bayesian Functional Principal Components Analysis.
#'   Journal of Computational and Graphical Statistics.
#'
#' @export
#'
bfmm <- function(form, data, alpha = 0.05, id = NULL, visit = NULL,
                 K1 = NULL, K2 = NULL, spline_basis = "OB", spline_dim = 25,
                 args = NULL, out_args = NULL, method = "MCMC", n_iter = 2000,
                 n_burnin = floor(n_iter/2), n_chains = 1, n_cores = 1,
                 treedepth = NULL, adapt_delta = NULL){

  inputs_const = model_input(form, data, id, visit, K1, K2, spline_basis,
                                 spline_dim, args, out_args)
  B_out = inputs_const$new_B
  Type = inputs_const$type
  Mod_File = inputs_const$file
  model = instantiate::stan_package_model(name = Mod_File, package = "BFun")
  VarNames = inputs_const$x_names
  message("Inputs formatted, running model")

  # Run the model
  if(method == "MCMC"){
    if(is.null(n_cores)){
      n_cores = n_chains
    }

    model_fit = model$sample(
      data = inputs_const$input_data,
      refresh = floor((n_iter + n_burnin)/10),
      chains = n_chains,
      parallel_chains = n_chains,
      iter_warmup = n_burnin,
      iter_sampling = n_iter,
      max_treedepth = treedepth,
      adapt_delta = adapt_delta,
      show_exceptions = F
    )
  }
  else if(method == "MF"){
    model_fit = model$variational(
      data = inputs_const$input_data,
      draws = n_iter,
      algorithm = "meanfield",
      show_exceptions = F
    )
  }
  else{
    stop("Method not supported")
  }
  message("Model run, formatting outputs")

  # Extract outputs and align
  Samples = cmd_extract_permuted(model_fit)
  Formatted = extract_permuted(Samples, inputs_const$input_data, B_out, Type)
  Aligned = Align_Posterior(Formatted, inputs_const$input_data, B_out, Type)

  # Posterior Estimates and inference
  Estimates = Posterior_Estimates(Aligned, B_out, inputs_const$input_data, Type)
  Inferences = Posterior_Infer(Aligned, B_out, inputs_const$input_data, Type,
                               Alpha = alpha)
  Summary = model_output(Estimates, Inferences, Type, VarNames)
  message("Outputs produced, running convergence diagnostics")

  # Model diagnostics
  Var_Exp = Var_Exp(Estimates, inputs_const$input_data, Type)
  Diagnostics = list(Variance_Explained = Var_Exp)
  if(method == "MCMC"){
    chain_samples = cmd_extract_chains(model_fit)
    format_chains = extract_chain(chain_samples, inputs_const$input_data, B_out, Type)
    converge = RHat_FullMod(format_chains, Estimates, inputs_const$input_data,
                            VarNames, Type)
    Diagnostics$RHat = converge
  }
  Summary$Diagnostics = Diagnostics
  message("Done with diagnostics")

  return(Summary)
}


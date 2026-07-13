# Wrapper that initialises theta and runs the penalised (LASSO) MSAR EM fit with
# bounded retries on failure. Returns the fitted model, or NULL if all retries fail.
init_and_fit_msar_lasso <-
  function(data,
           M,
           order,
           MaxIter = 200,
           retry = 4,
           eps = 1e-5,
           verbose = FALSE,
           lambda_fuse_A = 0,
           lambda_fuse_sigma = 0) {
    result <- list(fit = NULL, error = NULL)
    
    
    for (try in 1:(retry + 1)) {
      if (try > 1) {
        message("retry: ", try, " of ", retry)
      }
      

      # try (initialization is inside the tryCatch so that a degenerate
      # init -- e.g. a kmeans split that starves a regime and yields a
      # singular/NA covariance in NHMSAR:::Mstep.classif (det(Cxx) == NA) --
      # triggers a retry with a fresh initialization instead of aborting the
      # whole run. model_init is recomputed every loop iteration, so each
      # retry still gets a new random initialization.)
      result <- tryCatch({

        # Initialize model for NHMSAR
        model_init <-
          init_theta_msar(
            data = data,
            M = M,
            order = order,
            label = "HH",
            verbose = verbose
          )

        # Fit MSAR model with NHMSAR
        model_fit <-
          fit_msar(
            data = data,
            theta = model_init,
            penalty = "LASSO",
            MaxIter = MaxIter,
            eps = eps,
            verbose = verbose,
            lambda_fuse_A = lambda_fuse_A,
            lambda_fuse_sigma = lambda_fuse_sigma
          )

        list(fit = model_fit, error = NULL)

      }, error = function(e) {
        message("init/fit_msar: can't fit: \n\tError: ", e$message)
        # Prefix tags the failure kind so the diagnostic failure_log /
        # summarize_failures() can distinguish genuine errors from warnings
        # that were promoted to fatal by this handler.
        list(fit = NULL, error = paste0("error: ", e$message))
        # brower()
      }, warning = function(w) {
        message("init/fit_msar: can't fit: \n\tWarning: ", w$message)
        # NOTE: any warning is currently promoted to a fatal failure (fit = NULL).
        # The "warning:" tag lets us measure how many discards are warnings only,
        # i.e. how much of the missingness this strictness manufactures.
        list(fit = NULL, error = paste0("warning: ", w$message))
        # brower()
      })
      
      fit <- result[["fit"]]
      if (!is.null(fit)) {
        # we have a fit! -> break
        break
      }
      
    } # retry
    
    return(result)
  }
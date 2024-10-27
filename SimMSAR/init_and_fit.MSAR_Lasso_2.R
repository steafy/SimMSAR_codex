init_and_fit.MSAR_Lasso_2 <-
  function(data,
           M,
           order,
           MaxIter = 200,
           retry = 4,
           eps = 1e-5,
           verbose = FALSE) {
    result <- list(fit = NULL, error = NULL)
    
    # browser()
    
    for (try in 1:(retry + 1)) {
      if (try > 1) {
        message("retry: ", try, " of ", retry)
      }
      
      
      # Initialize model for NHMSAR
      model_init <-
        init.theta.MSAR_revised_2(
          data = data,
          M = M,
          order = order,
          label = "HH",
          verbose = verbose
        )
      
      # try
      result <- tryCatch({
        # browser()
        
        # Fit MSAR model with NHMSAR
        model_fit <-
          fit.MSAR_revised_2(
            data = data,
            theta = model_init,
            penalty = "LASSO",
            MaxIter = MaxIter,
            eps = eps,
            verbose = verbose
          )
        
        list(fit = model_fit, error = NULL)
        
      }, error = function(e) {
        message("fit.MSAR_revised_2: can't fit: \n\tError: ", e$message)
        list(fit = NULL, error = e$message)
        # brower()
      }, warning = function(w) {
        message("fit.MSAR_revised_2: can't fit: \n\tWarning: ", w$message)
        list(fit = NULL, error = w$message)
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
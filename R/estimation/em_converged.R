# -----------------------------------------------------------------------------
# Adapted from the NHMSAR package (Valerie Monbet), obtained from the CRAN
# archive. Original license: GPL. See the repository LICENSE (GPL-3) and the
# README section "License & attribution".
# Adaptation: renamed to em_converged(); logic unchanged apart from formatting.
# -----------------------------------------------------------------------------

# EM convergence test: TRUE once the relative log-likelihood increase falls below
# `threshold`; also flags a (numerically) decreasing step. Adapted from the
# NHMSAR package.
em_converged <-
function(loglik, previous_loglik, threshold = 1e-4, verbose = FALSE) {

  converged = 0;
  decrease = 0;
  if(!(previous_loglik==-Inf)){
  if (loglik - previous_loglik < -1e-2) # allow for a little imprecision 
    {
      if( verbose ) { 
        print(paste("******likelihood decreased from ",previous_loglik," to ", loglik,sep=""),quote = FALSE)
      }
      decrease = 1;
    }

  delta_loglik = abs(loglik - previous_loglik);
  avg_loglik = (abs(loglik) + abs(previous_loglik) + threshold)/2;
  bb = ((delta_loglik/avg_loglik) < threshold)
  if (bb) {converged = 1}
  }

  res <- NULL
  res$converged <- converged
  res$decrease <- decrease
  return(res)

}

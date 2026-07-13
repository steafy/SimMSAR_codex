#' M-Step with Cross-Validated LASSO Penalty for MSAR Models
#'
#' Maximization step of the EM algorithm with a *genuinely penalized*,
#' cross-validated LASSO for sparse network estimation. Updates autoregressive
#' parameters, intercepts, and covariance matrices using regime-weighted
#' observations.
#'
#' @param data 3D array of time series data (time × samples × variables).
#' @param theta Current parameter object (thetaMSAR).
#' @param FB Forward-backward output from E-step containing regime probabilities.
#' @param verbose Logical. Print progress messages. Default: FALSE.
#'
#' @return List of updated parameters:
#'   \describe{
#'     \item{A}{List of lag-1 coefficient matrices (one per regime)}
#'     \item{A0}{Matrix of intercepts}
#'     \item{sigma}{List of covariance matrices}
#'     \item{prior}{Updated regime prior probabilities}
#'     \item{transmat}{Updated transition matrix}
#'   }
#'
#' @details
#' For each regime m and each response node id, a weighted LASSO is fit with
#' \code{glmnet::cv.glmnet()} on the RAW (unweighted) pseudo-observations, with
#' the E-step regime-membership probabilities \code{gamma} passed via the native
#' \code{weights=} argument (correct WLS: minimises \eqn{\sum_t w_t (y_t - x_t
#' b)^2}, unlike the old lars path which multiplied x and y by \code{w} and thus
#' weighted by \code{w^2}). The penalty parameter \eqn{\lambda} is chosen by
#' cross-validation, so the estimate is sparse on its own -- \emph{without}
#' relying on the downstream \code{min_edg_val} threshold.
#'
#' Two design knobs (see \code{docs/MSTEP_LASSO_CV_PENALIZATION.md}), read from
#' options so they can be A/B tested without changing call sites:
#' \itemize{
#'   \item \code{getOption("simmsar_lasso_lambda", "1se")}: \code{"1se"} (sparser)
#'     or \code{"min"} (lower CV error) rule for the chosen lambda.
#'   \item \code{getOption("simmsar_lasso_refit", TRUE)}: if \code{TRUE},
#'     coefficients on the selected support are re-estimated by unpenalized
#'     weighted OLS ("relaxed LASSO", removes shrinkage bias in the edge weights);
#'     if \code{FALSE}, the shrunk glmnet coefficients are used directly.
#' }
#'
#' Intercepts (\code{A0}) and covariances (\code{sigma}) are computed from the
#' correctly \code{w}-weighted sufficient statistics, exactly as before.
#'
#' Typically used only in the first EM iteration; subsequent iterations inherit
#' the selected support via \code{\link{mstep_hh_reduct_msar}}.
#'
#' @note Dependencies loaded centrally via R/dependencies.R (glmnet required).
#'
#' @seealso
#' \code{\link{fit_msar}} which calls this function
#' \code{\link{mstep_hh_reduct_msar}} for the reduced M-step (iterations 2+)
#'
#' @keywords internal
#' @export
# Dependencies are loaded centrally via R/dependencies.R
# Required packages: glmnet

mstep_hh_lasso_msar <-
function(data,theta,FB,verbose = FALSE)  {
  if (!exists("repmat", mode = "function")) {
    repmat <- function(x, m, n) {
      matrix(rep(x, m * n), nrow = m, ncol = n, byrow = TRUE)
    }
  }

  # ---- design knobs (A/B-testable via options; see docs) --------------------
  lambda_rule <- match.arg(getOption("simmsar_lasso_lambda", "1se"),
                           c("1se", "min"))
  do_refit    <- isTRUE(getOption("simmsar_lasso_refit", TRUE))
  do_adaptive <- isTRUE(getOption("simmsar_lasso_adaptive", FALSE))
  lambda_s    <- if (lambda_rule == "1se") "lambda.1se" else "lambda.min"
  # cv.glmnet cost knobs (see docs/MSTEP_LASSO_CV_PENALIZATION.md §5.2). With the
  # default relaxed refit (do_refit) cv.glmnet is only used to pick the SUPPORT
  # (coefficients are re-estimated exactly by WLS), so a cheaper CV (5 folds) and
  # a coarser lambda grid (50) barely move the selected support while cutting
  # per-call cost -- the dominant cost under per-iteration re-selection. Defaults
  # are the validated fast-and-lossless setting (5 / 50); 3 folds or <=30 lambda
  # start to lose recovery.
  cv_nfolds   <- getOption("simmsar_lasso_nfolds", 5)
  cv_nlambda  <- getOption("simmsar_lasso_nlambda", 50)
  # Fixed (deterministic) CV folds -- ON by default. Under per-iteration
  # re-selection, RANDOM folds add CV-sampling noise that makes the selected
  # support flicker between EM iterations, which slows convergence badly (e.g.
  # 5-fold random took 28 EM iters vs 7 with fixed folds) and eats any per-call
  # saving. A deterministic interleaved partition (identical every iteration,
  # no RNG touched) removes that noise source, so the support only changes when
  # the E-step genuinely changes -- turning the cheaper CV into a ~2.5x net
  # speedup with unchanged recovery.
  cv_fixfolds <- isTRUE(getOption("simmsar_lasso_fixedfolds", TRUE))

  T=dim(data)[1]
  N.samples = dim(as.array(data))[2]
  d = dim(as.array(data))[3]
  if(is.null(d)|is.na(d)) {d = 1}
  M <- attributes(theta)$NbRegimes
  p <- attributes(theta)$order
  order <- max(p,1)
  data = array(data,c(T,N.samples,d))
  data2 = array(0,c(order*d,T-order+1,N.samples))
  cpt=1
  for (o in order:1) {
    for (kd in 1:d) {
      data2[cpt,,] = data[o:(T-order+o),,kd]
      cpt =cpt+1
    }
  }
  T = length(o:(dim(data)[1]-order+o))
  exp_num_trans = 0
  exp_num_visit = 0
  exp_num_visits1 = 0
  postmix = 0
  m = matrix(0,d*order,M) ; m_1 = m ; c = matrix(0,M,1) ; s=c ;
  op = array(0,c(d*order,d*order,M) ); op_1 = op ; op_2 = op_1 ;

  # RAW (unweighted) design + weight vector per regime, for cv.glmnet's native
  # weights= (correct WLS). Predictors: d*order lagged values; responses: d.
  rawx = array(0,c(N.samples*(T-1),d*order,M))
  rawy = array(0,c(N.samples*(T-1),d,M))
  rw   = matrix(0,N.samples*(T-1),M)

  for (ex in 1:N.samples) {
    obs = array(data2[,,ex],c(d*order,T))
    xit  = array(FB$probSS[,,,ex],c(M,M,T-2))
    gamma = matrix(FB$probS[ex,1:(T-1),],T-1,M)
    exp_num_trans = exp_num_trans+apply(xit,c(1,2),sum)
    exp_num_visits1 = exp_num_visits1+gamma[1,]
    postmix = postmix+apply(gamma,2,sum)
    obs = t(obs)
    idx = (ex-1)*(T-1)+(1:(T-1))

    for (j in 1:M) {
      w = matrix(gamma[,j], 1,T-1)
      wobs = obs[1:(T-1),] * repmat(t(w),1,d*order)
      wobs_1= obs[2:T,] * repmat(t(w),1,d*order)

      # RAW (unweighted) pseudo-observations + weights for glmnet
      x_raw = obs[1:(T-1),]
      y_raw = obs[2:T,, drop = FALSE]
      wj    = gamma[,j]
      if (is.na(sum(obs))) {
        # zero out NA rows and give them zero weight (they must not enter the fit)
        na_x = !is.finite(rowSums(x_raw)); na_y = !is.finite(rowSums(y_raw))
        bad  = na_x | na_y
        x_raw[!is.finite(x_raw)] = 0
        y_raw[!is.finite(y_raw)] = 0
        wj[bad] = 0
        wobs[is.na(wobs)] = 0
        wobs_1[is.na(wobs_1)] = 0
        obs[is.na(obs)] = 0
      }
      rawx[idx,,j] = x_raw
      rawy[idx,,j] = y_raw[,1:d, drop = FALSE]
      rw[idx,j]    = wj

      m[,j] = m[,j] + apply(wobs,2,sum)
      m_1[,j] = m_1[,j] + apply(wobs_1, 2,sum)
      op[,,j] = op[,,j] + t(wobs) %*% obs[1:(T-1),]
      op_1[,,j] = op_1[,,j] + t(wobs) %*% obs[2:T,]
      op_2[,,j] = op_2[,,j] + t(wobs_1) %*% obs[2:T,]
    }
  }

  # Markov chain
  prior = NHMSAR:::normalise(exp_num_visits1)
  prior = matrix(prior,M,1)
  transmat = NHMSAR:::mk_stochastic(exp_num_trans)

  if ( any(is.nan(postmix)) || any(postmix < 1e-6, na.rm = TRUE) ) {
    stop("error : smoothing probabilities are to small, in one regime at least. You should revise initialisation.")
  }

  moy <- array(0,c(M,d))
  sigma <- list()
  A2 <-list()

  for (j in 1:M) {
    A2[[j]] = list()
    A2.lasso = matrix(0,d,d)
    Cxx = postmix[j]*op[,,j] - m[,j]%*%t(m[,j])
    Cxy = postmix[j]*op_1[,,j] - m[,j]%*%t(m_1[,j])

    Xj = rawx[,,j]
    wj = rw[,j]

    for (id in 1:d){
      yj = rawy[,id,j]

      # cross-validated LASSO with correct WLS weights. glmnet needs >=2
      # predictors and some variation; guard degenerate cases.
      support <- integer(0)
      A_shrunk <- rep(0, d*order)
      ok <- FALSE
      if (ncol(Xj) >= 2 && sum(wj > 0) > 2 && stats::sd(yj[wj > 0]) > 0) {
        # Adaptive LASSO (Zou 2006): per-predictor penalty.factor = 1/|b_init|,
        # with b_init the weighted-OLS solution (already available as
        # solve(Cxx) %*% Cxy). Large true coefficients get a small penalty so
        # they are NOT dropped/over-shrunk, which is what plain 1se does.
        pf <- rep(1, ncol(Xj))
        if (do_adaptive) {
          binit <- tryCatch(as.numeric(solve(Cxx, Cxy[,id])),
                            error = function(e) rep(1, ncol(Xj)))
          pf <- 1 / (abs(binit) + 1e-4)
        }
        # deterministic interleaved folds (1,2,..,K,1,2,..) -> identical every
        # EM iteration, no RNG touched; else glmnet's default random folds.
        foldid <- if (cv_fixfolds) ((seq_len(nrow(Xj)) - 1L) %% cv_nfolds) + 1L else NULL
        cvfit <- tryCatch(
          if (is.null(foldid)) {
            glmnet::cv.glmnet(Xj, yj, weights = wj, alpha = 1,
                              intercept = TRUE, standardize = TRUE,
                              penalty.factor = pf,
                              nfolds = cv_nfolds, nlambda = cv_nlambda)
          } else {
            glmnet::cv.glmnet(Xj, yj, weights = wj, alpha = 1,
                              intercept = TRUE, standardize = TRUE,
                              penalty.factor = pf,
                              foldid = foldid, nlambda = cv_nlambda)
          },
          error = function(e) NULL)
        if (!is.null(cvfit)) {
          # coef.cv.glmnet interprets the "lambda.1se"/"lambda.min" string and
          # returns coefficients at that CV-chosen lambda (verified identical to
          # glmnet::coef.glmnet on the underlying fit). First row is the intercept.
          cf <- as.numeric(stats::coef(cvfit, s = lambda_s))[-1]
          A_shrunk <- cf
          support <- which(abs(cf) > 1e-12)
          ok <- TRUE
        }
      }

      if (!ok) {
        # Degenerate fallback: unpenalized weighted OLS on the full support
        # (matches the old behaviour only in this rare edge case).
        support <- seq_len(d*order)
      }

      if (length(support) > 0) {
        if (do_refit || !ok) {
          # relaxed LASSO: unpenalized weighted-OLS refit on the selected support
          Cxx_w = Cxx[support,support]
          Cxy_w = Cxy[support,id]
          A2.lasso[id,support] = t(Cxy_w) %*% solve(Cxx_w)
        } else {
          # use the shrunk glmnet coefficients directly
          A2.lasso[id,support] = A_shrunk[support]
        }
      }
    }

    # Intercept and covariance from the correctly w-weighted sufficient stats
    tmp = (m_1[,j]-(A2.lasso)%*%m[,j])/postmix[j]
    op_1[,,j] = t(as.matrix(op_1[,,j]))
    tmp2 = (op_2[, , j] + (A2.lasso) %*% op[, , j] %*% t(A2.lasso) -
           ((A2.lasso) %*% t(op_1[, , j]) + t((A2.lasso) %*% t(op_1[, , j]))))/postmix[j] -
           tmp %*% t(tmp)

    if (verbose) {
      print(c("........Regime ",j,", nnz ",sum(A2.lasso!=0),"/",d*d),quote = FALSE)
    }

    A2[[j]][[1]] = A2.lasso
    # Optional in-EM stabilization of the residual covariance (default OFF; see
    # docs/SIGMA_KAPPA_DEGENERACY_DIAGNOSIS.md). Keeps Sigma dense; no-op unless
    # simmsar_sigma_stab is set. Applied here so the NEXT E-step's likelihood and
    # the stored est_Sigma both see a well-conditioned matrix.
    stab <- if (exists("stabilize_sigma", mode = "function")) stabilize_sigma else identity
    sigma[[j]]=stab(tmp2[1:d,1:d])
    moy[j,1:d] = tmp[1:d]
  }

  if (p>0) {
    list(A=A2,A0=moy,sigma=sigma,prior=prior,transmat=transmat)
  }
}

fit.MSAR <- function (data, theta, MaxIter = 100, eps = 1e-05, verbose = FALSE, 
  covar.emis = NULL, covar.trans = NULL, method = NULL, constraints = FALSE, 
  reduct = FALSE, K = NULL, d.y = NULL, ARfix = FALSE, penalty = FALSE, 
  sigma.diag = FALSE, sigma.equal = FALSE, lambda1 = 0.1, 
  lambda2 = 0.1, a = 3.7, ...) 
{
  cl <- match.call()
  now <- Sys.time()
  if (missing(theta)) {
    stop("can not fit a MSAR model without initial value for theta")
  }
  att <- attributes(data)
  att.theta <- attributes(theta)
  data <- as.array(data)
  T <- dim(data)[1]
  if (is.null(T) || is.na(T)) {
    T <- length(data)
  }
  N.samples <- dim(data)[2]
  if (is.null(N.samples) || is.na(N.samples)) {
    N.samples <- 1
  }
  d <- att.theta$NbComp
  if (is.null(d) | is.na(d)) {
    d <- 1
  }
  data <- array(data, c(T, N.samples, d))
  order <- att.theta$order
  label <- att.theta$label
  M <- att.theta$NbRegimes
  if ((missing(covar.trans) && substr(label, 1, 1) == "N") || 
    (missing(covar.emis) && substr(label, 2, 2) == "N")) {
    print("Can not fit a non homogeneous MSAR without covariable")
  }
  if (length(covar.trans) == 1) {
    Lag = covar.trans + 1
    covar.trans = array(data[(1):(T - Lag + 1), , ], c(T - 
      Lag + 1, N.samples, d))
    data = array(data[Lag:T, , ], c(T - Lag + 1, N.samples, 
      d))
  }
  if (!is.null(covar.emis)) {
    if (is.null(dim(covar.emis))) {
      ncov.emis = 1
    }
    else if (is.na(dim(covar.emis)[3])) {
      ncov.emis = 1
    }
    else if (!is.na(dim(covar.emis)[3])) {
      ncov.emis = dim(covar.emis)[3]
    }
    covar.emis = array(covar.emis, c(T, N.samples, ncov.emis))
  }
  else {
    ncov.emis = 0
  }
  if (!is.null(covar.trans)) {
    ncov.trans = dim(covar.trans)[3]
  }
  else {
    ncov.trans = 0
  }
  BIC = NULL
  Npar = NULL
  cnt <- 0
  FB <- Estep.MSAR(data, theta, covar.trans = covar.trans, 
    covar.emis = covar.emis)
  loglik = FB$loglik
  previous_loglik <- FB$loglik - 1000
  ll_history = NULL
  converged = EM_converged(0, 2 * eps, eps)
  par = NULL
  while (converged[1] == 0 && cnt < MaxIter) {
    cnt <- cnt + 1
    if (verbose) {
      print(c("iteration ", cnt, "  loglik = ", loglik), 
        quote = FALSE)
    }
    if (label == "HH") {
      if (constraints == FALSE & penalty == FALSE & reduct == 
        FALSE) {
        par = Mstep.hh.MSAR(data, theta, FB, sigma.diag = sigma.diag, 
          sigma.equal = sigma.equal)
      }
      else if (constraints) {
        par = Mstep.hh.MSAR.with.constraints(data, theta, 
          FB, K = K, d.y = d.y)
        attributes(theta)$n_par = M + M * (M - 1) + 
          2 * M * d
      }
      else if (penalty == "ridge") {
        par = Mstep.hh.SCAD.cw.MSAR(data, theta, FB, 
          penalty = "ridge", lambda1 = 0, lambda2 = lambda2, 
          par = par)
      }
      else if (penalty == "LASSO") {
        if (cnt > 1) {
          par = Mstep.hh.reduct.MSAR(data, theta, FB, 
            sigma.diag = sigma.diag)
        }
        else {
          par = Mstep.hh.lasso.MSAR(data, theta, FB)
        }
      }
      else if (penalty == "SCAD") {
        par = Mstep.hh.SCAD.MSAR(data, theta, FB, penalty = "SCAD", 
          lambda1 = lambda1, lambda2 = lambda2, par = par)
      }
      else if (reduct) {
        par = Mstep.hh.reduct.MSAR(data, theta, FB, 
          sigma.diag = sigma.diag)
      }
      if (order > 0) {
        theta = list(par$A, par$A0, par$sigma, par$prior, 
          par$transmat)
      }
      else {
        theta = list(par$A0, par$sigma, par$prior, par$transmat)
      }
    }
    else if (label == "HN") {
      par = Mstep.hn.MSAR(data, theta, FB, covar = covar.emis, 
        verbose = verbose)
      if (order > 0) {
        theta = list(par$A, par$A0, par$sigma, par$prior, 
          par$transmat, par$par_emis)
      }
      else {
        theta = list(par$A0, par$sigma, par$prior, par$transmat, 
          par$par_emis)
      }
    }
    else if (label == "NH") {
      par = Mstep.nh.MSAR(data, theta, FB, covar = covar.trans, 
        method = method, ARfix = ARfix, reduct = reduct, 
        sigma.diag = sigma.diag, sigma.equal = sigma.equal, 
        penalty = penalty, lambda1 = lambda1, lambda2 = lambda2, 
        par = par)
      if (order > 0) {
        theta = list(par$A, par$A0, par$sigma, par$prior, 
          par$transmat, par$par.trans)
      }
      else {
        theta = list(par$A0, par$sigma, par$prior, par$transmat, 
          par$par.trans)
      }
    }
    else if (label == "NN") {
      par = Mstep.nn.MSAR(data, theta, FB, covar.emis = covar.emis, 
        covar.trans = covar.trans, method = method)
      if (order > 0) {
        theta = list(par$A, par$A0, par$sigma, par$prior, 
          par$transmat, par$par.trans, par$par.emis)
      }
      else {
        theta = list(par$A0, par$sigma, par$prior, par$transmat, 
          par$par.trans, par$par.emis)
      }
    }
    ll_history[cnt] = loglik
    converged = EM_converged(loglik, previous_loglik, eps)
    previous_loglik = loglik
    attributes(theta) = att.theta
    theta = as.thetaMSAR(theta, label = label, ncov.emis = ncov.emis, 
      ncov.trans = ncov.trans)
    FB = Estep.MSAR(data, theta, covar.emis = covar.emis, 
      covar.trans = covar.trans)
    loglik = FB$loglik
  }
  if (M > 1) {
    tr = NULL
    for (m in 1:M) {
      if (d == 1) {
        tr[m] = theta$sigma[m]
      }
      else {
        tr[m] = sum(diag(theta$sigma[[m]]))
      }
    }
    i.tr = order(tr)
    theta$A0 = theta$A0[i.tr, ]
    sigma.tmp = NULL
    A.tmp = theta$A
    for (m in 1:M) {
      sigma.tmp[[m]] = theta$sigma[[i.tr[m]]]
      if (order > 0) {
        if (d > 1) {
          for (o in 1:order) {
            A.tmp[[m]][[o]] = theta$A[[i.tr[m]]][[o]]
          }
        }
        else {
          A.tmp[m, ] = theta$A[i.tr[m], ]
        }
      }
    }
    theta$A = A.tmp
    theta$sigma = sigma.tmp
    theta$prior = theta$prior[i.tr, ]
    temp = theta$transmat[i.tr, i.tr]
    attributes(temp)$dimnames[[2]] <- attributes(theta$transmat)$dimnames[[2]]
    attributes(temp)$dimnames[[2]] <- attributes(theta$transmat)$dimnames[[2]]
    theta$transmat = temp
    if (substr(label, 2, 2) == "N") {
      tmp = theta$par.emis
      for (j in 1:M) {
        theta$par.emis[[j]] = tmp[[i.tr[j]]]
      }
    }
    if (substr(label, 1, 1) == "N") {
      theta$par.trans = theta$par.trans[i.tr, ]
    }
    theta = as.thetaMSAR(theta, label = label, ncov.emis = ncov.emis, 
      ncov.trans = ncov.trans)
    FB$probS = FB$probS[, , i.tr]
  }
  if (penalty != "SCAD") {
    lambda1 = rep(0, M)
  }
  npar = M * d + M * (M - 1)
  if (substr(label, 1, 1) == "N") {
    npar = npar + M * length(theta$par.trans[1, ])
  }
  if (substr(label, 2, 2) == "N") {
    npar = npar + M * length(theta$par.emis[[1]])
  }
  for (m in 1:M) {
    npar = npar + sum(abs(theta$A[[m]][[1]]) > 0)
    if (penalty != "SCAD" | max(abs(lambda1)) == 0) {
      npar = npar + sum(abs(theta$sigma[[m]][upper.tri(theta$sigma[[m]], 
        diag = TRUE)]) > 0)
    }
    else {
      npar = npar + sum(abs(par$sigma.inv[[m]][upper.tri(par$sigma.inv[[m]], 
        diag = TRUE)]) > 1e-05)
    }
  }
  if (sigma.diag) {
    npar = (M - 1) + M * (M - 1) + M * d + M * d
  }
  if (sigma.equal) {
    npar = (M - 1) + M * (M - 1) + M * d + d * (d + 1)/2
  }
  if (sigma.diag & sigma.equal) {
    npar = (M - 1) + M * (M - 1) + M * d + d
  }
  attributes(theta)$n_par = npar
  BIC = -2 * ll_history[cnt] + attributes(theta)$n_par * log(length(c(data)))
  ll.pen = NULL
  if (penalty == "SCAD") {
    a = 3.7
    if (length(lambda1) == 1) {
      lambda1 = matrix(lambda1, 1, M)
    }
    if (length(lambda2) == 1) {
      lambda2 = matrix(lambda2, 1, M)
    }
    pen = 0
    for (m in 1:M) {
      w = matrix(0, d, d)
      if (lambda1[m] > 0) {
        wi = solve(theta$sigma[[m]])
        abs.S = abs(theta$sigma[[m]])
        w = lambda1[m] * abs.S
        wA = which(abs.S > lambda1[m] & abs.S <= a * 
          lambda1[m])
        w[wA] = -(abs.S[wA]^2 - 2 * a * lambda1[m] * 
          abs.S[wA] + lambda1[m]^2)/2/(a - 1)
        wA = which(abs.S > a * lambda1[m])
        w[wA] = (a + 1)^2 * lambda1[m]^2/2
        w = matrix(w, d, d)
      }
      pen = pen + sum((w - diag(diag(w))) * abs(theta$sigma[[m]]))
      omega = matrix(0, d, d)
      if (lambda2[m] > 0) {
        abs.A = abs(theta$A[[m]][[1]])
        omega = lambda2[m] * abs.A
        wA = which(abs.A > lambda2[m] & abs.A <= a * 
          lambda2[m])
        omega[wA] = -(abs.A[wA]^2 - 2 * a * lambda2[m] * 
          abs.A[wA] + lambda2[m]^2)/2/(a - 1)
        wA = which(abs.A > a * lambda2[m])
        omega[wA] = (a + 1)^2 * lambda2[m]^2/2
        omega = matrix(omega, d, d)
      }
      pen = pen + sum(omega * abs(theta$A[[m]][[1]]))
    }
    ll.pen = (FB$loglik - ((T - 1) * N.samples) * pen)
  }
  res = list(theta = theta, ll_history = ll_history, Iter = cnt, 
    Npar = Npar, BIC = BIC, smoothedprob = FB$probS, ll.pen = ll.pen)
  class(res) <- "MSAR"
  res$call = cl
  res
}


Estep.MSAR <- function (data, theta, smth = FALSE, verbose = FALSE, covar.emis = covar.emis, 
  covar.trans = covar.trans) 
{
  go <- Sys.time()
  if (verbose) {
    message("Starting LWS_pfilt")
  }
  att.theta = attributes(theta)
  label <- att.theta$label
  p <- att.theta$order
  if (verbose) 
    print(theta)
  M <- att.theta$NbRegimes
  d <- att.theta$NbComp
  if (is.null(d) || is.na(d)) {
    d = 1
  }
  n_par <- att.theta$n_par
  order = att.theta$order
  data <- as.array(data)
  T = dim(data)[1]
  if (is.null(T)) {
    T = length(data)
  }
  if (abs(length(data)/T - trunc(length(data)/T)) > 1e-05) {
    stop("error : size of data should be nxT.sample with n integer")
  }
  N.samples = dim(data)[2]
  if (is.null(N.samples) || is.na(N.samples)) {
    N.samples <- 1
  }
  data <- array(data, c(T, N.samples, d))
  gamma <- array(0, c(N.samples, T - p, M))
  xi = array(0, c(M, M, T - (p + 1), N.samples))
  loglik = 0
  if (substr(label, 2, 2) == "N") {
    ncov.emis = dim(covar.emis)[3]
    if (is.null(ncov.emis) || is.na(ncov.emis)) {
      ncov.emis = 1
    }
    covar.emis = array(covar.emis, c(T, N.samples, ncov.emis))
  }
  if (substr(label, 1, 1) == "H") {
    tmat = as.matrix(theta$transmat)
    prior = as.matrix(theta$prior)
    for (ex in 1:N.samples) {
      if (verbose) {
        print(ex)
      }
      g <- emisprob.MSAR(data[, ex, ], theta = theta, 
        covar = covar.emis[, ex, ])
      FB = forwards_backwards(prior, tmat, g)
      gamma[ex, , ] = t(FB$gamma)
      xi[, , , ex] = FB$xi
      loglik = loglik + FB$loglik
    }
  }
  else {
    if (missing(covar.trans)) {
      stop("error : covariable is missing")
    }
    if (length(covar.trans) == 1) {
      Lag = covar.trans + 1
      covar.trans = array(data[(1):(T - Lag + 1), , ], 
        c(T - Lag + 1, N.samples, d))
      data = array(data[Lag:T, , ], c(T - Lag + 1, N.samples, 
        d))
    }
    ncov.trans = dim(covar.trans)[3]
    if (is.null(ncov.trans) || is.na(ncov.trans)) {
      ncov.trans = 1
    }
    ct = array(0, c(T, N.samples, ncov.trans))
    ct[1:min(dim(ct)[1], dim(covar.trans)[1]), , ] = covar.trans[1:min(dim(ct)[1], 
      dim(covar.trans)[1]), , ]
    covar.trans = ct
    for (ex in 1:N.samples) {
      if (verbose) {
        print(ex)
      }
      g <- emisprob.MSAR(data[, ex, ], theta = theta, 
        covar = covar.emis[, ex, ])
      transmat = theta$transmat
      par.trans = theta$par.trans
      nh_transition = attributes(theta)$nh.transitions
      inp = covar.trans[(order + 1):T, ex, ]
      transmat.t = nh_transition(array(inp, c((T - order), 
        1, ncov.trans)), par.trans, transmat)
      FB = nhforwards_backwards(theta$prior, transmat.t, 
        g)
      gamma[ex, , ] = t(FB$gamma)
      xi[, , , ex] = FB$xi
      loglik = loglik + FB$loglik
    }
  }
  list(loglik = loglik, probS = gamma, probSS = xi, M = FB$M)
}

emisprob.MSAR <- function (data, theta, covar = NULL) 
{
  w = NULL
  d <- attributes(theta)$NbComp
  M <- attributes(theta)$NbRegimes
  order <- attributes(theta)$order
  label <- attributes(theta)$label
  data = as.matrix(data)
  A = NULL
  if (d == 1) {
    if (order > 0) {
      for (j in 1:M) {
        A[[j]] = sapply(theta$A[j, ], as.numeric)
      }
    }
    else {
      A[[1]] <- array(0, c(1, 1, M))
    }
    Sigma = matrix(theta$sigma, M, 1)
    A0 = matrix(theta$A0, M, 1)
  }
  else {
    if (order > 0) {
      A = theta$A
    }
    else {
      A <- matrix(0, M, d)
    }
    Sigma = theta$sigma
    A0 = theta$A0
  }
  T <- dim(as.matrix(data))[1]
  prob <- matrix(1, T - order, M)
  for (j in 1:M) {
    A0.emis = array(0, c(d, T - order))
    if (order > 0) {
      for (o in (1:order)) {
        A0.emis = A0.emis + matrix(A[[j]][[o]], d, d) %*% 
          t(data[((order + 1):T) - o, ])
      }
    }
    A0.emis = A0.emis + A0[j, ]
    if (substr(label, 2, 2) == "N") {
      nh.emissions <- attributes(theta)$nh.emissions
      par.emis = theta$par.emis
      ncov = dim(covar)[2]
      if (is.null(ncov) & length(covar) > 0) {
        covar = matrix(covar, length(covar), 1)
        ncov = 1
      }
      for (i in 1:d) {
        femis = nh.emissions(matrix(covar[(order + 1):T, 
          ], T - order, ncov), as.matrix(par.emis[[j]]))
        A0.emis[i, ] = A0.emis[i, ] + femis[i, ]
      }
    }
    if (!is.na(sum(data))) {
      prob[, j] = pdf.norm(t(data[((order + 1):T), ]), 
        A0.emis, matrix(Sigma[[j]], d, d))
    }
    else {
      moy = t(data[((order + 1):T), ]) - A0.emis
      w = which(is.na(apply(moy, 2, sum)))
      prob[-w, j] = pdf.norm(matrix(moy[, -w], d, T - 
        order - length(w)), 0, matrix(Sigma[[j]], d, 
        d))
    }
  }
  prob = t(prob)
  prob
}

nhforwards_backwards <- function (prior, transition, obslik, filter_only = 0) 
{
  T = dim(obslik)[2]
  Q1 = length(prior)
  scale = matrix(1, 1, T)
  loglik = 0
  alpha = matrix(0, Q1, T)
  gamma = matrix(0, Q1, T)
  xi = array(0, c(Q1, Q1, T - 1))
  t = 1
  alpha[, 1] = prior * obslik[, t]
  scale[t] = sum(alpha[, t])
  alpha[, t] = normalise(alpha[, t])
  for (t in 2:T) {
    tmp = (t(transition[, , t]) %*% alpha[, t - 1]) * obslik[, 
      t]
    scale[t] = sum(tmp)
    alpha[, t] = normalise(tmp)
  }
  loglik = sum(log(scale))
  beta = matrix(0, Q1, T)
  gamma = matrix(0, Q1, T)
  beta[, T] = matrix(1, Q1, 1)
  gamma[, T] = normalise(alpha[, T] * beta[, T])
  for (t in seq(T - 1, 1, -1)) {
    b = beta[, t + 1] * obslik[, t + 1]
    beta[, t] = normalise((transition[, , t + 1] %*% b))
    gamma[, t] = normalise(alpha[, t] * beta[, t])
    xi[, , t] = normalise((transition[, , t + 1] * (alpha[, 
      t] %*% t(b))))
  }
  FB = NULL
  FB$gamma = gamma
  FB$xi = xi
  FB$loglik = loglik
  FB$M = Q1
  FB$alpha = alpha
  FB$beta = beta
  return(FB)
}


EM_converged <- function (loglik, previous_loglik, threshold = 1e-04) 
{
  converged = 0
  decrease = 0
  if (!(previous_loglik == -Inf)) {
    if (loglik - previous_loglik < -0.01) {
      print(paste("******likelihood decreased from ", 
        previous_loglik, " to ", loglik, sep = ""), 
        quote = FALSE)
      decrease = 1
    }
    delta_loglik = abs(loglik - previous_loglik)
    avg_loglik = (abs(loglik) + abs(previous_loglik) + threshold)/2
    bb = ((delta_loglik/avg_loglik) < threshold)
    if (bb) {
      converged = 1
    }
  }
  res <- NULL
  res$converged <- converged
  res$decrease <- decrease
  return(res)
}


Mstep.hh.reduct.MSAR <- function (data, theta, FB, sigma.diag = FALSE) 
{
  T = dim(data)[1]
  N.samples = dim(as.array(data))[2]
  d = dim(as.array(data))[3]
  if (is.null(d) | is.na(d)) {
    d = 1
  }
  M <- attributes(theta)$NbRegimes
  p <- attributes(theta)$order
  order <- max(p, 1)
  data = array(data, c(T, N.samples, d))
  data2 = array(0, c(order * d, T - order + 1, N.samples))
  cpt = 1
  for (o in order:1) {
    for (kd in 1:d) {
      data2[cpt, , ] = data[o:(T - order + o), , kd]
      cpt = cpt + 1
    }
  }
  T = length(o:(dim(data)[1] - order + o))
  exp_num_trans = 0
  exp_num_visit = 0
  exp_num_visits1 = 0
  postmix = 0
  m = matrix(0, d * order, M)
  m_1 = m
  c = matrix(0, M, 1)
  s = c
  op = array(0, c(d * order, d * order, M))
  op_1 = op
  op_2 = op_1
  for (ex in 1:N.samples) {
    obs = array(data2[, , ex], c(d * order, T))
    xit = array(FB$probSS[, , , ex], c(M, M, T - 2))
    gamma = matrix(FB$probS[ex, 1:(T - 1), ], T - 1, M)
    exp_num_trans = exp_num_trans + apply(xit, c(1, 2), 
      sum)
    exp_num_visits1 = exp_num_visits1 + gamma[1, ]
    postmix = postmix + apply(gamma, 2, sum)
    obs = t(obs)
    for (j in 1:M) {
      w = matrix(gamma[, j], 1, T - 1)
      wobs = obs[1:(T - 1), ] * repmat(t(w), 1, d * order)
      wobs_1 = obs[2:T, ] * repmat(t(w), 1, d * order)
      if (is.na(sum(obs))) {
        wobs[is.na(wobs)] = 0
        wobs_1[is.na(wobs_1)] = 0
        obs[is.na(obs)] = 0
      }
      m[, j] = m[, j] + apply(wobs, 2, sum)
      m_1[, j] = m_1[, j] + apply(wobs_1, 2, sum)
      op[, , j] = op[, , j] + t(wobs) %*% obs[1:(T - 1), 
        ]
      op_1[, , j] = op_1[, , j] + t(wobs) %*% obs[2:T, 
        ]
      op_2[, , j] = op_2[, , j] + t(wobs_1) %*% obs[2:T, 
        ]
    }
  }
  prior = normalise(exp_num_visits1)
  prior = matrix(prior, M, 1)
  transmat = mk_stochastic(exp_num_trans)
  if (min(postmix) < 1e-06) {
    stop("error : smoothing probabilities are to small, in one regime at least. You should revise initialisation.")
  }
  moy <- array(0, c(M, d))
  sigma <- list()
  A2 <- list()
  for (j in 1:M) {
    A2[[j]] = list()
    Cxx = (postmix[j] * op[, , j] - m[, j] %*% t(m[, j]))/postmix[j]^2
    Cxy = (postmix[j] * op_1[, , j] - m[, j] %*% t(m_1[, 
      j]))/postmix[j]^2
    Cyy = (postmix[j] * op_2[, , j] - m_1[, j] %*% t(m_1[, 
      j]))/postmix[j]^2
    A.th = theta$A[[j]][[1]]
    S.th = theta$sigma[[j]]
    SA.th = Cyy - Cxy %*% A.th - A.th %*% Cxy + A.th %*% 
      Cxx %*% t(A.th)
    ll0 = -sum(diag(SA.th %*% solve(S.th)))
    cnt = 0
    wA = which(A.th != 0)
    lA = length(wA)
    A = matrix(0, d * d, d * d)
    b = matrix(Cxy, d * d, 1)
    lwi = 0
    lwj = 0
    for (id in 1:d) {
      wi = which(A.th[id, ] != 0)
      for (jd in 1:d) {
        cnt = cnt + 1
        A[cnt, lwi + (1:length(wi))] = Cxx[jd, wi]
        A[cnt, lA + lwj + (1:(d - length(wi)))] = -S.th[-wi, 
          jd]/2
      }
      lwj = lwj + (d - length(wi))
      lwi = lwi + length(wi)
    }
    A2.lasso = matrix(0, d, d)
    tmp = solve(A, b)[1:lA]
    lwi = 0
    for (id in 1:d) {
      wi = which(A.th[id, ] != 0)
      A2.lasso[id, wi] = tmp[lwi + (1:length(wi))]
      lwi = lwi + length(wi)
    }
    tmp = (m_1[, j] - (A2.lasso) %*% m[, j])/postmix[j]
    tmp2 = Cyy + A2.lasso %*% Cxx %*% t(A2.lasso) - (A2.lasso %*% 
      Cxy + t(A2.lasso %*% Cxy))
    S2.lasso = Cyy - Cxy %*% A2.lasso - A2.lasso %*% Cxy + 
      A2.lasso %*% Cxx %*% t(A2.lasso)
    ll1 = -sum(diag(S2.lasso %*% solve(tmp2)))
    tmp = (m_1[, j] - (A2.lasso) %*% m[, j])/postmix[j]
    tmp2 = Cyy + A2.lasso %*% Cxx %*% t(A2.lasso) - (A2.lasso %*% 
      Cxy + t(A2.lasso %*% Cxy))
    A2[[j]][[1]] = A2.lasso
    if (sigma.diag) {
      sigma[[j]] = diag(diag(tmp2[1:d, 1:d]), d)
    }
    else {
      w = which(abs(theta$sigma[[j]]) > 0)
      sigma[[j]] = matrix(0, d, d)
      sigma[[j]][w] = tmp2[1:d, 1:d][w]
    }
    moy[j, 1:d] = tmp[1:d]
  }
  if (p > 0) {
    list(A = A2, A0 = moy, sigma = sigma, prior = prior, 
      transmat = transmat)
  }
}


normalise <- function (M) 
{
  c = sum(M)
  d = c + (c == 0)
  M = M/d
  return(M)
}


mk_stochastic <- function (T) 
{
  if (is.vector(T)) {
    T = normalise(T)
  }
  else {
    n = length(dim(T))
    normaliser = apply(T, 1, sum)
    normaliser = matrix(normaliser, dim(T)[1], dim(T)[2])
    normaliser = normaliser + (normaliser == 0)
    T = T/normaliser
  }
  return(T)
}



Mstep.hh.lasso.MSAR <- function (data, theta, FB) 
{
  T = dim(data)[1]
  N.samples = dim(as.array(data))[2]
  d = dim(as.array(data))[3]
  if (is.null(d) | is.na(d)) {
    d = 1
  }
  M <- attributes(theta)$NbRegimes
  p <- attributes(theta)$order
  order <- max(p, 1)
  data = array(data, c(T, N.samples, d))
  data2 = array(0, c(order * d, T - order + 1, N.samples))
  cpt = 1
  for (o in order:1) {
    for (kd in 1:d) {
      data2[cpt, , ] = data[o:(T - order + o), , kd]
      cpt = cpt + 1
    }
  }
  T = length(o:(dim(data)[1] - order + o))
  exp_num_trans = 0
  exp_num_visit = 0
  exp_num_visits1 = 0
  postmix = 0
  m = matrix(0, d * order, M)
  m_1 = m
  c = matrix(0, M, 1)
  s = c
  op = array(0, c(d * order, d * order, M))
  op_1 = op
  op_2 = op_1
  wx = array(0, c(N.samples * (T - 1), d, M))
  wy = array(0, c(N.samples * (T - 1), d, M))
  for (ex in 1:N.samples) {
    obs = array(data2[, , ex], c(d * order, T))
    xit = array(FB$probSS[, , , ex], c(M, M, T - 2))
    gamma = matrix(FB$probS[ex, 1:(T - 1), ], T - 1, M)
    exp_num_trans = exp_num_trans + apply(xit, c(1, 2), 
      sum)
    exp_num_visits1 = exp_num_visits1 + gamma[1, ]
    postmix = postmix + apply(gamma, 2, sum)
    obs = t(obs)
    for (j in 1:M) {
      w = matrix(gamma[, j], 1, T - 1)
      wx[(ex - 1) * (T - 1) + (1:(T - 1)), , j] = obs[1:(T - 
        1), ] * repmat(t(w), 1, d * order)
      wy[(ex - 1) * (T - 1) + (1:(T - 1)), , j] = obs[2:T, 
        ] * repmat(t(w), 1, d * order)
      wobs = obs[1:(T - 1), ] * repmat(t(w), 1, d * order)
      wobs_1 = obs[2:T, ] * repmat(t(w), 1, d * order)
      if (is.na(sum(obs))) {
        wobs[is.na(wobs)] = 0
        wobs_1[is.na(wobs_1)] = 0
        obs[is.na(obs)] = 0
      }
      m[, j] = m[, j] + apply(wobs, 2, sum)
      m_1[, j] = m_1[, j] + apply(wobs_1, 2, sum)
      op[, , j] = op[, , j] + t(wobs) %*% obs[1:(T - 1), 
        ]
      op_1[, , j] = op_1[, , j] + t(wobs) %*% obs[2:T, 
        ]
      op_2[, , j] = op_2[, , j] + t(wobs_1) %*% obs[2:T, 
        ]
    }
  }
  prior = normalise(exp_num_visits1)
  prior = matrix(prior, M, 1)
  transmat = mk_stochastic(exp_num_trans)
  if (min(postmix) < 1e-06) {
    stop("error : smoothing probabilities are to small, in one regime at least. You should revise initialisation.")
  }
  moy <- array(0, c(M, d))
  sigma <- list()
  A2 <- list()
  for (j in 1:M) {
    ll.old = ll_gauss.MSAR(data, theta, FB$probS, regime = j)
    bic.old = -2 * ll.old + log(length(data)) * (d + sum(theta$A[[j]][[1]] != 
      0) + d * d)
    A2[[j]] = list()
    A2.lasso = matrix(0, d, d)
    residuals = matrix(0, N.samples * (T - 1), d)
    Cxx = postmix[j] * op[, , j] - m[, j] %*% t(m[, j])
    Cxy = postmix[j] * op_1[, , j] - m[, j] %*% t(m_1[, 
      j])
    for (id in 1:d) {
      lars.1 = lars(wx[, , j], wy[, id, j])
      coef = predict(lars.1, wx[, , j], typ = "coef")$coefficients
      BIC.lm = NULL
      mylm = list()
      w = list()
      for (kst in 2:length(lars.1$lambda)) {
        w[[kst]] = which(coef[kst, ] > 1e-05)
        mylm[[kst]] = lm(wy[, id, j] ~ wx[, w[[kst]], 
          j])
        BIC.lm[kst] = -2 * sum(log(pdf.norm(matrix(c(wy[, 
          id, 1]), 1, length(wy)), matrix(c(mylm[[kst]]$fitted.values), 
          1, length(wy)), as.matrix(var(wy[, id, j] - 
          mylm[[kst]]$fitted.values))))) + log(length(wy)) * 
          (length(w[[kst]]) + 2)
      }
      w.bic = which.min(BIC.lm)
      Cxx_w = Cxx[w[[w.bic]], w[[w.bic]]]
      Cxy_w = Cxy[w[[w.bic]], id]
      A2.lasso[id, w[[w.bic]]] = t(Cxy_w) %*% solve(Cxx_w)
    }
    tmp = (m_1[, j] - (A2.lasso) %*% m[, j])/postmix[j]
    op_1[, , j] = t(as.matrix(op_1[, , j]))
    tmp2 = (op_2[, , j] + (A2.lasso) %*% op[, , j] %*% t(A2.lasso) - 
      ((A2.lasso) %*% t(op_1[, , j]) + t((A2.lasso) %*% 
        t(op_1[, , j]))))/postmix[j] - tmp %*% t(tmp)
    theta.tmp = theta
    theta.tmp$A[[j]][[1]] = A2.lasso
    theta.tmp$sigma[[j]] = tmp2[1:d, 1:d]
    theta.tmp$A0[j, ] = tmp[1:d]
    ll.new = ll_gauss.MSAR(data, theta.tmp, FB$probS, regime = j)
    bic.new = -2 * ll.new + log(length(data)) * (d + sum(A2.lasso != 
      0) + d * d)
    print(c("........Regime ", j, ", ll.old ", ll.old, ", ll.new ", 
      ll.new), quote = FALSE)
    print(c("                , bic.old ", bic.old, ", bic.new ", 
      bic.new), quote = FALSE)
    if (bic.new > bic.old) {
      for (id in 1:d) {
        w = which(theta$A[[j]][[1]][id, ] != 0)
        Cxx_w = Cxx[w, w]
        Cxy_w = Cxy[w, id]
        A2.lasso[id, w] = t(Cxy_w) %*% solve(Cxx_w)
      }
      tmp = (m_1[, j] - (A2.lasso) %*% m[, j])/postmix[j]
      op_1[, , j] = t(as.matrix(op_1[, , j]))
      tmp2 = (op_2[, , j] + (A2.lasso) %*% op[, , j] %*% 
        t(A2.lasso) - ((A2.lasso) %*% t(op_1[, , j]) + 
        t((A2.lasso) %*% t(op_1[, , j]))))/postmix[j] - 
        tmp %*% t(tmp)
    }
    A2[[j]][[1]] = A2.lasso
    sigma[[j]] = tmp2[1:d, 1:d]
    moy[j, 1:d] = tmp[1:d]
  }
  if (p > 0) {
    list(A = A2, A0 = moy, sigma = sigma, prior = prior, 
      transmat = transmat)
  }
}

ll_gauss.MSAR <- function (data, theta, gamma, regime) 
{
  T = dim(data)[1]
  if (is.null(T)) {
    T = length(data)
  }
  N.samples = dim(data)[2]
  if (is.null(N.samples) || is.na(N.samples)) {
    N.samples = 1
  }
  d = dim(as.array(data))[3]
  if (is.null(d) || is.na(d)) {
    d = 1
  }
  data = array(data, c(T, N.samples, d))
  order = attributes(theta)$order
  gamma = matrix(gamma[, , regime], N.samples, dim(gamma)[2])
  A0 = theta$A0[regime, ]
  A = list()
  A[[1]] = theta$A[[regime]][[1]]
  sigma = theta$sigma[[regime]]
  f = 0
  for (ex in 1:N.samples) {
    ar = array(0, c(d, T - order))
    A0.emis = matrix(0, d, T - order)
    for (o in 1:order) {
      ar = ar + A[[o]] %*% t(data[((order + 1):T) - o, 
        ex, ])
    }
    for (i in 1:d) {
      A0.emis[i, ] = ar[i, ] + A0[i]
    }
    dn = pdf.norm(t(data[(order + 1):T, ex, ]), A0.emis, 
      sigma)
    dn[dn < 1e-10] = 1e-10
    B = log(dn)
    f = f + sum(gamma[ex, ] * B[1:length(gamma[ex, ])])
  }
  f = -f
  return(f)
}


as.thetaMSAR <- function (x, label = "HH", regime_names = NULL, ncov.emis = ncov.emis, 
  ncov.trans = ncov.trans) 
{
  if (!is.thetaMSAR(x)) {
    stop("as.thetaMSAR: your input is not like a theta at all")
  }
  dimname <- dimnames(x)
  att <- attributes(x)
  att$dimnames <- NULL
  att$names <- NULL
  att$class <- NULL
  label = att$label
  M = att$NbRegimes
  order = att$order
  d = att$NbComp
  dimname = c("A0", "sigma", "prior", "transmat")
  n_par = NULL
  if (order > 0) {
    dimname = c("A", dimname)
  }
  if (substr(label, 1, 1) == "N") {
    dimname = c(dimname, "par.trans")
    n_par = n_par + length(c(x$par.trans))/M
  }
  if (substr(label, 2, 2) == "N") {
    dimname = c(dimname, "par.emis")
    n_par = n_par + length(c(x$par.emis))/M
  }
  names(x) = dimname
  x$A0 = matrix(x$A0, M, d)
  rownames(x$A0) = c(paste("Regime", 1:M, sep = ""))
  colnames(x$A0) = c(paste("A0", 1:d, sep = ""))
  x$prior = matrix(x$prior, M, d)
  x$prior = matrix(x$prior, M, 1)
  rownames(x$prior) = c(paste("Regime", 1:M, sep = ""))
  colnames(x$prior) = ""
  rownames(x$transmat) = c(paste("Regime", 1:M, sep = ""))
  colnames(x$transmat) = c(paste("Regime", 1:M, sep = ""))
  if (substr(label, 2, 2) == "N") {
    names(x$par.emis) = c(paste("Regime", 1:M, sep = ""))
    for (i in 1:M) {
      if (!is.null(dim(x$par.emis[[i]]))) {
        rownames(x$par.emis[[i]]) = rep("", d)
      }
      if (!is.null(dim(x$par.emis[[i]])[2])) {
        colnames(x$par.emis[[i]]) = c(paste("coef.emis", 
          1:dim(x$par.emis[[i]])[2], sep = ""))
      }
    }
  }
  if (substr(label, 1, 1) == "N") {
    rownames(x$par.trans) = c(paste("Regime", 1:M, sep = ""))
    colnames(x$par.trans) = c(paste("coef.trans", 1:max(2, 
      ncov.trans + 1), sep = ""))
  }
  if (d == 1) {
    x$sigma = matrix(x$sigma, M, d)
    rownames(x$sigma) = c(paste("Regime", 1:M, sep = ""))
    colnames(x$sigma) = ""
  }
  else {
    for (i in 1:M) {
      colnames(x$sigma[[i]]) = rep("", d)
      rownames(x$sigma[[i]]) = rep("", d)
    }
  }
  if (order > 0 && d == 1) {
    rownames(x$A) = c(paste("Regime", 1:M, sep = ""))
    colnames(x$A) = c(paste("A", 1:order, sep = ""))
  }
  else if (order > 0) {
    names(x$A) = c(paste("Regime", 1:M, sep = ""))
    for (i in 1:M) {
      names(x$A[[i]]) = c(paste("A", 1:order, sep = ""))
      for (j in 1:order) {
        colnames(x$A[[i]][[j]]) = rep("", d)
        rownames(x$A[[i]][[j]]) = rep("", d)
      }
    }
    names(x$sigma) = c(paste("Regime", 1:M, sep = ""))
  }
  if (is.null(att$order)) {
    att$order <- order
  }
  if (is.null(att$n_par)) {
    att$n_par <- order * d^2 + d + d^2 + 1 + M
  }
  for (a in names(att)) {
    attr(x, a) <- att[[a]]
  }
  class(x) <- "MSAR"
  return(x)
}

is.thetaMSAR <- function (x) 
{
  if (!is.list(x)) {
    return(FALSE)
  }
  if (length(x) < 4) {
    FALSE
  }
  else {
    TRUE
  }
}




library(pracma)
library(lars)

Mstep.hh.lasso.MSAR_patched <- function (data, theta, FB) 
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
  prior = NHMSAR:::normalise(exp_num_visits1)
  prior = matrix(prior, M, 1)
  transmat = NHMSAR:::mk_stochastic(exp_num_trans)
  if (min(postmix) < 1e-06) {
    stop("error : smoothing probabilities are to small, in one regime at least. You should revise initialisation.")
  }
  moy <- array(0, c(M, d))
  sigma <- list()
  A2 <- list()
  for (j in 1:M) {
    ll.old = NHMSAR:::ll_gauss.MSAR(data, theta, FB$probS, regime = j)
    bic.old = -2 * ll.old + log(length(data)) * (d + sum(theta$A[[j]][[1]] != 
                                                           0) + d * d)
    A2[[j]] = list()
    A2.lasso = matrix(0, d, d)
    residuals = matrix(0, N.samples * (T - 1), d)
    Cxx = postmix[j] * op[, , j] - m[, j] %*% t(m[, j])
    Cxy = postmix[j] * op_1[, , j] - m[, j] %*% t(m_1[, 
                                                      j])
    for (id in 1:d) {
      lars.1 = lars(wx[, , j], wy[, id, j], type = "lasso", normalize = TRUE)
      coef = predict(lars.1, wx[, , j], typ = "coef")$coefficients
      BIC.lm = NULL
      mylm = list()
      w = list()
      
      for (kst in 2:length(lars.1$lambda)) {
        w[[kst]] = which(abs(coef[kst, ]) > 1e-05) # Changed to check the absolute values of coeffs
        
        if ( length( w[[kst]] ) == 0 ) { 
          # alle coef <= 1e-05 => skip case !!!
          # w[[kst]] = NA
        }
        else
        {
          mylm[[kst]] = lm(wy[, id, j] ~ wx[, w[[kst]], j])
          BIC.lm[kst] = -2 * sum(log(NHMSAR:::pdf.norm(
            matrix(c(wy[,id, 1]), 1, length(wy)),
            matrix(c(mylm[[kst]]$fitted.values), 1, length(wy)), 
            as.matrix(var(wy[, id, j] - mylm[[kst]]$fitted.values))
          ))) + log(length(wy)) *
            (length(w[[kst]]) + 2)
        }
      }
      
      if ( is.null(BIC.lm) ) {
        # alle coef <= 1e-05 => skip case !!!
      } else {
        w.bic = which.min(BIC.lm)
        Cxx_w = Cxx[w[[w.bic]], w[[w.bic]]]
        Cxy_w = Cxy[w[[w.bic]], id]
        A2.lasso[id, w[[w.bic]]] = t(Cxy_w) %*% solve(Cxx_w)
      }
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
    ll.new = NHMSAR:::ll_gauss.MSAR(data, theta.tmp, FB$probS, regime = j)
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

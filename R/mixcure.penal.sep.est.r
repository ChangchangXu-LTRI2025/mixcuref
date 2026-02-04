##########################################################################################
####  Mixcure model with two sets of covariates for FT-PL or ML via direct comparison  ###
##########################################################################################


mixcure.penal.sep.est <- function(formula, k.cure, k.surv, data, init, pl, iterlim = 200) {
require(splines)
require(survival)
require(abind)


  #########################################################################################

  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];

  design.matrix <- model.matrix(formula, data = design.matrix);

  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix);
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;

#############################################################

  loglik.mixture.part <- function(
    p, survt, design.matrix1, design.matrix0,
    index.cure.var = index.cure.v, index.surv.var = index.surv.v, pl) {
    t      <- survt[, 1]
    status <- survt[, 2]
    event  <- (status == 1L)
    cens   <- !event
    logt   <- log(t)

    design.mtx.comb <- cbind(design.matrix0, design.matrix1)

    # gamma (assumes index.gamma exists in caller env)
    gamma <- p[index.gamma - (length(k.cure) + length(k.surv))]
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    X0 <- as.matrix(design.matrix0)
    X1 <- as.matrix(design.matrix1)

    # ---- theta & eps (same branching logic) ----
    # lp_cure <- drop(as.matrix(design.matrix)[, index.cure.var, drop = FALSE] %*%
    #                     as.matrix(p[!(p %in% c(index.surv.var - length(k.cure),
    #                                            index.gamma - length(k.cure)))]))
    len_X0 <- length(index.cure.var)
    len_X1 <- length(index.surv.var)

    lp_cure <- drop(as.matrix(design.matrix)[, index.cure.var, drop = FALSE] %*%
                      as.matrix(p[1:len_X0]))
    theta <- plogis(lp_cure)

      lp_surv <- drop(design.mtx.comb[, index.surv.var, drop = FALSE] %*%
                        as.matrix(p[(len_X0+1):(len_X0+len_X1)]))
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)

    # Stable rewrite: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    eta   <- emeps / D
    delta <- (1 - theta) * eta

    # IMPORTANT: keep kap exactly (note the signs)
    kap <- -theta * (1 - theta) * (1 - eta) + (1 - theta)^2 * eta * (1 - eta)

    # pi = exp(eps)*eps*eta^2
    pi <- eps * emeps / (D * D)

    # ---- negative log-likelihood ----
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))
    if (isFALSE(pl)) return(loglikelihood)

    # =========================================================================================
    # Fast Fisher blocks
    # =========================================================================================
    max.len <- ncol(design.matrix)

    # Mirror your loop indexing: i,j treated as 1..max.len columns of design.matrix0/1
    X0A <- X0[, 1:max.len, drop = FALSE]
    X1A <- X1[, 1:max.len, drop = FALSE]

    e <- as.numeric(event)
    c <- 1 - e

    # --- Block A (subset to index.cure.var)
    wA <- e * (theta * (1 - theta)) + c * kap
    A_full <- crossprod(X0A, X0A * wA)
    info.a <- A_full[index.cure.var, index.cure.var, drop = FALSE]

    # --- Block B: your loop uses j in (index.surv.var-max.len, max.len+1)
    # weight: eps*(1-delta)*delta on cens only
    Xt0 <- cbind(X0A, logt)  # (n x (max.len+1))
    wB <- c * (eps * (1 - delta) * delta)

    B_full <- -crossprod(X1A, Xt0 * wB)  # (max.len x (max.len+1))
    cols_B <- c(index.surv.var - max.len, index.gamma - max.len)
    info.b <- B_full[index.cure.var, cols_B, drop = FALSE]

    # --- Block D: design.xt1 = [design.matrix1, logt]
    Xt1 <- cbind(X1A, logt)  # (n x (max.len+1))

    # cens term in your code: eps*delta - eps^2*(delta*(1-delta))
    wd2 <- eps * delta - (eps^2) * (delta * (1 - delta))
    wD <- e * eps + c * wd2

    D_full <- crossprod(Xt1, Xt1 * wD)
    D_full[max.len + 1, max.len + 1] <- D_full[max.len + 1, max.len + 1] + sum(event) / (gamma * gamma)

    rowscols_D <- c(index.surv.var - max.len, index.gamma - max.len)
    info.d <- D_full[rowscols_D, rowscols_D, drop = FALSE]

    # Schur complement without explicit inverse
    tmp <- solve(info.d, t(info.b))

    # The original has a special-case orientation when length(index.cure.v) < 3 & part.cure == TRUE
    # if (length(index.cure.v) < 3 && isTRUE(part.cure)) {
    #   info.set0 <- info.a - t(info.b) %*% tmp  # t(B) D^{-1} B
    # } else {
      info.set0 <- info.a - info.b %*% tmp    # B D^{-1} t(B)
    #}

    # Stable log(det(info.set0)*det(info.d))
    detS <- determinant(info.set0, logarithm = TRUE)
    detD <- determinant(info.d, logarithm = TRUE)
    if (detS$sign * detD$sign <= 0) return(Inf)

    logdet <- as.numeric(detS$modulus) + as.numeric(detD$modulus)
    loglikelihood - 0.5 * logdet
  }

  #################################################################

  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for each of cure part variables;

    is = k.surv + length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p = init[-c(k.cure,is)],
      survt = survt, design.matrix1 = design.matrix,
      design.matrix0=design.matrix,
      index.cure.var=index.cure.v[-k.cure],
      index.surv.var=index.surv.v[-k.surv],
      pl=pl, iterlim = iterlim, hessian=T
    );

  hessmat <- maximizer$hessian
  var.mat <- solve(hessmat)
  alpha.hat <- maximizer$estimate[length(maximizer$estimate)];
  loglik <- -maximizer$minimum
  z.score <- maximizer$estimate / sqrt(diag(var.mat));

  index.cure.v1 = index.cure.v[-k.cure]
  index.surv.v1 = index.surv.v[-k.surv]
  len_cure = length(index.cure.v1)
  len_surv = length(index.surv.v1)
  coef.table.cure <- cbind(
    'coef'        = maximizer$estimate[1:len_cure],
    'exp(coef)'   = exp(maximizer$estimate[1:len_cure]),
    'se(coef)'    = sqrt(diag(var.mat)[1:len_cure]),
    'z'           = z.score[1:len_cure],
    'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[1:len_cure]))),
    'LCI.95%' = maximizer$estimate[1:len_cure] - 1.96 * sqrt(diag(var.mat)[1:len_cure]),
    'UCI.95%' = maximizer$estimate[1:len_cure] + 1.96 * sqrt(diag(var.mat)[1:len_cure])
  );
  rownames(coef.table.cure) <- colnames(design.matrix)[index.cure.v1];

  coef.table.surv <- cbind(
    'coef'        = maximizer$estimate[(len_cure+1):(len_cure+len_surv)],
    'exp(coef)'   = exp(maximizer$estimate[(len_cure+1):(len_cure+len_surv)]),
    'se(coef)'    = sqrt(diag(var.mat)[(len_cure+1):(len_cure+len_surv)]),
    'z'           = z.score[(len_cure+1):(len_cure+len_surv)],
    'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[(len_cure+1):(len_cure+len_surv)]))),
    'LCI.95%' = maximizer$estimate[(len_cure+1):(len_cure+len_surv)] - 1.96 * sqrt(diag(var.mat)[(len_cure):(len_cure+len_surv)]),
    'UCI.95%' = maximizer$estimate[(len_cure+1):(len_cure+len_surv)] + 1.96 * sqrt(diag(var.mat)[(len_cure):(len_cure+len_surv)])
  );
  rownames(coef.table.surv) <- colnames(design.matrix)[index.surv.v1];

  coef.table.alpha <- cbind(
    'coef'     = alpha.hat,
    'se(coef)' = sqrt(diag(var.mat)[(len_cure+len_surv+1)]),
    'z'        = z.score[(len_cure+len_surv+1)],
    'Pr(>|z|)' = 2 * (1 - pnorm(abs(z.score[(len_cure+len_surv+1)]))),
    'LCI.95%'  = maximizer$estimate[(len_cure+len_surv+1)] - 1.96 * sqrt(diag(var.mat)[(len_cure+len_surv+1)]),
    'UCI.95%'  = maximizer$estimate[(len_cure+len_surv+1)] + 1.96 * sqrt(diag(var.mat)[(len_cure+len_surv+1)]),
    'loglik' = -maximizer$minimum
  );
  rownames(coef.table.alpha) <- 'alpha';


  #######################################
  ## Output tables from either method; ##
  #######################################

  colnames(var.mat) <- c(
    paste('cure.', colnames(design.matrix)[index.cure.v1]),
    paste('surv.', colnames(design.matrix)[c(index.surv.v1-ncol(design.matrix))]),
    'alpha'
  );
  rownames(var.mat) <- colnames(var.mat);

  out <- list(
    coefficients = list(
      cure = coef.table.cure,
      surv = coef.table.surv,
      alpha = coef.table.alpha
      #   run.time
    ),
    cov = var.mat
  );
  class(out) <- c('mixcure', 'list');

  return(out);

}



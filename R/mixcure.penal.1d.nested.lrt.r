#####################################################################
### Function for ikelihood ratio test (df=1) mixcure model         ##
#### via the nested deviance method under penalized loglikelihoods ##
#####################################################################
###################################################
#### previously 'mixcure.penal.ESTV.nested1d.r'  ##
###################################################

mixcure.penal.1d.nested.lrt <- function(formula, data, init, pl, loglik, iterlim = 200) {
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

  samp.s <- nrow(design.matrix)


#############################################################
#############################################################
#############################################################

  #create CI for profile likelihood, this option only outputs estimates and PL under specified model;


  ##loglik function for testing parameters of cure or surv part;
  loglik.mixture.part <- function(
    p, survt, design.matrix1, design.matrix0,
    index.cure.var = index.cure.v,
    index.surv.var = index.surv.v,
    pl
  ) {
    t      <- survt[, 1]
    status <- survt[, 2]
    event  <- (status == 1L)
    cens   <- !event
    logt   <- log(t)

    design.mtx.comb <- cbind(design.matrix0, design.matrix1)

    # gamma (assumes index.gamma exists in caller env)
    gamma <- p[index.gamma - 1]
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    X0 <- as.matrix(design.matrix0)
    X1 <- as.matrix(design.matrix1)

    # ---- theta & eps ----
    if (k <= length(index.cure.v)) {
      lp_cure <- drop(design.mtx.comb[, index.cure.var, drop = FALSE] %*%
                        as.matrix(p[index.cure.v[-length(index.cure.v)]]))
      theta <- plogis(lp_cure)

      lp_surv <- drop(design.mtx.comb[, index.surv.var, drop = FALSE] %*%
                        as.matrix(p[index.surv.var - 1]))
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    } else {
      lp_cure <- drop(design.mtx.comb[, index.cure.var, drop = FALSE] %*%
                        as.matrix(p[index.cure.var]))
      theta <- plogis(lp_cure)

      lp_surv <- drop(design.mtx.comb[, index.surv.var, drop = FALSE] %*%
                        as.matrix(p[index.surv.v[-length(index.surv.v)]]))
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    }

    # Stable rewrite: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    eta   <- emeps / D
    delta <- (1 - theta) * eta

    # IMPORTANT: keep kap exactly as in your function
    kap <- theta * (1 - theta) * (1 - eta) - (1 - theta)^2 * eta * (1 - eta)

    # pi = exp(eps)*eps*eta^2 (rewrite safely)
    pi <- eps * emeps / (D * D)

    # ---- negative log-likelihood ----
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))
    if (!isTRUE(pl)) return(loglikelihood)

    # =========================================================================================
    # Fast Fisher blocks
    # =========================================================================================
    max.len <- max(length(index.cure.var), length(index.surv.var))

    # Mirror loop indexing: i,j treated as 1..max.len columns
    X0A <- X0[, 1:max.len, drop = FALSE]
    X1A <- X1[, 1:max.len, drop = FALSE]

    e <- as.numeric(event)
    c <- 1 - e

    # --- Block A: full (max.len x max.len)
    wA <- e * (theta * (1 - theta)) + c * kap
    info.a <- crossprod(X0A, X0A * wA)

    # --- Block B: full (max.len x (max.len+1)), weight theta(1-theta)*pi on censored only
    Xt0 <- cbind(X0A, logt)  # (n x (max.len+1))
    wB <- c * (theta * (1 - theta) * pi)
    info.b <- -crossprod(X1A, Xt0 * wB)

    # --- Block D: full ((max.len+1) x (max.len+1))
    Xt1 <- cbind(X1A, logt)  # (n x (max.len+1))
    wd2 <- eps * delta - (eps^2) * delta + (eps^2) * (delta^2)
    wD <- e * eps + c * wd2

    info.d <- crossprod(Xt1, Xt1 * wD)
    info.d[max.len + 1, max.len + 1] <- info.d[max.len + 1, max.len + 1] + sum(event) / (gamma * gamma)

    # Schur complement without explicit inverse
    tmp <- solve(info.d, t(info.b))
    info.set0 <- info.a - info.b %*% tmp

    # Stable log(det(info.set0)*det(info.d))
    detS <- determinant(info.set0, logarithm = TRUE)
    detD <- determinant(info.d, logarithm = TRUE)
    if (detS$sign * detD$sign <= 0) return(Inf)

    logdet <- as.numeric(detS$modulus) + as.numeric(detD$modulus)

    loglikelihood - 0.5 * logdet
  }


  #################################################################

  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for each cure part variable;

  dim.v <- ncol(design.matrix)
  ll.cure <- rep(0,dim.v)
  llr.cure <- rep(0,dim.v)
  pval.cure <- rep(0,dim.v)
# index.cure.v[-1] for no intercept calculation
  for (k in index.cure.v) {
    maximizer <- nlm(
      f = loglik.mixture.part, p = init[-k],
      survt = survt, design.matrix0 = design.matrix,
      design.matrix1=design.matrix,
      index.cure.var=index.cure.v[-k], pl=pl,
      iterlim = iterlim, hessian=T
    );
    loglik.part = -maximizer$minimum;
    dif.ll = -2*(loglik.part-loglik);  #loglik is ll under Ha;
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.cure[k]<- loglik.part
    llr.cure[k]<- dif.ll
    pval.cure[k]<- pval
    if (det(maximizer$hessian) < 1e-05)
      diag(maximizer$hessian) <- diag(maximizer$hessian) + 1e-06

  }



  ### loglikelihood calculation for each surv part variable;

  ll.surv <- rep(0,ncol(design.matrix))
  llr.surv <- rep(0,ncol(design.matrix))
  pval.surv <- rep(0,ncol(design.matrix))


  for (k in index.surv.v) {
    is=k-length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p = init[-k],
      survt = survt, design.matrix1 = design.matrix,
      design.matrix0=design.matrix,
      index.surv.var=index.surv.v[-is], pl=pl,
      iterlim = iterlim, hessian=FALSE
    );

    loglik.part = -maximizer$minimum;
    dif.ll = -2*(loglik.part-loglik);
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.surv[is]<- loglik.part
    llr.surv[is]<-dif.ll
    pval.surv[is]<-pval
  }



  coef.table.cure <- cbind(
    'LL.cure' = ll.cure,
    'LLR'         = llr.cure,
    'Pr(>chisq)'  = pval.cure
  );
  rownames(coef.table.cure) <- colnames(design.matrix);

  coef.table.surv <- cbind(
      'LL.surv' = ll.surv,
    'LLR'         = llr.surv,
    'Pr(>chisq)'  = pval.surv
  );
  rownames(coef.table.surv) <- colnames(design.matrix);

  coef.table.alpha <- "NA";


#run.time = proc.time() - init.time

out <- list(
  coefficients = list(
    cure = coef.table.cure,
    surv = coef.table.surv,
    alpha = coef.table.alpha
 #   run.time
  )
);
class(out) <- c('mixcure.1d.nested.lrt', 'list');

return(out);

}


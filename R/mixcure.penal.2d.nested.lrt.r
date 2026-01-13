#####################################################################
### Function for ikelihood ratio test (df=2) mixcure model         ##
#### via the nested deviance method under penalized loglikelihoods ##
#####################################################################
#### previously 'mixcure.penal.ESTV.nested.r' ##
################################################

mixcure.penal.2d.nested.lrt <- function(formula, data, init, pl, loglik, iterlim = 200) {
require(splines)
require(survival)
require(abind)
  # require(foreach)
  # require(parallel)
  # require(doSNOW)

  #########################################################################################

  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];

  design.matrix <- model.matrix(formula, data = design.matrix);

  # index ranges of coefficients of glm and cox models
  index.cure.vt <- 1 : ncol(design.matrix);
  index.surv.vt <- (ncol(design.matrix) + 1) : (2*length(index.cure.vt))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.vt)+1;

  #samp.s <- nrow(design.matrix)


  ####################################################
  ## nonlinear minimization algoritm to solve       ##
  ## penalized mixture cure loglikelihood functions ##
  ####################################################

  loglik.mixture <- function(
    p, survt, design.matrix,
    index.cure.var, index.cure.v,
    index.surv.var, index.surv.v,
    pl
  ) {
    t      <- survt[, 1]
    status <- survt[, 2]
    event  <- (status == 1L)
    cens   <- !event
    logt   <- log(t)

    X <- as.matrix(design.matrix)

    # gamma (as in your code)
    gamma <- p[index.gamma]
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    # ---- theta & eps for the (unpenalised) likelihood (same branching intent) ----
    # Original uses design.matrix[,-k] and scalar vs %*%; preserve that behaviour.
    Xlik <- X[, -k, drop = FALSE]

    if (ncol(design.matrix) <= 2) {
      lp_cure <- drop(Xlik * p[index.cure.var])
      theta <- plogis(lp_cure)

      lp_surv <- drop(Xlik * p[index.surv.var])
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    } else {
      lp_cure <- drop(Xlik %*% p[index.cure.var])
      theta <- plogis(lp_cure)

      lp_surv <- drop(Xlik %*% p[index.surv.var])
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    }

    # Stable rewrite for likelihood: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    # negative log-likelihood (same expression, stabilised)
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))

    if (!isTRUE(pl)) return(loglikelihood)

    # =========================================================================================
    # Penalised part: keep the logic of zeroing two parameters, then recompute on full design
    # =========================================================================================
    p2 <- p
    p2[c(k, k + length(index.cure.vt))] <- 0  # uses your the k and index.cure.vt

    # linear predictors on FULL design.matrix
    lp_cure_pl <- drop(X %*% p2[index.cure.v])
    theta <- plogis(lp_cure_pl)

    lp_surv_pl <- drop(X %*% p2[index.surv.v])
    logeps <- gamma * logt + lp_surv_pl
    eps <- exp(logeps)

    # Stable rewrite for eta/delta/pi
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    eta   <- emeps / D
    delta <- (1 - theta) * eta

    # kap and pi exactly as your pl-block intends (kap for est/PLCI; B uses eps*(1-delta)*delta)
    kap <- (1 - eta) * (1 - theta) * (theta + eta)
    pi  <- eps * emeps / (D * D)  # stable form of exp(eps)*eps*eta^2

    # =========================================================================================
    # Fast Fisher blocks (replace nested loops with crossprod)
    # =========================================================================================
    kc <- length(index.cure.v)
    ks <- length(index.surv.v)

    e <- as.numeric(event)
    c <- 1 - e

    # Xc and Xt for block calculations
    Xc <- X[, index.cure.v, drop = FALSE]          # columns referenced by i in your loops
    Xt <- cbind(Xc, logt)                          # design.xt (note: uses cure cols + logt)

    # --- Block A (kc x kc): event weight theta(1-theta), cens weight kap
    wA <- e * (theta * (1 - theta)) + c * kap
    info.a <- crossprod(Xc, Xc * wA)

    # --- Block B (kc x (kc+1)): cens weight eps*(1-delta)*delta
    wB <- c * (eps * (1 - delta) * delta)
    info.b <- -crossprod(Xc, Xt * wB)

    # --- Block D ((kc+1) x (kc+1)): event weight eps; cens weight eps*delta - eps^2*delta + eps^2*delta^2
    wd2 <- eps * delta - (eps^2) * delta + (eps^2) * (delta^2)
    wD <- e * eps + c * wd2
    info.d <- crossprod(Xt, Xt * wD)

    # Add gamma-gamma term (bottom-right element corresponds to logt column)
    info.d[kc + 1, kc + 1] <- info.d[kc + 1, kc + 1] + sum(event) / (gamma * gamma)

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

  ######END of loglik.mixture####################################

  dim.v <- ncol(design.matrix)
  est.list <- matrix(0,ncol=dim.v, nrow = (2*dim.v+2))


  for (k in index.cure.vt) {

  maximizer <- nlm(
    f = loglik.mixture,
    p = init,
    survt=survt, design.matrix=design.matrix,
    index.cure.var=index.cure.vt[-k],
    index.surv.var=index.surv.vt[-k],
    index.cure.v=index.cure.vt,
    index.surv.v=index.surv.vt,
    pl = pl,
    iterlim = iterlim, hessian=F);


loglik.part <- -maximizer$minimum #in loglik function loglik was calculated as minus of actual loglik value
est.value <- maximizer$estimate
llr <- 2*(loglik - loglik.part)
pval <- pchisq(llr, df =2 , lower.tail = F)
est.list[,k] <- c(est.value[-c(k,(k + dim.v))], loglik.part, llr, pval)
  }

rownames(est.list) <- c(colnames(design.matrix)[-k],colnames(design.matrix)[-k],"shape","loglik.part", "llr", "pval")
coef.table <-  est.list[-c(1:(2*dim.v-1)),]
colnames(coef.table) <- colnames(design.matrix);


out <- list(
  coefficients = coef.table
  #cov = var.mat
);
class(out) <- c('mixcure.2d.nested.lrt', 'list');

return(out);
}


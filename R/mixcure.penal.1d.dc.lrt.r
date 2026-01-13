#####################################################################
####      Function for LRT of mixcure model for              ########
####      FT-PL or ML via direct likelihood ratio comparison ########
#####################################################################
#### previously 'mixcure.penal.lrt.test.r' ##########
#####################################################

mixcure.penal.1d.dc.lrt <- function(formula, data, init, pl, iterlim = 200) {
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
#### loglik function of full model

  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl, PLCI=F) {

    ####  parameter and variable dep parameters;
    t <- survt[, 1];  status <- survt[, 2];  event <- (status == 1L);  cens  <- !event;  logt  <- log(t)
    # Sub-matrices for the two linear predictors
    Xc <- design.matrix[, index.cure.var, drop = FALSE]
    k.vec <- ncol(Xc)                       # number of covariates in each part
    Xs <- design.matrix[, (index.surv.var - k.vec), drop = FALSE]
    index.gamma <- 2 * k.vec + 1            # matches parameter layout
    gamma <- p[index.gamma]

    # Guard (optional but helps nlm avoid NaNs)
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    # Cure part: theta
    lp_cure <- drop(Xc %*% p[index.cure.var])
    theta <- plogis(lp_cure)
    # Survival part: eps = t^gamma * exp(Xs beta) = exp(gamma*logt + Xs beta)
    lp_surv <- drop(Xs %*% p[index.surv.var])
    logeps <- gamma * logt + lp_surv
    eps <- exp(logeps)
    # Common denominator: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    # Negative log-likelihood
    # event term: log(1-theta) + log(gamma) - log(t) + log(eps) - eps
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))

    if (!pl) return(loglikelihood)

    # Quantities for penalty
    eta   <- emeps / D
    delta <- (1 - theta) * eta
    kap   <- (1 - eta) * (1 - theta) * (theta + eta)
    pi    <- eps * emeps / (D * D)      # equals exp(eps)*eps*eta^2 but overflow-safe
    # Build Xt once: survival covariates + log(t) for gamma column
    Xt <- cbind(Xs, logt)
    # Convert event/cens to 0/1 for weighting without subsetting
    e <- as.numeric(event)
    c <- 1 - e

    # --- Block A (k x k): for cure by cure parameters
    w1 <- theta * (1 - theta)
    wA <- e * w1 + c * kap
    info.a <- crossprod(Xc, Xc * wA)

    # --- Block B (k x (k+1)): for cure by surv parameters, only censored contribute
    wB <- c * (w1 * pi)
    info.b <- -crossprod(Xc, Xt * wB)

    # --- Block D ((k+1) x (k+1)): for surv by surv parameters
    wd2 <- eps * delta - eps^2 * delta + eps^2 * delta^2
    wD <- e * eps + c * wd2
    info.d <- crossprod(Xt, Xt * wD)

    # Add your extra gamma-gamma term
    info.d[k.vec + 1, k.vec + 1] <- info.d[k.vec + 1, k.vec + 1] + sum(event) / (gamma * gamma)

    # Schur complement without explicit inverse:
    # S = A - B %*% D^{-1} %*% t(B)
    tmp <- solve(info.d, t(info.b))     # (k+1) x k
    S <- info.a - info.b %*% tmp

    # log|det| in a stable way
    detS <- determinant(S, logarithm = TRUE)
    detD <- determinant(info.d, logarithm = TRUE)

    if (!PLCI) {if (detS$sign * detD$sign <= 0) return(Inf)}
    logdet <- as.numeric(detS$modulus) + as.numeric(detD$modulus)
    loglikelihood - 0.5 * logdet
  }

  ######END of loglik.mixture####################################


  # Parameter estimation under Ha (non-restricted likelihood)
  # maximize penalized or unpenalized loglikelihood by nlm;
  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix,
    pl = pl,
    iterlim = iterlim, hessian=F);

  loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value


#############################################################


  ##loglik function for testing parameters of cure or surv part;
  #design.matrix1-surv part, design.matrix0-cure part;

  loglik.mixture.part <- function(
    p, survt, design.matrix1, design.matrix0,
    index.cure.var = index.cure.v, index.surv.var = index.surv.v,
    pl, part.cure = FALSE
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

    # ---- theta & eps (same branching logic as your function) ----
    if (k > length(index.cure.v)) {
      lp_cure <- drop(as.matrix(design.matrix)[, index.cure.var, drop = FALSE] %*%
                        as.matrix(p[index.cure.var]))
      theta <- plogis(lp_cure)

      lp_surv <- drop(design.mtx.comb[, index.surv.var, drop = FALSE] %*%
                        as.matrix(p[-c(index.cure.var, index.gamma - 1)]))
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    } else {
      lp_cure <- drop(as.matrix(design.matrix)[, index.cure.var, drop = FALSE] %*%
                        as.matrix(p[-c(index.surv.var - 1, index.gamma - 1)]))
      theta <- plogis(lp_cure)

      lp_surv <- drop(design.mtx.comb[, index.surv.var, drop = FALSE] %*%
                        as.matrix(p[index.surv.var - 1]))
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    }

    # Stable rewrite: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    eta   <- emeps / D
    delta <- (1 - theta) * eta

    # IMPORTANT: keep your kap exactly (note the signs)
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
    max.len <- max(length(index.cure.var), length(index.surv.var))

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

    # Your original has a special-case orientation when length(index.cure.v) < 3 & part.cure == TRUE
    if (length(index.cure.v) < 3 && isTRUE(part.cure)) {
      info.set0 <- info.a - t(info.b) %*% tmp  # t(B) D^{-1} B
    } else {
      info.set0 <- info.a - info.b %*% tmp    # B D^{-1} t(B)
    }

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

  dim.v <- ncol(design.matrix)
  ll.cure <- rep(0,dim.v)
  llr.cure <- rep(0,dim.v)
  pval.cure <- rep(0,dim.v)

  for (k in index.cure.v) {
      maximizer <- nlm(
      f = loglik.mixture.part,
      p = init[-k],
      survt = survt, design.matrix0 = design.matrix,
      design.matrix1=design.matrix, part.cure = T,
      index.cure.var=index.cure.v[-k],
      pl=pl,
      iterlim = iterlim, hessian=F
    );

    loglik.part = -maximizer$minimum;
    dif.ll = -2*(loglik.part-loglik);  #loglik is ll under non-restricted model;
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.cure[k]<- loglik.part
    llr.cure[k]<- dif.ll
    pval.cure[k]<- pval

  }


  ### loglikelihood calculation for each surv part variable;

  ll.surv <- rep(0,ncol(design.matrix))
  llr.surv <- rep(0,ncol(design.matrix))
  pval.surv <- rep(0,ncol(design.matrix))


  for (k in index.surv.v) {
    is=k-length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p =  init[-k],
      survt = survt, design.matrix1 = design.matrix,
      design.matrix0=design.matrix, part.cure = F,
      index.surv.var=index.surv.v[-is],
      pl=pl,
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


#######################################
## Output tables from FT-PLE or ML; ##
#######################################

out <- list(
  coefficients = list(
    cure = coef.table.cure,
    surv = coef.table.surv,
    alpha = coef.table.alpha
 #  run.time
  )
);
class(out) <- c('mixcure.dc.lrt', 'list');

return(out);

}



################################
### Function for mixcure model##
#### penalized loglikelihoods ##
###################################################
#### Last modified Dec31 2018 for one x variable ##
###################################################

########### NOTE on Oct22:
########## Error checking for cure uppper points iter2;
########## maximizer$temp1 is not consistent in and out of iter2 loop;

mixcure.penal.profile.CI <- function(formula, data, init, pl, apct = 0.05, LRT.pval = F, iterlim = 200) {
require(splines)
require(survival)
require(abind)
require(R.utils)


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


  ####################################################
  ## nonlinear minimization algoritm to solve       ##
  ## penalized mixture cure loglikelihood functions ##
  ####################################################

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

  #############################################################################################

  ## loglik function for constructing profile likelihood of cure or surv part ##
  loglik.mixture.profile <- function(p, survt, k=k, design.matrix1=design.matrix, design.matrix0=design.matrix, param.est, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {
    t      <- survt[, 1]
    status <- survt[, 2]
    event  <- (status == 1L)
    cens   <- !event
    logt   <- log(t)

    # Combined design used in eps construction
    design.mtx.comb <- cbind(design.matrix0, design.matrix1)
    ik <- k - length(index.cure.var)

    # Parameter indices used in your original code
    index.gamma <- 2 * ncol(design.matrix0) + 1  # (only used as "index.gamma-1" here)
    gamma <- p[index.gamma - 1]
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    # ---- theta (profile handling exactly as your original logic) ----
    if (k > length(index.cure.v)) {
      lp_cure <- drop(design.matrix0[, index.cure.var, drop = FALSE] %*%
                        as.matrix(p[index.cure.var]))
    } else {
      lp_cure <- drop(design.matrix0[, index.cure.var[-k], drop = FALSE] %*%
                        as.matrix(p[-c(index.surv.var - 1, index.gamma - 1)]) -
                        design.mtx.comb[, k] * param.est)
    }
    theta <- plogis(lp_cure)

    # ---- eps (profile handling exactly as your original logic) ----
    if (k > length(index.cure.v)) {
      lp_surv <- drop(design.mtx.comb[, index.surv.var[-ik], drop = FALSE] %*%
                        as.matrix(p[-c(index.cure.var, index.gamma - 1)]) +
                        design.mtx.comb[, k] * param.est)
    } else {
      lp_surv <- drop(design.mtx.comb[, index.surv.var, drop = FALSE] %*%
                        as.matrix(p[index.surv.var - 1]))
    }

    logeps <- gamma * logt + lp_surv
    eps <- exp(logeps)

    # Stable common terms
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps   # theta + (1-theta)*exp(-eps)
    eta   <- emeps / D
    delta <- (1 - theta) * eta
    # kap is same as your code (for est, PLCI)
    kap <- (1 - eta) * (1 - theta) * (theta + eta)
    # pi = exp(eps)*eps*eta^2  (rewrite safely)
    pi <- eps * emeps / (D * D)
    # Unpenalised negative log-likelihood (same as your code)
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))

    if (!isTRUE(pl)) return(loglikelihood)

    # =========================================================================================
    # Penalisation piece: build Fisher blocks using weighted crossproducts
    # =========================================================================================

    max.len <- max(length(index.cure.var), length(index.surv.var))

    # We will embed the actual matrices into max.len-sized placeholders exactly like your loops did,
    # then subset to the same rows/cols you were using.
    # Cure design (block A) uses design.matrix0
    X0 <- as.matrix(design.matrix0)
    # Survival design for block D uses design.matrix1 and log(t)
    X1 <- as.matrix(design.matrix1)

    # --- Block A (max.len x max.len)
    # weights: event -> theta*(1-theta); cens -> kap
    wA <- as.numeric(event) * (theta * (1 - theta)) + as.numeric(cens) * kap
    # Use first max.len columns as in your indexing regime
    X0A <- X0[, 1:max.len, drop = FALSE]
    A_full <- crossprod(X0A, X0A * wA)
    info.a <- A_full[index.cure.var, index.cure.var, drop = FALSE]

    # --- Block B (max.len x (max.len+1)),
    Xt0 <- cbind(X0A, logt)  # (n x (max.len+1))
    wB <- as.numeric(cens) * (eps * (1 - delta) * delta)
    # Left side uses design.matrix1[, i] with i in 1:max.len
    X1B <- X1[, 1:max.len, drop = FALSE]
    B_full <- -crossprod(X1B, Xt0 * wB)  # (max.len x (max.len+1))
    # Column mapping from your original code
    cols_B <- c(index.surv.var - max.len, index.gamma - max.len)
    info.b <- B_full[index.cure.var, cols_B, drop = FALSE]

    # --- Block D (max.len+1 x max.len+1),
    # design.xt1 = [design.matrix1, logt]
    Xt1 <- cbind(X1B, logt)  # (n x (max.len+1))
    wd2 <- eps * delta - (eps^2) * delta + (eps^2) * (delta^2)
    wD <- as.numeric(event) * eps + as.numeric(cens) * wd2
    D_full <- crossprod(Xt1, Xt1 * wD)

    # Add your gamma-gamma term the FULL block (position max.len+1)
    D_full[max.len + 1, max.len + 1] <- D_full[max.len + 1, max.len + 1] + sum(event) / (gamma * gamma)
    rowscols_D <- c(index.surv.var - max.len, index.gamma - max.len)
    info.d <- D_full[rowscols_D, rowscols_D, drop = FALSE]

    # Schur complement without explicit inverse:
    # info.set0 = A - B D^{-1} B^T
    tmp <- solve(info.d, t(info.b))
    info.set0 <- info.a - info.b %*% tmp

    # log(det(info.set0)*det(info.d)) safely
    detS <- determinant(info.set0, logarithm = TRUE)
    detD <- determinant(info.d, logarithm = TRUE)

    # If signs go bad, return Inf so nlm backs off
    if (detS$sign * detD$sign <= 0) return(Inf)
    logdet <- as.numeric(detS$modulus) + as.numeric(detD$modulus)
    loglikelihood - 0.5 * logdet
  }


  #################################################################
  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for each cure part variable;

  dim.v <- ncol(design.matrix)

  if (LRT.pval == T) {
    ll.cure <- rep(0,dim.v)
    llr.cure <- rep(0,dim.v)
    pval.cure <- rep(0,dim.v)
    #ll.est.cure <- array(0,c((2*dim.v+1),(2*dim.v+1)))
    #varmat.A.cure <-sample(0,3*dim.v,replace = TRUE)
    #varmat.A.cure <-array(0,c((2*dim.v+1),(2*dim.v+1),(2*dim.v+1))) #vcov matrix of each single parameter (cure.var1,cure.var2) likelihood ratio test;

    # dim(varmat.A.cure) = c(dim.v,dim.v,dim.v)
   # score.U.cure <- array(0,c((2*dim.v+1),(2*dim.v+1))) #gradient vector for all variables of each single parameter(cure.var1,cure.var2) likelihood ratio test;
    # init = c(0.866409603,	-0.334994711,	0.004421239,	-0.070689979,	-8.11466538,	-0.349162116,	0.143310052,	-0.049019853,	0.1)

    for (k in index.cure.v[-1]) {
   # mle under the reduced (null) model for cure parameter;
    maximizer <- nlm(
      f = loglik.mixture.part, p = init,
      survt = survt, design.matrix0 = design.matrix,
      design.matrix1=design.matrix,
      index.cure.var=index.cure.v[-k], pl=pl,
      iterlim = iterlim, hessian=F
    );
   # ll.est.cure[,k] <- maximizer$estimate
    loglik.part = -maximizer$minimum;
    dif.ll = -2*(loglik.part-loglik);  #loglik is ll under Ha;
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.cure[k]<- loglik.part
    llr.cure[k]<- dif.ll
    pval.cure[k]<- pval
    if (det(maximizer$hessian) < 1e-05)
      diag(maximizer$hessian) <- diag(maximizer$hessian) + 1e-06
   # varmat.A.cure[,,k] <- solve(maximizer$hessian)  #if hessian matrix contains row/col of zeros, add small values;
   # score.U.cure [,k] <- maximizer$gradient

  }
}
    ###################################
    # Profile likelihood CI endpoint  #
    ###################################

  # Note:
  # loglik     -- loglikelihood of all the parameters under MLE of full likelihood;
  # l.up       -- loglik corresponds to upper or lower CI bounds B0=B+delta or B0=B-delta;
  # tol        -- tolerance level for defining loglik difference of estimated and actual parameters are converged;
  # lambda     -- lambda quantity for profile likelhood calculation as referenced;
  # delta.up   -- delta quantity for upper endpoint of profile likelihood calculation as referenced;
  # delta.lo   -- delta quantity for lower endpoint of profile likelihood calculation as referenced;
  # l.temp/l0.b-- loglik of estimated endpoint parameter values under profile likelihood for the corresponding parameter;

  ######################################################
  ## By parameter upper or lower endpoint calculation ##
  ######################################################
  # apct=0.05

  l.null = loglik - 0.5 * qchisq(1-apct,df=1,ncp = 0,lower.tail=T)
  ni = 1

  ######## Cure part variable CI endpoints ########
  #################################################

  upper.cure <- rep(0,dim.v)
  lower.cure <- rep(0,dim.v)
 # for (k in index.cure.v) {

    for (k in index.cure.v) {

  ##################upper endpoint##########################
   # tol = 0.2
    tol = 0.1
 l.temp <- loglik

 n=ni+1
 # iter0 <- 1; converge = FALSE;l0.b.up = 0;
 # while(l.temp > l.null & (!converge| is.nan(l0.b.up)) & iter0<=30) {

   #assign initial values to parameter estimates
   param.est.up <- maximizer0$estimate

   converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.up = 0;delta.up = 0
    while ((!converge| is.nan(l0.b.up)) & iter1 <= 25 & !EXIT1 & !is.nan(delta.up)) {

             # calculate log-lik, score and hessian under l0.b;
      maximizer.temp <-  nlm(
        f = loglik.mixture, p = param.est.up, survt=survt, design.matrix=design.matrix,
        pl = pl, iterlim = 1, hessian=TRUE)
      score.temp = maximizer.temp$gradient
      hessian.temp = maximizer.temp$hessian
      if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06

      #### Approach 1: lambda = (2*(l0.b-l.up+e*A^-1*U)/(e*A^-1*e))^0.5
      # l0.b.up <- -loglik.mixture(p=param.est.up, survt, design.matrix,
      #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
      l0.b.up <- -maximizer.temp$minimum
      inv.hessian.temp <- solve(hessian.temp)

      if ((l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k]<0) {
      #   lambda <- (0.5*(l0.b.up - l.null) + (-inv.hessian.temp %*% score.temp)[k])/inv.hessian.temp[k,k]
        lambda <- (inv.hessian.temp %*% score.temp)[k]/inv.hessian.temp[k,k]
       # #define increment for estimated value: delta=-A^-1(U-lambda*e)
      delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
       } else{
      lambda <- (2*(l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
      #define increment for estimated value: delta=-A^-1(U-lambda*e)
      delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)
      };

      # maximizing loop for unpenalized estimates;
      #if (pl == F) {
      inside <- FALSE; iter2 <- 1;
      while (!inside & iter2 <= 100 & !is.nan(delta.up)) {

        # add increment to stepwise parameter value;
        param.est.temp.up <- param.est.up
        param.est.temp.up[k] <- param.est.temp.up[k]+delta.up
       # if (k==2) {param.est.temp.up[1] <- -1;param.est.temp.up[9] <- 0.1}

        #compute loglikelihood function using updated parameter values;

        maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.up[-k], survt=survt,
                                param.est = param.est.temp.up[k], k = k,
                                pl = pl, iterlim = iterlim, hessian=TRUE)
        l.temp.up = -maximizer.temp1$minimum

        #if (!is.nan(l.temp.up))

         #compare to see if updated l is still
          inside <- (l.temp.up > (l.null - 0.05)) #l.null - 0.05 for all others, 0.2 for k=3 of high rate H0
          #diff.up = l.temp.up - l.null
          #converge0 <- (abs(diff.up) <= tol)
          alevel.up <- pchisq(2*(l.temp-l.temp.up),df=1,ncp=0,lower.tail = T)
         # print(c(delta.up, alevel.up, n,l.temp.up,k,iter1,iter2))
          if (!inside) {delta.up <- delta.up/((n+1)/n);iter2 <- iter2 + 1}  #(n+0.1)/n for low rate H0;
          if (is.nan(delta.up)) {param.est.temp.up[k] <- NA}
       } #for iter2

        #}
      #Using converged increment for parameter to get corresponding score and variance expressions;
      param.est.up <- insert(maximizer.temp1$estimate, ats=k, values=param.est.temp.up[k])

      l0.b.up = l.temp.up

      diff.up = l0.b.up - l.null
      converge <- (abs(diff.up) <= tol)
      if ((!converge| is.nan(l0.b.up)) & !is.nan(delta.up)) {iter1 <- iter1 + 1; n = n + 1} else {EXIT1 = T;}
      if (is.nan(delta.up)==T) {param.est.up[k] <- NA}
    } #for iter1
  #} #for iter0
    upper.cure[k] <- param.est.up[k]


  ###############lower endpoint#####################

    n=ni
 #iter0 <- 1; converge = FALSE; l0.b.lo=0
 #while(l.temp > l.null & (!converge| is.nan(l0.b.lo)) & iter0<=30) {

   #assign initial values to parameter estimates
   param.est.lo <- maximizer0$estimate

   converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.lo=0; delta.lo=0
   while ((!converge| is.nan(l0.b.lo)) & iter1 <= 25 & !EXIT1 & !is.nan(delta.lo)) {

     # calculate log-lik, score and hessian under l0.b;
     maximizer.temp <-  nlm(
       f = loglik.mixture, p = param.est.lo, survt=survt, design.matrix=design.matrix,
       pl = pl, iterlim = 1, hessian=TRUE)
     score.temp = maximizer.temp$gradient
     hessian.temp = maximizer.temp$hessian
     if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06

     #### Approach 1: lambda = (2*(l0.b-l.null+e*A^-1*U)/(e*A^-1*e))^0.5
     # l0.b.lo <- -loglik.mixture(p=param.est.lo, survt, design.matrix,
     #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
     l0.b.lo <- -maximizer.temp$minimum
     inv.hessian.temp <- solve(hessian.temp)
     if ((l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k]<0) {
     #if ((l0.b.lo < l.null + 0.5 * score.temp %*% inv.hessian.temp %*% score.temp)|(l0.b.lo - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k]<0) {
       lambda <- ((-inv.hessian.temp %*% score.temp)[k])/inv.hessian.temp[k,k]
       delta.lo <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
       #define increment for estimated value: delta=-A^-1(U-lambda*e)
     } else{
       lambda <- -(2*(l0.b.lo - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
       #define increment for estimated value: delta=-A^-1(U-lambda*e)
       delta.lo <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)};

       # maximizing loop for unpenalized estimates;
       #if (pl == F) {
       inside <- FALSE; iter2 <- 1;
       while (!inside & iter2 <= 100 & !is.nan(delta.lo)) {

         # add increment to stepwise parameter value;
         param.est.temp.lo <- param.est.lo
         param.est.temp.lo[k] <- param.est.temp.lo[k] + delta.lo
        # param.est.temp.lo[9] <- 0.1

         #compute loglikelihood function using lodated parameter values;
         maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.lo[-k], survt=survt,
                                 param.est = param.est.temp.lo[k], k = k,
                                 pl = pl, iterlim = iterlim, hessian=TRUE)
         l.temp.lo = -maximizer.temp1$minimum

        if (k==3) {dt.lo=0.1} else {dt.lo=0.5}
           inside <- (l.temp.lo > l.null - 0.1)
          # diff.lo = l.temp.lo - l.null
          # converge0 <- (abs(diff.lo) <= tol)
           alevel.lo <- pchisq(2*(l.temp-l.temp.lo),df=1,ncp=0,lower.tail = T)
          # print(c(delta.lo, alevel.lo, n,l.temp.lo,k,iter1,iter2))
           if (!inside) {delta.lo <- delta.lo/((n+dt.lo)/n);iter2 <- iter2 + 1} #for variables other than LuminalA (n+0.5)/n;
           if (is.nan(delta.lo)) {param.est.temp.lo[k] <- NA}

       } # for iter2;
       param.est.lo <- insert(maximizer.temp1$estimate, ats=k, values=param.est.temp.lo[k])
          l0.b.lo = l.temp.lo

       diff.lo = l0.b.lo - l.null
       converge <- (abs(diff.lo) <= tol)
       if ((!converge | is.nan(l0.b.lo)) & !is.nan(delta.lo)) {iter1 <- iter1 + 1; n = n + 2} else {EXIT1 = T}
       if (is.nan(delta.lo)==T) {param.est.lo[k] <- NA}

   } #for iter1

   lower.cure[k] <- param.est.lo[k]

 }



  ### loglikelihood calculation for each surv part variable;


  if (LRT.pval == T) {
    ll.surv <- rep(0,ncol(design.matrix))
    llr.surv <- rep(0,ncol(design.matrix))
    pval.surv <- rep(0,ncol(design.matrix))

     for (k in index.surv.v[-1]) {
    # mle under the reduced (null) model for surv parameter;
    is=k-length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p = init,
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
  }

  ######## Surv part variable CI endpoints ########
  #################################################

  upper.surv <- rep(0,dim.v)
  lower.surv <- rep(0,dim.v)
  for (k in index.surv.v) {

    is=k-length(index.cure.v)
    ##################upper endpoint##########################
    l.temp <- loglik
    tol = 0.1

    n=ni
    # iter0 <- 1; converge = FALSE;l0.b.up = 0;
    # while(l.temp > l.null & (!converge| is.nan(l0.b.up)) & iter0<=30) {

    #assign initial values to parameter estimates
    param.est.up <- maximizer0$estimate

    converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.up = 0;delta.up = 0
    while ((!converge| is.nan(l0.b.up)) & iter1 <= 25 & !EXIT1 & !is.nan(delta.up)) {

      # calculate log-lik, score and hessian under l0.b;
      maximizer.temp <-  nlm(
        f = loglik.mixture, p = param.est.up, survt=survt, design.matrix=design.matrix,
        pl = pl, iterlim = 1, hessian=TRUE)
      score.temp = maximizer.temp$gradient
      hessian.temp = maximizer.temp$hessian
      if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06

      #### Approach 1: lambda = (2*(l0.b-l.up+e*A^-1*U)/(e*A^-1*e))^0.5
      # l0.b.up <- -loglik.mixture(p=param.est.up, survt, design.matrix,
      #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
      l0.b.up <- -maximizer.temp$minimum
      inv.hessian.temp <- solve(hessian.temp)
      if ((l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k]<0) {
        lambda <- (inv.hessian.temp %*% score.temp)[k]/inv.hessian.temp[k,k]
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
      } else{
        lambda <- (2*(l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)};

      # maximizing loop for unpenalized estimates;
      #if (pl == F) {
      inside <- FALSE; iter2 <- 1;
      while (!inside & iter2 <= 100 & !is.nan(delta.up)) {

        # add increment to stepwise parameter value;
        param.est.temp.up <- param.est.up
        param.est.temp.up[k] <- param.est.temp.up[k]+delta.up

        #compute loglikelihood function using updated parameter values;

        maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.up[-k], survt=survt,
                                param.est = param.est.temp.up[k], k = k,
                                pl = pl, iterlim = iterlim, hessian=TRUE)
        l.temp.up = -maximizer.temp1$minimum

        #if (!is.nan(l.temp.up))

        #compare to see if updated l is still
        inside <- (l.temp.up > (l.null - 0.1))
        #diff.up = l.temp.up - l.null
        #converge0 <- (abs(diff.up) <= tol)
        alevel.up <- pchisq(2*(l.temp-l.temp.up),df=1,ncp=0,lower.tail = T)
        # print(c(delta.up, alevel.up, n,l.temp.up,k,iter1,iter2))
        if (!inside) {delta.up <- delta.up/((n+1)/n);iter2 <- iter2 + 1}
        if (is.nan(delta.up)) {param.est.temp.up[k] <- NA}
      } #for iter2

      #}
      #Using converged increment for parameter to get corresponding score and variance expressions;
      param.est.up <- insert(maximizer.temp1$estimate, ats=k, values=param.est.temp.up[k])
      l0.b.up = l.temp.up

      diff.up = l0.b.up - l.null
      converge <- (abs(diff.up) <= tol)
      if ((!converge| is.nan(l0.b.up)) & !is.nan(delta.up)) {iter1 <- iter1 + 1; n = n + 1} else {EXIT1 = T;}
      if (is.nan(delta.up)==T) {param.est.up[k] <- NA}
    } #for iter1

    upper.surv[is] <- param.est.up[k]


    ###############lower endpoint#####################
    n=ni
    #iter0 <- 1; converge = FALSE; l0.b.lo=0
    #while(l.temp > l.null & (!converge| is.nan(l0.b.lo)) & iter0<=30) {

    #assign initial values to parameter estimates
    param.est.lo <- maximizer0$estimate

    converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.lo=0; delta.lo=0
    while ((!converge| is.nan(l0.b.lo)) & iter1 <= 25 & !EXIT1 & !is.nan(delta.lo)) {

      # calculate log-lik, score and hessian under l0.b;
      maximizer.temp <-  nlm(
        f = loglik.mixture, p = param.est.lo, survt=survt, design.matrix=design.matrix,
        pl = pl, iterlim = 1, hessian=TRUE)
      score.temp = maximizer.temp$gradient
      hessian.temp = maximizer.temp$hessian
      if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06

      #### Approach 1: lambda = (2*(l0.b-l.null+e*A^-1*U)/(e*A^-1*e))^0.5
      # l0.b.lo <- -loglik.mixture(p=param.est.lo, survt, design.matrix,
      #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
      l0.b.lo <- -maximizer.temp$minimum
      inv.hessian.temp <- solve(hessian.temp)
      if (l0.b.lo < l.null + 0.5 * score.temp %*% inv.hessian.temp %*% score.temp) {
        lambda <- (-inv.hessian.temp %*% score.temp)[k]/inv.hessian.temp[k,k]
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.lo <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
      } else{
        lambda <- -(2*(l0.b.lo - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.lo <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)
        };

      # maximizing loop for unpenalized estimates;
      #if (pl == F) {
      inside <- FALSE; iter2 <- 1;
      while (!inside & iter2 <= 100 & !is.nan(delta.lo)) {

        # add increment to stepwise parameter value;
        param.est.temp.lo <- param.est.lo
        param.est.temp.lo[k] <- param.est.temp.lo[k] + delta.lo

        #compute loglikelihood function using lodated parameter values;
        maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.lo[-k], survt=survt,
                                param.est = param.est.temp.lo[k], k = k,
                                pl = pl, iterlim = iterlim, hessian=TRUE)
        l.temp.lo = -maximizer.temp1$minimum


        inside <- (l.temp.lo > l.null - 0.1)
        # diff.lo = l.temp.lo - l.null
        # converge0 <- (abs(diff.lo) <= tol)
        alevel.lo <- pchisq(2*(l.temp-l.temp.lo),df=1,ncp=0,lower.tail = T)
        # print(c(delta.lo, alevel.lo, n,l.temp.lo,k,iter1,iter2))
        if (!inside) {delta.lo <- delta.lo/((n+1)/n);iter2 <- iter2 + 1}
        if (is.nan(delta.lo)) {param.est.temp.lo[k] <- NA}

      } # for iter2;

      #}
      #Using converged increment for parameter to get corresponding score and variance expressions;
      param.est.lo <- insert(maximizer.temp1$estimate, ats=k, values=param.est.temp.lo[k])
      l0.b.lo = l.temp.lo

      diff.lo = l0.b.lo - l.null
      converge <- (abs(diff.lo) <= tol)
      if ((!converge | is.nan(l0.b.lo)) & !is.nan(delta.lo)) {iter1 <- iter1 + 1; n = n + 1} else {EXIT1 = T}
      if (is.nan(delta.lo)==T) {param.est.lo[k] <- NA}
    } #for iter1

    lower.surv[is] <- param.est.lo[k]
  }


  coef.table.cure <- cbind(
    'coef'        = maximizer0$estimate[index.cure.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.cure.v]),
    # 'LL.cure' = ll.cure,
    # 'LLR'         = llr.cure,
    # 'Pr(>chisq)'  = pval.cure,
    'LCI.95%' = lower.cure,
    'UCI.95%' = upper.cure
  );
  rownames(coef.table.cure) <- colnames(design.matrix);

  coef.table.surv <- cbind(
    'coef'        = maximizer0$estimate[index.surv.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.surv.v]),
    # 'LL.surv' = ll.surv,
    # 'LLR'         = llr.surv,
    # 'Pr(>chisq)'  = pval.surv,
    'LCI.95%' = lower.surv,
    'UCI.95%' = upper.surv
  );
  rownames(coef.table.surv) <- colnames(design.matrix);

  coef.table.alpha <- 0;


#run.time = proc.time() - init.time


#######################################
## Output tables from either method; ##
#######################################


out <- list(
  coefficients = list(
    cure = coef.table.cure,
    surv = coef.table.surv,
    alpha = coef.table.alpha
  )
);
class(out) <- c('mixcure', 'list');

return(out);

}


#### print.mixcure #############################################################
# DESCRIPTION
#   To print a mixcure object.
# INPUT
#   object : a mixcure object, which is an outcome of function mixcure.
#   digits : number of digits for printing, passed to print.default.
#   ...    : other parameters passed to print.default.
# OUTPUT
#   NULL.

print.mixcure <- function(object, digits = 3, ...) {
  sep.line.cure   <- paste(c(rep('-', 37),  ' CURE ' , rep('-', 37)), collapse = '');
  sep.line.surv   <- paste(c(rep('-', 37),  ' SURVIVAL ' , rep('-', 37)), collapse = '');
  sep.line.alpha <- paste(c(rep('-', 36), ' ALPHA ', rep('-', 36)), collapse = '');

  message(sep.line.cure);
  print.default(object$coefficients$cure,   digits = digits, ...);

  message(sep.line.surv);
  print.default(object$coefficients$surv,   digits = digits, ...);

  message(sep.line.alpha);
  print.default(object$coefficients$alpha, digits = digits, ...);

  return(NULL);
};


#### coef.mixcure ##############################################################
coef.mixcure <- function(object) {
  coefs <- c(
    object$coefficients$cure[, 'coef'],
    object$coefficients$surv[, 'coef'],
    object$coefficients$alpha[, 'coef']
  );
  names(coefs) <- c( paste('cure.', rownames(object$coefficients$cure)),
                     paste('surv.', rownames(object$coefficients$surv)),
                     rownames(object$coefficients$alpha) );

  return(coefs);
}


vcov.mixcure <- function(object) {
  return(object$cov);
}



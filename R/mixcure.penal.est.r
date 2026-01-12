################################
### Function for mixcure model##
#### penalized loglikelihoods ##
######################################################################
#### Last modified Jan 12 2026 for multiple variables as covariates ##
######################################################################

mixcure.penal.est <- function(formula, data, init, pl, iterlim = 200) { 
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
  
  #samp.s <- nrow(design.matrix)
  
  
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
    iterlim = iterlim, hessian=TRUE);

  hessmat <- maximizer0$hessian
  # if (det(hessmat) < 1e-05) 
  #   diag(hessmat) <- diag(hessmat) + 1e-06
  
var.mat <- solve(hessmat)

alpha.hat <- maximizer0$estimate[index.gamma];

loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value

######### Wald statistic inference ##########
  # confidence intervals for estimated coefficients
  z.score <- maximizer0$estimate / sqrt(diag(var.mat));
  
coef.table.cure <- cbind(
  'coef'        = maximizer0$estimate[index.cure.v],
  'exp(coef)'   = exp(maximizer0$estimate[index.cure.v]),
  'se(coef)'    = sqrt(diag(var.mat)[index.cure.v]),
  'z'           = z.score[index.cure.v],
  'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.cure.v]))),
  'LCI.95%' = maximizer0$estimate[index.cure.v] - 1.96 * sqrt(diag(var.mat)[index.cure.v]),
  'UCI.95%' = maximizer0$estimate[index.cure.v] + 1.96 * sqrt(diag(var.mat)[index.cure.v])
);
rownames(coef.table.cure) <- colnames(design.matrix);

coef.table.surv <- cbind(
  'coef'        = maximizer0$estimate[index.surv.v],
  'exp(coef)'   = exp(maximizer0$estimate[index.surv.v]),
  'se(coef)'    = sqrt(diag(var.mat)[index.surv.v]),
  'z'           = z.score[index.surv.v],
  'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.surv.v]))),
  'LCI.95%' = maximizer0$estimate[index.surv.v] - 1.96 * sqrt(diag(var.mat)[index.surv.v]),
  'UCI.95%' = maximizer0$estimate[index.surv.v] + 1.96 * sqrt(diag(var.mat)[index.surv.v])
);
rownames(coef.table.surv) <- colnames(design.matrix);

coef.table.alpha <- cbind(
  'coef'     = alpha.hat,
  'se(coef)' = sqrt(diag(var.mat)[index.gamma]),
  'z'        = z.score[index.gamma],
  'Pr(>|z|)' = 2 * (1 - pnorm(abs(z.score[index.gamma]))),
  'LCI.95%'  = maximizer0$estimate[index.gamma] - 1.96 * sqrt(diag(var.mat)[index.gamma]),
  'UCI.95%'  = maximizer0$estimate[index.gamma] + 1.96 * sqrt(diag(var.mat)[index.gamma]),
  'loglik' = -maximizer0$minimum 
);
rownames(coef.table.alpha) <- 'alpha';

#######################################
## Output tables from either method; ##
#######################################

colnames(var.mat) <- c(
  paste('cure.', colnames(design.matrix)), 
  paste('surv.', colnames(design.matrix)),
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



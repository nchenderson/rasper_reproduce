AngleRidge <- function(y, X, alpha, lambda, beta.prior) {
  ## The DistTL approach
  p <- ncol(X)
  XtX <- crossprod(X, X)
  diag(XtX) <- diag(XtX) + alpha
  Xty <- crossprod(X, y) + lambda*beta.prior
  beta.hat <- solve(XtX, Xty)

  fitted.vals <- as.numeric(X%*%beta.hat)
  return(list(coef=beta.hat, fitted.values=fitted.vals))
}

LooCVAngleRidge <- function(y, X, alpha.seq, lambda.seq, beta.prior) {

  nobs <- length(y)
  nlam <- length(lambda.seq)
  nalpha <- length(alpha.seq)
  loo_cv_vals <- alpha_vals <- lam_vals <- rep(NA, nalpha*nlam)
  count <- 1
  for(h in 1:nlam) {
    for(j in 1:nalpha) {
      mse.val <- rep(NA, nobs)
      for(k in 1:nobs) {
        y.train <- y[-k]
        X.train <- X[-k,]
        X.test <- matrix(X[k,], nrow=1, ncol=ncol(X))
        y.test <- y[k]

        ridge_fit <- AngleRidge(y=y.train, X=X.train, alpha=alpha.seq[j], lambda=lambda.seq[h], beta.prior=beta.prior)
        rfit <- as.numeric(X.test%*%ridge_fit$coef)
        mse.val[k] <- (y.test - rfit)*(y.test - rfit)
      }
      loo_cv_vals[count] <- mean(mse.val)
      alpha_vals[count] <- alpha.seq[j]
      lam_vals[count] <- lambda.seq[h]
      count <- count + 1
    }
  }
  loo_dframe <- data.frame(LooScore=loo_cv_vals, lambda=lam_vals, alpha=alpha_vals)
  return(loo_dframe)
}

InformativeRidge <- function(y, X, alpha, beta.prior) {
  ## The DistTL approach
  p <- ncol(X)
  XtX <- crossprod(X, X)
  diag(XtX) <- diag(XtX) + alpha
  Xty <- crossprod(X, y) + alpha*beta.prior
  beta.hat <- solve(XtX, Xty)

  fitted.vals <- as.numeric(X%*%beta.hat)
  return(list(coef=beta.hat, fitted.values=fitted.vals))
}

LooCVInformativeRidge <- function(y, X, alpha.seq, beta.prior) {

  nobs <- length(y)
  nalpha <- length(alpha.seq)
  loo_cv_vals <- rep(NA, nalpha)
  for(j in 1:nalpha) {
    mse.val <- rep(NA, nobs)
    for(k in 1:nobs) {
      y.train <- y[-k]
      X.train <- X[-k,]
      X.test <- matrix(X[k,], nrow=1, ncol=ncol(X))
      y.test <- y[k]

      ridge_fit <- InformativeRidge(y=y.train, X=X.train, alpha=alpha.seq[j], beta.prior=beta.prior)
      rfit <- as.numeric(X.test%*%ridge_fit$coef)

      mse.val[k] <- (y.test - rfit)*(y.test - rfit)
    }
    loo_cv_vals[j] <- mean(mse.val)
  }
  return(loo_cv_vals)
}

NonInformativeRidge <- function(y, X, alpha) {
  p <- ncol(X)
  XtX <- crossprod(X, X)
  Xty <- crossprod(X, y)
  diag(XtX) <- diag(XtX) + alpha
  beta.hat <- solve(XtX, Xty)

  fitted.vals <- as.numeric(X%*%beta.hat)
  return(list(coef=beta.hat, fitted.values=fitted.vals))
}


LooCVNoninformativeRidge <- function(y, X, alpha.seq) {

  nobs <- length(y)
  nalpha <- length(alpha.seq)
  loo_cv_vals <- rep(NA, nalpha)
  for(j in 1:nalpha) {
    mse.val <- rep(NA, nobs)
    for(k in 1:nobs) {
      y.train <- y[-k]
      X.train <- X[-k,]
      X.test <- matrix(X[k,], nrow=1, ncol=ncol(X))
      y.test <- y[k]

      ridge_fit <- NonInformativeRidge(y=y.train, X=X.train, alpha=alpha.seq[j])
      rfit <- as.numeric(X.test%*%ridge_fit$coef)

      mse.val[k] <- (y.test - rfit)*(y.test - rfit)
    }
    loo_cv_vals[j] <- mean(mse.val)
  }
  return(loo_cv_vals)
}

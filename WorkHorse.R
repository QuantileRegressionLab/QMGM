fit.cat = function(midFit, x, tau, offset, binary, lambda, rho, thresh) {
  n <- nrow(x)
  k <- length(midFit$yo)	
  up <- apply(tau - midFit$G, 1, function(x) which(x < 0)[1])
  FLAG <- sum(up == 1, na.rm = TRUE) + sum(is.na(up))
  low <- up - 1
  low[low == 0] <- 1
  low[is.na(low)] <- k
  up[is.na(up)] <- k
  Z <- cbind(midFit$yo[low], midFit$yo[up])
  PI <- t(apply(midFit$G, 1, function(x, p){
    up <- which(p - x < 0)[1]
    low <- up - 1
    low[low == 0] <- 1
    low[is.na(low)] <- length(x)
    up[is.na(up)] <- length(x)
    x[c(low, up)]
  }, p = tau))
  gamma <- (tau - PI[,1])/(PI[,2] - PI[,1])
  gamma[!is.finite(gamma)] <- 0
  B <- gamma*(Z[,2] - Z[,1]) + Z[,1]
  if(!is.null(lambda)){
    if(binary){
      B <- ao(B, lambda) - offset
    } else {
      B <- bc(B, lambda) - offset
    }
  } else {B <- B - offset}
  
  # is tau outside range?
  r <- attr(midFit$G, "range")
  if(any(tau < r[1]) | any(tau > r[2])){
    warning("tau = ", tau, " is outside mid-probabilities range ", "[", round(r[1], 3), ", " , round(r[2], 3), "] for ", FLAG, " out of ", n, " observations. See details for ?midrq")
  }
  
  # penalized LASSO
  fit.par <- as.numeric(coef(glmnet(x = x[,-1], y = B, family = "gaussian", lambda = rho, alpha = 1, 
                                     thresh = thresh))) #  penalty.factor = apply(x[,-1], 2, sd),
  
  return(list(beta = fit.par, B = B))
}



midrqLoss = function(b, xx, midFit, tau, rho) {
  
  xb = cbind(1, xx)%*%b
  condist = npcdist(bws = midFit, exdat = xx, eydat = xb)$condist
  out = mean((condist - tau)^2) + rho * sum(abs(b[-1]))

  return(out)
}

midCDF_est = function(Obs) {
  
  n = nrow(Obs)
  p = ncol(Obs)
  Obs.names = names(Obs)

  # Estimate mid-CDF
  midFit = foreach(j = 1:p, .packages = c("quantreg", "Qtools", "np")) %dopar% {
    ff = as.formula(paste0(Obs.names[j], " ~ ", paste(Obs.names[-j], collapse = " + ")))
    y = Obs[,j]
    x = model.matrix(ff, Obs)
    # is the response binary?
    binary = setequal(sort(unique(y)), c(0,1))
    # is the response categorical?
    # categorical = length(unique(y)) == n
    categorical = T
    
    if(categorical) {
      cmidecdf.fit(x = x, y = y, ecdf_est = "logit")
    } else {
      list(dist = npcdistbw(xdat = x[,-1], ydat = y)) #, cykertype = "epanechnikov" ftol = .1, tol = .1)
           # dens = npcdensbw(xdat = x[,-1], ydat = y, cykertype = "epanechnikov", nmulti = 2, method = "cv.ls"))
    }
  }
  
  return(midFit)
}


QMGM = function(Obs, tau, midFit, rho, offset.var, type.var, thresh) {
  
  n = nrow(Obs)
  p = ncol(Obs)
  Obs.names = names(Obs)
  if(missing(offset.var)) offset.var = rep(0, n)
  
  # fit model
  fit = list()
  for (j in 1:p) {
    ff = as.formula(paste0(Obs.names[j], " ~ ", paste(Obs.names[-j], collapse = " + ")))
    y = Obs[,j]
    x = model.matrix(ff, Obs)
    # is the response continuous?
    if(type.var[j] == "g") lambda = NULL else lambda = 0
    # check offset
    offset = offset.var[,j]
    # is the response binary?
    binary = setequal(sort(unique(y)), c(0, 1))
    # is the response categorical?
    # categorical = length(unique(y)) == n
    categorical = T
    
    if(categorical) {
      fit[[j]] = fit.cat(midFit = midFit[[j]], x = x, tau = tau, offset = offset, binary = binary, lambda = lambda, rho = rho, thresh = thresh)
      yhat = x%*%fit[[j]]$beta + offset
      if(!is.null(lambda)) {
        if(binary) {
          fit[[j]]$yhat = apply(yhat, 2, invao, lambda = lambda)
        } else {
          fit[[j]]$yhat = apply(yhat, 2, invbc, lambda = lambda)
        }
      } else {
        fit[[j]]$yhat = yhat
      }
    } else {
      b0 = as.numeric(coef(glmnet(x = x[,-1], y = y, alpha = 1, lambda = rho)))
      q0 = npqreg(bws = midFit[[j]]$dist, tau = tau, itmax = 1e3)$quantile
      fit[[j]] = list(beta = as.numeric(coef(glmnet(x = x[,-1], y = q0, alpha = 1, lambda = rho))))
      # fit[[j]] = proxGD(x = x[,-1], y = y, tau = tau, b0 = b0, bw = midFit[[j]]$dist, bwd = midFit[[j]]$dens, 
      #              lambda = rho, alpha = 100, iter = 1e2, conv = 1e-04)
      # fit[[j]] = optim(par = b0, fn = midrqLoss, midFit = midFit[[j]], x = x[,-1], tau = tau, rho = rho, 
      #                  method = "BFGS", control = list(abstol = 1e-06, reltol = 1e-06))
      # fit[[j]]$par = ifelse(abs(fit[[j]]$par) < thresh, 0, fit[[j]]$par)
      fit[[j]]$yhat = x%*%fit[[j]]$beta
    }
  }
  
  betahat = lapply(1:p, function(j) fit[[j]]$beta)
  yhat = sapply(1:p, function(j) fit[[j]]$yhat)
  names(betahat) = names(yhat) = Obs.names
  score = -sum(.5 * (Obs - yhat)^2)
  # llscore = sum(log(colSums((Obs - yhat)^2)))
  gAIC = -2 * score + 2 * sum(unlist(betahat) != 0)
  gBIC = -2 * score + log(n) * sum(unlist(betahat) != 0)
  # gAIC = n * llscore + 2 * sum(unlist(betahat) != 0)
  # gBIC = n * llscore + log(n) * sum(unlist(betahat) != 0)
  gcrit = c(gAIC, gBIC)
  
  # Adjacency matrix
  A = matrix(NA, p, p)
  for(i in 1:p) {
    A[i,] = as.numeric(append(x = betahat[[i]][-1], values = 0, i - 1) != 0)
  }
  # A = 1 * ((A + t(A)) != 0) # OR
  A = 1 * ((A + t(A)) == 2) # AND
  
  return(list(rho = rho, betahat = betahat, yhat = yhat, Adj.mat = A, Fitted = yhat, score = score, gcrit = gcrit))
}


# EBIC
extendedBIC <- function(gamma,omegahat,S,n) {
    p = nrow(omegahat)
    es = sum(omegahat[upper.tri(omegahat)]!=0)
    return(-log(det(omegahat)) + sum(diag(omegahat %*% S)) + es * (log(n)/n) + es * gamma * (4*log(p)/n))
  }


# Graph performance
Graph.performance = function(est.A, true.adj) {
  
  TP = sum(true.adj[upper.tri(true.adj)] != 0 & est.A[upper.tri(est.A)] != 0) 
  TN = sum(true.adj[upper.tri(true.adj)] == 0 & est.A[upper.tri(est.A)] == 0) 
  FP = sum(true.adj[upper.tri(true.adj)] == 0 & est.A[upper.tri(est.A)] != 0)
  FN = sum(true.adj[upper.tri(true.adj)] != 0 & est.A[upper.tri(est.A)] == 0)
  
  out = c(TP / (TP + FP), 
          TP / (TP + FN), 
          2 * TP / (2 * TP + FP + FN), 
          (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
          (TP + TN) / (TP + TN + FP + FN))
  names(out) = c("Precision", "Recall", "F1-score", "Matthews CC", "Accuracy")

  return(out)
}

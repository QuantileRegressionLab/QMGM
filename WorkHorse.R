check <- function(x, tau = .5){
  x * (tau - (x<0))
}

cmidecdf.fit.pen = function(x, y, intercept, ecdf_est, npc_args = list(), theta = NULL){
  p <- ncol(x)
  n <- length(y)
  if(length(unique(y)) == n) {
    x.er <- extendrange(y, f = 1)
    # yo = unique(round(sort(y), digits = 1))
    yo <- unique(na.omit(sort(c(seq(x.er[1], x.er[2], length = 2 * round(sqrt(n) * 3.9)),
                                quantile(y, c(0.01,0.05,0.25,0.50,0.75,0.95,0.99), na.rm = TRUE)))))
  } else {
    yo <- sort(unique(y))
    # x.er <- extendrange(yo, f = c(0, .1))
    # yo <- round(sort(unique(seq(0, x.er[2], length = round(2 * sqrt(n) * 3.9)))))
  }
  K <- length(yo)
  Z <- mapply(function(t, y) (y <= t), (yo), MoreArgs = list(y = (y)))
  if(missing(intercept) & all(x[,1] == 1)) intercept <- TRUE
  
  
  if(ecdf_est %in% c("logit", "probit", "cloglog")){
    fitbin <- apply(Z, 2, function(z, x, link) suppressWarnings(glm(z ~ x - 1, epsilon = 1e-10, maxit = 1e3, family = binomial(link))), x = x, link = ecdf_est)
    Fhat <- sapply(fitbin, predict, type = "response")
    # Fse <- sapply(fitbin, function(x) predict(x, type = "response", se.fit = TRUE)$se.fit)
    
    # fitbin <- apply(Z, 2, function(z, x, link) suppressWarnings(cv.glmnet(x = x, y = z, family = binomial(link),
    #                 lambda = lambda)), x = x[,-1], link = ecdf_est)
    # Fhat <- sapply(fitbin, predict, type = "response", newx = x[,-1])
    
    bhat <- sapply(fitbin, coef) # p x K
    linkinv <- family(fitbin[[1]])$linkinv
  }
  
  
  # rearrange if CDF is not monotone
  for(j in 1:n){
    tmp <- Fhat[j,]
    if(any(diff(tmp) < 0)){
      sf <- rearrange(stepfun(yo, c(-Inf, tmp)))
      Fhat[j,] <- sf(yo)
    }
  }
  
  M <- apply(Fhat, 1, diff)
  if(ncol(Fhat) > 2) M <- t(M)
  G <- Fhat[,-1] - 0.5*M
  G <- cbind(Fhat[,1]/2, G)
  if(length(unique(y)) == n) G <- Fhat
  r <- c(max(G[,1]), min(G[,ncol(G)]))
  
  attr(G, "range") <- r
  ecdf_fit <- if(ecdf_est == "npc") bw else list(coef = bhat, linkinv = linkinv)
  
  ans <- list(G = G, Fhat = Fhat, yo = yo, ecdf_fit = ecdf_fit, ecdf_est = ecdf_est) # Fse = Fse
  class(ans) <- "cmidecdf"
  
  return(ans)
}

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
  if (length(unique(B)) == 1) {
    fit.par <- as.numeric(c(B[1], rep(0, ncol(x) - 1)))
    # fit.par <- as.numeric(coef(glmnet(x = x[,-1], y = jitter(B, factor = 1e-06), family = "gaussian", lambda = rho * 1, alpha = 1,
    #                                   thresh = thresh)))
  } else {
    fit.par <- as.numeric(coef(glmnet(x = x[,-1], y = B, family = "gaussian", lambda = rho * 1, alpha = 1,
                                     thresh = thresh))) # penalty.factor = apply(x[,-1], 2, sd), , standardize.response = F, standardize = T
  }
  
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
  midFit = foreach(j = 1:p, .packages = c("quantreg", "Qtools")) %dopar% {
    ff = as.formula(paste0(Obs.names[j], " ~ ", paste(Obs.names[-j], collapse = " + ")))
    y = Obs[,j]
    x = model.matrix(ff, Obs)
    # x.SIR = cbind(1, components(SIRasymp(X = x[,-1], h = 4, y = y, k = 3), which = "k"))
    
    cmidecdf.fit.pen(x = x, y = y, ecdf_est = "logit") #, cykertype = "epanechnikov" ftol = .1, tol = .1)
  }
  
  return(midFit)
}


QMGM = function(Obs, tau, midFit, rho, offset.var, type.var, ruleReg, thresh) {
  
  n = nrow(Obs)
  p = ncol(Obs)
  Obs.names = names(Obs)
  if(missing(offset.var)) offset.var = matrix(0, n, p)
  

  # fit model
  fit = list()
  for (j in 1:p) {
    ff = as.formula(paste0(Obs.names[j], " ~ ", paste(Obs.names[-j], collapse = " + ")))
    y = Obs[,j]
    x = model.matrix(ff, Obs)
    # x = as.data.frame(x)
    # for(j in which(1*apply(x, 2, function(r) length(table(r)) == 2) == 1)) {
    #   x[,j] = factor(x[,j])
    # }

    # is the response continuous?
    if(type.var[j] == "g") lambda = NULL else lambda = 0
    # check offset
    offset = offset.var[,j]
    # is the response binary?
    binary = setequal(sort(unique(y)), c(0, 1))
    # binary = ifelse(length(table(y)) == 2, T, F)
    
    # response on the scale of linear predictor
    if(!is.null(lambda)){
      if(binary){
        hy = ao(y, lambda)
      } else {
        if(any(y <= 0)) stop("The response must be strictly positive")
        hy = bc(y, lambda)
      }
    } else {hy = y}
    
    fit[[j]] = fit.cat(midFit = midFit[[j]], x = x, tau = tau, offset = offset, binary = binary, lambda = lambda, rho = rho, thresh = thresh)
    fit[[j]]$hy = hy
    yhat = x%*%fit[[j]]$beta + offset
    fit[[j]]$yhat = yhat
    if(!is.null(lambda)) {
      if(binary) {
        fit[[j]]$fitted.values = apply(yhat, 2, invao, lambda = lambda)
      } else {
        fit[[j]]$fitted.values = apply(yhat, 2, invbc, lambda = lambda)
      }
    } else {
      fit[[j]]$fitted.values = yhat
    }
  }
  
  betahat = lapply(1:p, function(j) fit[[j]]$beta)
  Hy = sapply(1:p, function(j) fit[[j]]$hy)
  Yhat = sapply(1:p, function(j) fit[[j]]$yhat)
  Fitted = sapply(1:p, function(j) fit[[j]]$fitted.values)
  Residuals = Obs - Fitted
  Std.Residuals = apply(Residuals, 2, scale, center = F)
  names(betahat) = names(Fitted) = names(Hy) = names(Yhat) = names(Residuals) = Obs.names
  score = sum(log(colSums(check(Residuals, tau = tau))))
  gAIC = score + 2 * 1 * (sum(unlist(betahat) != 0) - 0) / (2 * n) # log(p - 1)
  gBIC = score + log(n) * 1 * (sum(unlist(betahat) != 0) - 0) / (2 * n)
  gBICp = score + log(n) * log(p - 1) * (sum(unlist(betahat) != 0) - 0) / (2 * n)
  gBIC2p = score + log(n) * log(p - 1)/2 * (sum(unlist(betahat) != 0) - 0) / (2 * n)
  gBIC3p = score + log(n) * log(p - 1)/3 * (sum(unlist(betahat) != 0) - 0) / (2 * n)
  
  gEBIC0 = score + (log(n) + 4 * 0 * log(p - 1)) * (sum(unlist(betahat) != 0) - 0) / (2 * n) 
  gEBIC5 = score + (log(n) + 4 * 0.5 * log(p - 1)) * (sum(unlist(betahat) != 0) - 0) / (2 * n) 
  gEBIC1 = score + (log(n) + 4 * 1 * log(p - 1)) * (sum(unlist(betahat) != 0) - 0) / (2 * n) 
  gcrit = c(gAIC, gBIC, gBICp, gBIC2p, gBIC3p, gEBIC0, gEBIC5, gEBIC1)
  
  # Adjacency matrix (weighted and signs matrices)
  A = wA = signsA = matrix(NA, p, p)
  for(i in 1:p) {
    A[i,] = as.numeric(append(x = betahat[[i]][-1], values = 0, i - 1) != 0)
    wA[i,] = abs(append(x = betahat[[i]][-1], values = 0, i - 1))
    signsA[i,] = sign(append(x = betahat[[i]][-1], values = 0, i - 1))
  }
  wA = 0.5 * (wA + t(wA))
  signsA = signsA + t(signsA)
  
  if(ruleReg == "OR") {
    A = 1 * ((A + t(A)) != 0) # OR
  } else if (ruleReg == "AND") {
    A = 1 * ((A + t(A)) == 2) # AND
  }
  wA = wA * A
  signsA = signsA * A
  signsA[A == 0] = NA
  sign_colors = matrix("darkgrey", p, p)
  sign_colors[signsA == 1] = sign_colors[signsA == 2] = "darkgreen"
  sign_colors[signsA == -1] = sign_colors[signsA == -2] = "red"
  
  # save in output
  out_obj = list()
  out_obj$tau = tau
  # out_obj$x = x
  # out_obj$y = y
  out_obj$offset = offset
  out_obj$lambda = lambda
  out_obj$binary = binary
  out_obj$intercept
  out_obj$type.var = type.var
  out_obj$rho = rho
  out_obj$betahat = betahat
  out_obj$Hy = Hy
  out_obj$Yhat = Yhat
  # out_obj$Fitted = Fitted
  # out_obj$residuals = Residuals
  out_obj$Adj = A
  out_obj$wAdj = wA
  out_obj$signsAdj = signsA
  out_obj$edgecolorAdj = sign_colors
  out_obj$gcrit = gcrit
 
  return(out_obj)
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

# Graph summary
Graph.summary = function(est.graph) {
  
  p = length(V(est.graph))
  centrality <- data.frame(row.names   = V(est.graph)$name,
                           Degree      = degree(graph = est.graph, mode = "all"),
                           Closeness   = closeness(graph = est.graph, mode = "all"),
                           Betweenness = betweenness(graph = est.graph, directed = F),
                           Eigenvector = eigen_centrality(graph = est.graph)$vector,
                           Eccentricity = eccentricity(graph = est.graph, mode = "all"))
  
  out_obj = list()
  out_obj$centrality = centrality
  out_obj$n.edges = 0.5 * sum(as_adjacency_matrix(est.graph))
  out_obj$fr.edges = 0.5 * sum(as_adjacency_matrix(est.graph)) / choose(p, 2)
  
  return(out_obj)
}

# Compute the likelihood for MGM
calcLL <- function(X,
                   y,
                   beta_vector,
                   type,
                   level,
                   v,
                   weights,
                   LLtype = 'model')
  
  
{
  
  if(missing(level)) stop('No levels passed to calcLL !')
  
  # This function calculates three different LL:
  # 1) LLtype = 'model': The LL of a given model via fit
  # 2) LLtype = 'nullmodel': The LL of the Null (Intercept) Model
  # 3) LLtype = 'saturated': The LL of the saturated model
  
  X <- as.matrix(X)
  y <- as.vector(y)
  
  n <- nrow(X)
  if(missing(weights)) weights <- rep(1, n)
  
  if(LLtype == 'model') {
    
    if(type[v] == 'g') {
      beta_vector <- matrix(beta_vector, ncol = 1)
      predicted_mean <- cbind(rep(1, n), X) %*% as.vector(beta_vector)
      sd_residual <- sd(y-predicted_mean)
      LL_model <- dnorm(y, mean = predicted_mean, sd = sd_residual, log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    
    if(type[v] == 'p') {
      beta_vector <- matrix(beta_vector, ncol = 1)
      predicted_mean <- cbind(rep(1, n), X) %*% as.vector(beta_vector)
      LL_model <- dpois(y, exp(predicted_mean), log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    
    if(type[v] == 'c') {
      
      n_cats <- level[v] # number of levels
      
      ## Compute LL (see http://www.stanford.edu/~hastie/Papers/glmnet.pdf, equation 22)
      m_respdum <- matrix(NA, n, n_cats) # dummy for data
      m_coefs <- matrix(NA, n, n_cats) # dummy for coefficients
      cats <- unique(y)
      
      LL_n <- rep(NA, n) # Storage
      m_LL_parts <- matrix(NA, nrow = n, ncol=n_cats+1)
      
      for(catIter in 1:n_cats) {
        m_respdum[,catIter] <- (y==cats[catIter])*1 # dummy matrix for categories
        m_coefs[,catIter] <- cbind(rep(1, n), X) %*% as.vector(beta_vector[[catIter]])
        m_LL_parts[,catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
      }
      
      m_LL_parts[, n_cats+1] <- - log(rowSums(exp(m_coefs))) # the log part, see eq (22)
      LL_n <- rowSums(m_LL_parts) # sum up n_cat + 1 parts
      mean_LL_model <- sum(LL_n * weights) # apply weighting
      
    }
    
  }
  
  return(mean_LL_model)
  
}

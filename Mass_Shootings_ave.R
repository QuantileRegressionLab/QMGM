rm(list = ls())
gc()
graphics.off()
library(readxl)
library(multiUS)
library(stringr)

Offender = read_excel("Mass_Shootings.xlsx", sheet = "Offender", skip = 1)
Crime = read_excel("Mass_Shootings.xlsx", sheet = "Crime", skip = 1)
Traumas = read_excel("Mass_Shootings.xlsx", sheet = "Traumas", skip = 1)
Crisis = read_excel("Mass_Shootings.xlsx", sheet = "Crisis", skip = 1)
Mental = read_excel("Mass_Shootings.xlsx", sheet = "Mental", skip = 1)
# Mental = Mental[,-9]
Social = read_excel("Mass_Shootings.xlsx", sheet = "Social", skip = 1)
Motivation = read_excel("Mass_Shootings.xlsx", sheet = "Motivation", skip = 1)
Victims = read_excel("Mass_Shootings.xlsx", sheet = "Victims", skip = 0, col_types = rep("numeric", 2))

Crime = rowMeans(Crime != "0")
Traumas = rowMeans(Traumas != "0")
Crisis = rowMeans(Crisis != "0")
Social = rowMeans(Social != "0")
Mental = rowMeans(Mental != "0")
Motivation = rowMeans(Motivation != "0")

VictimsAge = aggregate(Victims$Age, by = list(Victims$`Case #`), FUN = median, na.rm = T)
for(j in 1:max(Victims$`Case #`)) {
  if(any(is.na(Victims$Age[Victims$`Case #` == j]))) {
    Victims$Age[Victims$`Case #` == j][is.na(Victims$Age[Victims$`Case #` == j])] = VictimsAge[j,2]
  }
}
VictimsAge = VictimsAge[-c(146),2]

Offender$Race[Offender$Race != 0] = 1
Offender$`Military Service`[Offender$`Military Service` == 2] = 0
Offender$`Relationship Status`[Offender$`Relationship Status` == 2] = 1
Offender$`Relationship Status`[Offender$`Relationship Status` == 3] = 0

DATA = data.frame(Offender$`Number Killed`, Offender$`Number Injured`, Offender$`Total Firearms Brought to the Scene`, 
                  Offender$Age, VictimsAge, Offender$`Insider / Outsider`, Offender$Immigrant, Offender$`Relationship Status`,
                  Social, Crime, Traumas, Crisis, Mental, Motivation, row.names = Offender$`Case #`, check.names = T)

DATA = DATA[-145,] # Wilkinsburg shooting (March 9, 2016)
DATA = head(DATA, -5) # 1966-2022 (2023 is incomplete)
100 *  sum(is.na(DATA)) / prod(dim(DATA)) # % of missing data

names(DATA) = as.character(c("Killed", "Injured", "Firearms.brought.to.the.scene", 
                "Age", "Victims.age", "Insider", "Immigrant", "Relationship.status",
                "Social", "Crime", "Traumas", "Crisis", "Mental", "Motivation"))


DATA = KNNimp(data = DATA, k = round(sqrt(nrow(na.omit(DATA)))), meth = "median")
apply(DATA[,1:8], 2, table)
# DATA %>% mutate_at(vars(names(DATA)[1:8]), list(~ round(., 0)))
# DATA = DATA[DATA$.Number.Injured. < 866,]

n = nrow(DATA)
p = ncol(DATA)
pc = ncol(DATA)


# QMGM --------------------------------------------------------------------
library(qgraph)
library(igraph)
library(mgm)
library(foreach)
library(parallel)
library(doParallel)
source("WorkHorse.R")

ncores = detectCores()
ncores
cl = makeCluster(ncores)
registerDoParallel(cl)
invisible(clusterEvalQ(cl = cl, source("WorkHorse.R")))


tau = 0.5
# tau = c(0.2, 0.5, 0.8)
# tau = seq(1/8, 7/8, by = 1/8)
rhol = exp(seq(log(1e-03), log(5), length.out = 1e2))  # seq(0, 1, length.out = 50)
thresh = 1e-07
# type.var = c(rep("p", 2), rep("g", 2), rep("p", 10))
type.var = c(rep("p", 3), rep("g", 2), rep("p", 3), rep("g", 6))
vertex.shape = ifelse(type.var == "g", "circle", "square")

B = 1
resample.QMGM = list()
for(b in 1:B) {
  print(b)
  set.seed(b)
  boot.id = sample(x = 1:n, size = n, replace = F)
  DATA_boot = DATA[boot.id,]
  DATA_boot[,type.var == "g"] = apply(DATA_boot[,type.var == "g"], 2, scale)
  DATA_boot$Injured = DATA_boot$Injured + 0.1
  # DATA_boot$Social = DATA_boot$Social + 0.1
  # DATA_boot$Traumas = DATA_boot$Traumas + 0.1
  # DATA_boot$Crisis = DATA_boot$Crisis + 0.1
  # DATA_boot$Mental = DATA_boot$Mental + 0.1
  # DATA_boot$Motivation = DATA_boot$Motivation + 0.1
  # DATA_boot$Crime = DATA_boot$Crime + 0.1
  
  midFit = midCDF_est(Obs = DATA_boot)
  output = foreach(t = 1:length(tau), .packages = c("foreach", "quantreg", "Qtools", "np", "glmnet")) %do% {
    foreach(r = 1:length(rhol), .packages = c("foreach", "quantreg", "Qtools", "np", "glmnet")) %dopar% {
      QMGM(Obs = DATA_boot, tau = tau[t], midFit = midFit, rho = rhol[r], type.var = type.var, ruleReg = "OR", thresh = thresh)
    }
  }
  resample.QMGM[[b]] = output
}

wA.opt.array = edgecolorA.opt.array = array(NA, dim = c(p, p, B))
for(b in 1:B) {
  output = resample.QMGM[[b]]
  est.A = lapply(1:length(tau), function(t) lapply(1:length(rhol), function(r) output[[t]][[r]]$Adj))
  
  est.A.grid = lapply(1:length(rhol), function(r) list())
  for(r in 1:length(rhol)) {
    .list = lapply(1:length(tau), function(t) est.A[[t]][[r]])
    est.A.grid[[r]] = 1 * (Reduce('+', .list) != 0) # Reduce('*', .list)
  }
  
  # graph selection
  crit.mat = sapply(1:length(tau), function(t) sapply(1:length(rhol), function(r) output[[t]][[r]]$gcrit), simplify = "array")
  crit.mat.grid = apply(crit.mat, c(1, 2), sum)
  rho.opt.index = apply(crit.mat.grid, 1, which.min)
  
  A.opt = est.A.grid[[rho.opt.index[2]]]
  wA.opt = A.opt * (Reduce("+", lapply(1:length(tau), function(t) output[[t]][[rho.opt.index[[2]]]]$wAdj)) / length(tau))
  wA.opt.array[,,b] = wA.opt
  
  edgecolor.array = array(unlist(lapply(1:length(tau), function(t) output[[t]][[rho.opt.index[[2]]]]$edgecolorAdj)), dim = c(p, p, length(tau)))
  edgecolorA.opt = matrix("darkgrey", p, p) 
  edgecolorA.opt[apply(edgecolor.array, 1:2, function(x) (!any(x == "red")) * any(x == "darkgreen")) == 1] = "darkgreen"
  edgecolorA.opt[apply(edgecolor.array, 1:2, function(x) (any(x == "red")) * (!any(x == "darkgreen"))) == 1] = "red"
  edgecolorA.opt.array[,,b] = edgecolorA.opt
  
  # sign.list = Reduce("+", lapply(1:length(tau), function(t) output[[t]][[rho.opt.index[[2]]]]$signsAdj))
  # edgecolorA.opt[sign.list == (1 * length(tau))] = edgecolorA.opt[sign.list == (2 * length(tau))] = "darkgreen"
  # edgecolorA.opt[sign.list == -(1 * length(tau))] = edgecolorA.opt[sign.list == -(2 * length(tau))] = "red"
  # edgecolorA.opt.array[,,b] = edgecolorA.opt
}

100*apply(wA.opt.array != 0, 1:2, mean)
apply(wA.opt.array, 1:2, mean)

wA.opt.ave = (100 * apply(wA.opt.array != 0, 1:2, mean) >= 0) * apply(wA.opt.array, 1:2, mean)
sum((wA.opt.ave != 0) / 2) / choose(p, 2)
wA.opt.ave[abs(wA.opt.ave) <= 1e-04] = 0

groupsV = list("Shooting" = 1:3, "Characteristics" = c(4,5,7,8), "Background" = c(6,9:p))
# ("Victims" = 1:, "Background" = 3:7, "Weapons" = 8,
#                "Social" = 9, "Crime" = 10, "Trauma" = 11, "Crisis" = 12, "Mental" = 13, "Motivation" = p)
namesV = gsub(".", " ", names(DATA), fixed = T)


jpeg("Mass_Shootings_ave_72.jpg", height=2*1000, width=2*1500, unit='px')
qgraph(wA.opt.ave,
       layout = "spring", repulsion = 1,
       edge.color = apply(edgecolorA.opt.array, 1:2, FUN = function(x) names(sort(-table(x)))[1]),
       nodeNames = namesV, 
       shape = ifelse(type.var == "g", "circle", "square"),
       color = c("#ED3939", "#53B0CF", "#FFB026"), # #72CF53
       # color = c("lightgreen", "violet", "yellow", "orange", "green", "#E69F00", "lightblue", "red", "steelblue"),
       groups = groupsV,
       legend.mode="style2", legend.cex= 2.25,
       vsize = 5, esize = 15)
dev.off()

est.graph.ave = graph_from_adjacency_matrix(adjmatrix = 1*(wA.opt.ave != 0), mode = "undirected")
A.opt.summary.ave = Graph.summary(est.graph = est.graph.ave)


# MGM ---------------------------------------------------------------------
type.var.mgm = type.var
type.var.mgm[(1:pc)[which(apply(DATA, 2, function(x) length(unique(x))) == 2)]] = "c"
level.mgm = ifelse(type.var.mgm == "c", 2, 1)

B = 1
wA.opt.mgm.array = edgecolorA.opt.mgm.array = array(NA, dim = c(p, p, B))
for(b in 1:B) {
  print(b)
  set.seed(b)
  boot.id = sample(x = 1:n, size = n, replace = F)
  DATA_mgm_boot = DATA[boot.id,]
  DATA_mgm_boot[,type.var.mgm == "g"] = apply(DATA_mgm_boot[,type.var.mgm == "g"], 2, scale)
  DATA_mgm_boot = as.matrix(DATA_mgm_boot)
  
  mgm.output = foreach(r = 1:length(rhol), .packages = c("mgm")) %dopar% {
    mgm(data = DATA_mgm_boot, type = type.var.mgm, lambdaSeq = rhol[r], lambdaSel = "EBIC",
        threshold = "none", thresholdCat = F, lambdaGam = 0, ruleReg = "OR", pbar = F, signInfo = F, k = 2)
  }
  
  est.A.mgm = lapply(1:length(rhol), function(r) mgm.output[[r]]$pairwise$wadj)
  
  LL_models = mgm.aic = mgm.bic = mgm.bicp = mgm.bic2p = mgm.bic3p = matrix(NA, p, length(rhol))
  for(j in 1:p) {
    for(r in 1:length(rhol)) {
      LL_models[j,r] = calcLL(X = DATA_mgm_boot[,-j], y = DATA_mgm_boot[,j], beta_vector = mgm.output[[r]]$nodemodels[[j]]$model, 
                              type = type.var.mgm, v = j, level = level.mgm)
      
      if(type.var.mgm[j] == "c") {
        beta_vector = sum(mgm.output[[r]]$nodemodels[[j]]$model[[1]] != 0)
      } else {
        beta_vector = sum(mgm.output[[r]]$nodemodels[[j]]$model != 0)
      }
      mgm.aic[j,r] = -2*LL_models[j,r] + 2*1*beta_vector # log(p-1)
      mgm.bic[j,r] = -2*LL_models[j,r] + log(n)*1*beta_vector
      mgm.bicp[j,r] = -2*LL_models[j,r] + log(n)*log(p-1)*beta_vector
      mgm.bic2p[j,r] = -2*LL_models[j,r] + log(n)*log(p-1)/2*beta_vector
      mgm.bic3p[j,r] = -2*LL_models[j,r] + log(n)*log(p-1)/3*beta_vector
    }
  }
  
  rho.opt.mgm.index = matrix(NA, 5, 1)
  rho.opt.mgm.index[1] = which.min(colSums(mgm.aic))
  rho.opt.mgm.index[2] = which.min(colSums(mgm.bic))
  rho.opt.mgm.index[3] = which.min(colSums(mgm.bicp))
  rho.opt.mgm.index[4] = which.min(colSums(mgm.bic2p))
  rho.opt.mgm.index[5] = which.min(colSums(mgm.bic3p))
  
  mgm.adj = est.A.mgm[[rho.opt.mgm.index[2]]]
  wA.opt.mgm.array[,,b] = mgm.adj
  edgecolorA.opt.mgm.array[,,b] = mgm.output[[rho.opt.mgm.index[2]]]$pairwise$edgecolor
  # mgm.ebic = sum(sapply(1:p, function(j) mgm.output$nodemodels[[j]]$EBIC))
}

100*apply(wA.opt.mgm.array != 0, 1:2, mean)
apply(wA.opt.mgm.array, 1:2, mean)

wA.opt.mgm.ave = (100 * apply(wA.opt.mgm.array != 0, 1:2, mean) >= 0) * apply(wA.opt.mgm.array, 1:2, mean)
sum((wA.opt.mgm.ave != 0) / 2) / choose(p, 2)
wA.opt.mgm.ave[abs(wA.opt.mgm.ave) <= 1e-04] = 0

jpeg("Mass_Shootings_ave_mgm2.jpg", height=2*1000, width=2*1500, unit='px')
qgraph(wA.opt.mgm.ave, 
       layout = "spring", repulsion = 1,
       edge.color = apply(edgecolorA.opt.mgm.array, 1:2, FUN = function(x) names(sort(-table(x)))[1]),
       nodeNames = namesV,
       shape = ifelse(type.var.mgm == "g", "circle", "square"),
       color = c("#ED3939", "#53B0CF", "#FFB026"), # #72CF53
       # color = c("lightgreen", "violet", "yellow", "orange", "green", "#E69F00", "lightblue", "red", "steelblue"),
       groups = groupsV,
       legend.mode="style2", legend.cex= 2.25,
       vsize = 5, esize = 15)
dev.off()

est.graph.mgm.ave = graph_from_adjacency_matrix(adjmatrix = 1*(wA.opt.mgm.ave != 0), mode = "undirected")
A.opt.summary.mgm.ave = Graph.summary(est.graph = est.graph.mgm.ave)

##########################
# compare the two graphs
# global centrality measures
centralization.degree(est.graph.ave, normalized = T)$centralization
centralization.betweenness(est.graph.ave, normalized = T)$centralization
centralization.closeness(est.graph.ave, normalized = T)$centralization

centralization.degree(est.graph.mgm.ave, normalized = T)$centralization
centralization.betweenness(est.graph.mgm.ave, normalized = T)$centralization
centralization.closeness(est.graph.mgm.ave, normalized = T)$centralization

# intersection
graph.intersection(est.graph.ave, est.graph.mgm.ave)
length(E(graph.intersection(est.graph.ave,est.graph.mgm.ave))) # number of common edges

# difference
difference(est.graph.ave, est.graph.mgm.ave)
graph.difference(est.graph.ave, est.graph.mgm.ave)

# Hamming distance
(sum(get.adjacency(est.graph.ave) != get.adjacency(est.graph.mgm.ave)) / 2) / choose(p, 2)
ecount(est.graph.ave) + ecount(est.graph.mgm.ave) - 2 * ecount(graph.intersection(est.graph.ave, est.graph.mgm.ave))

library(ggplot2)
g = centralityPlot(list(MQGM7 = est.graph.ave, MGM = est.graph.mgm.ave), labels = namesV,
               include = c("Degree","Betweenness","Closeness"))
g = g + geom_point(size = 5) + theme(text = element_text(size = 20))
g

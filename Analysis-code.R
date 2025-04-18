rm(list=ls())
library(lme4)
library(car)
library(plyr)
library(vegan)
library(dplyr)
library(picante)
library(ieggr)
library(ggcor)
library(psych)
library(corrplot)
library(devtools)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(reshape2)
library(tidyverse)
library(scales)
library(GGally)
library(cowplot)
library(betapart)
library(nlme)
library(emmeans)
library(tidyr) 
save.wd <- iwd(choose.dir())

###### Environmental variable analysis-d0 -------------------------------
env <- read.csv("Env.csv",sep=",",header = TRUE,row.names = 1)
env.d0 <- env[(env$treat=="CK"),]

### linear mixed-effects models
output <- c()
for (i in 6:ncol(env.d0))
{
  env.test <- env.d0[,c(1:5,i)]
  
  for (site in unique(env.test$site))
  {
    env.test1 <- env.test[env.test$site != site,]
    env.test1$treat <- 1
    env.test1[8:14,]$treat <- 0
    env.test1 <- na.omit(env.test1)
    
    fm1<-lmer(scale(env.test1[,6])~treat+(1|Date),data=env.test1)
      
      presult<-car::Anova(fm1,type=2)
      coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
      names(coefs)<-paste0(names(coefs),".mean")
      
      SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
      names(SEvalues)<-paste0(names(SEvalues),".se")
      
      tvalues<-coef(summary(fm1))[ , "t value"] ##t values
      names(tvalues)<-paste0(names(tvalues),".t")
      
      chisqP<-c(presult[,1],presult[,3])
      names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
      
      result <- as.data.frame(c(coefs,tvalues,SEvalues,chisqP))
      colnames(result) <- paste(colnames(env.d0)[i],unique(env.test1$site)[1],
                                unique(env.test1$site)[2],"CK",sep = "_")
      result <- as.data.frame(t(result))
      output <- rbind(output,result)
  }
}

### FDR correction 
output1 <- c()
for (i in 1:(nrow(output)/3)){
  data1 <- output[c((3*i-2):(3*i)),]
  data1$treat.P.fdr <- p.adjust(as.numeric(data1$treat.P),method="fdr")
  output1 <- rbind(output1,data1)
}

###### Treatment effects on environmental variables -------------------------------
treatment <- read.csv("Treatment.csv",sep=",",header = TRUE, row.names=1)
treatment <- treatment[3:86,]
Env <- read.csv("Env.csv",sep=",",header = TRUE, row.names=1)
spc <- match.name(rn.list = list(treatment=treatment,Env=Env))
Env <- spc$Env

### linear mixed-effects models
Env.test <- Env
output <- c()
for (site in unique(Env.test$site))
{
  Env.test1 <- Env.test[Env.test$site == site,]
  for (treat in unique(Env.test1[(Env.test1$treat!="CK"),]$treat))
  {
    Env.test2 <- Env.test1[(Env.test1$treat=="CK")|(Env.test1$treat==treat),]
    Env.test2$treatment <- 0
    Env.test2[(Env.test2$treat!="CK"),]$treatment <- 1
    for (i in 6:(ncol(Env.test2)-1))
    {
      fm1<-lmer(scale(Env.test2[,i])~treatment+(1|Date),data=Env.test2)
      
      presult<-car::Anova(fm1,type=2)
      coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
      names(coefs)<-paste0(names(coefs),".mean")
      
      SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
      names(SEvalues)<-paste0(names(SEvalues),".se")
      
      tvalues<-coef(summary(fm1))[ , "t value"] ##t values
      names(tvalues)<-paste0(names(tvalues),".t")
      
      chisqP<-c(presult[,1],presult[,3])
      names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
      
      result <- as.data.frame(c(coefs,tvalues,SEvalues,chisqP))
      colnames(result) <- paste(site,treat,colnames(Env.test2)[i], sep = "_")
      result <- as.data.frame(t(result))
      output <- rbind(output,result)
    }
  }
}

###### Treatment effects on alpha diversity -------------------------------
dat <- read.table("Bacteria-Fungi-zOTUtable-alphadiversity.csv",sep=",",header = TRUE, row.names=1)

### linear mixed-effects models
output <- c()
for (site in unique(dat$site))
{
  dat1 <- dat[dat$site == site,]
  for (treat in unique(dat1[(dat1$treat!="CK"),]$treat))
  {
    dat2 <- dat1[(dat1$treat=="CK")|(dat1$treat==treat),]
    dat2$treatment <- 0
    dat2[(dat2$treat!="CK"),]$treatment <- 1
    for (i in 6:(ncol(dat2)-1))
    {
      fm1<-lmer(scale(dat2[,i])~treatment+(1|Date),data=dat2)
      
      presult<-car::Anova(fm1,type=2)
      coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
      names(coefs)<-paste0(names(coefs),".mean")
      
      SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
      names(SEvalues)<-paste0(names(SEvalues),".se")
      
      tvalues<-coef(summary(fm1))[ , "t value"] ##t values
      names(tvalues)<-paste0(names(tvalues),".t")
      
      chisqP<-c(presult[,1],presult[,3])
      names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
      
      result <- as.data.frame(c(coefs,tvalues,SEvalues,chisqP))
      colnames(result) <- paste(site,treat,colnames(dat2)[i], sep = "_")
      result <- as.data.frame(t(result))
      output <- rbind(output,result)
    }
  }
}

###### Treatment effects on beta diversity - Adonis/NMDS -------------------------------
treatment <- read.csv("Treatment.csv",sep=",",header = TRUE, row.names=1)
treatment <- treatment[3:86,]
### Bacteria
otu <- read.table("Bacteria-zOTUtable.csv",sep=",",header = TRUE, row.names=1)
### Fungi
otu <- read.table("Fungi-zOTUtable.csv",sep=",",header = TRUE, row.names=1)

otu <- as.data.frame(otu[,1:86])
spc <- match.name(rn.list = list(treatment=treatment),cn.list = list(otu=otu))
comm <- spc$otu
comm = as.data.frame(t(comm))
dim(comm)
comm = comm[, colSums(comm) > 0]
dim(comm)

dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

### Adonis - whole community 
dist.bray = vegdist(comm,binary=FALSE)
adonis.bray = adonis2(dist.bray ~ treat+site*day, data = treatment, permutations = 999)

### Adonis - each site
dist.bray.site = vegdist(comm[treatment$site=="P",],binary=FALSE)
adonis.bray = adonis2(dist.bray.site ~ treat+day, 
                      data = treatment[treatment$site=="P",], permutations = 999)

###### Network construction - LSA-------------------------------
### Bacteria
otu <- read.table("Bacteria-zOTUtable.csv",sep=",",header = TRUE, row.names=1)
### Fungi
otu <- read.table("Fungi-zOTUtable.csv",sep=",",header = TRUE, row.names=1)
treatment <- read.csv("Treatment.csv",sep=",",header = TRUE, row.names=1)
treatment <- treatment[3:86,]

unique(treatment$all)

treat <- treatment[treatment$all=="PB",]
spc <- match.name(rn.list = list(treat=treat),cn.list = list(otu=otu))
comm <- spc$otu
dim(comm)
comm = comm[, colSums(comm) > 0]
dim(comm)
comm = comm[rowSums(comm) > 0,]
dim(comm)

counts<-rowSums(comm>0)
### Bacteria
comm<-comm[counts>=4,]
### Fungi
comm<-comm[counts>=3,]
write.table(comm, file = "Bacteria_PB_47.txt", sep = "\t", row.names = TRUE)

##### Input data reshape and basic information
data <- read.csv("Bacteria_PB.csv",sep=",",header = TRUE)
data <- data[c('X', 'Y', 'LS')]
data <- col3.dist(data,to.dist=FALSE)
data <- as.data.frame(data)
data[is.na(data)] <- 0

### The function "net.pro" can be used to calculate the main topological properties of a network
### including: Total nodes,Total links,Average degree (avgK),Average weighted degree (avgKw), Average clustering coefficient (avgCC)
### including: Average path distance (GD), Density (D),Transitivity (Trans),Krackhardt Connectedness (Con),Power-law distribution
net.pro = function (assmc)
{
  require(igraph)
  amg = abs(as.matrix(assmc))
  amg[is.na(amg)] = 0
  diag(amg) = 0
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  nod.name = V(ag)$name
  centr.deg = igraph::centr_degree(ag, mode = "all", loops = FALSE, 
                                   normalized = TRUE)
  nod.deg = centr.deg$res
  net.deg = mean(nod.deg)
  nod.wdeg = igraph::strength(ag, mode = "all", loops = FALSE)
  net.wdeg = mean(nod.wdeg)
  require(sna)
  nod.transit = igraph::transitivity(ag, type = "local")
  net.noden = nrow(assmc)
  net.edgen = sum(abs(amg[upper.tri(amg)]) > 0)
  net.edgen.pos = sum(rowSums(assmc>0))/2
  net.edgen.neg = sum(rowSums(assmc<0))/2
  net.edgen.posp = net.edgen.pos/net.edgen
  net.edgen.negp = net.edgen.neg/net.edgen
  net.meandis = igraph::mean_distance(ag, directed = FALSE)
  net.density = igraph::edge_density(ag, loops = FALSE)
  net.cc = mean(nod.transit, na.rm = TRUE)
  net.transit = igraph::transitivity(ag, type = "global")
  net.connect = sna::connectedness(amg)
  fitpowerlaw = function(graph) {
    d = igraph::degree(graph, mode = "all")
    dd = igraph::degree_distribution(graph, mode = "all", 
                                     loops = FALSE, cumulative = FALSE)
    degree = 1:max(d)
    probability = dd[-1]
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    lmijCS=Anova(reg,type = "II")
    p = lmijCS$`Pr(>F)`[[1]]
    list(alpha = alpha, R.square = R.square,p = p)
  }
  net.plr2 = fitpowerlaw(ag)$R.square
  net.alpha = fitpowerlaw(ag)$alpha
  net.p = fitpowerlaw(ag)$p

  net.att = data.frame(NetworkIndex = c("Total nodes","Total links","Total positive links","Total negative links","Total positive links proportion","Total negative links proportion", 
                                        "Average degree (avgK)", "Average weighted degree (avgKw)", 
                                        "Average clustering coefficient (avgCC)", "Average path distance (GD)", 
                                        "Density (D)","Transitivity (Trans)", "Krackhardt Connectedness (Con)",
                                        "alpha","r","P"), 
                       Value = c(net.noden, net.edgen,net.edgen.pos,net.edgen.neg,net.edgen.posp,
                                 net.edgen.negp,net.deg, net.wdeg,  
                                 net.cc, net.meandis, net.density,net.transit, net.connect,
                                 net.alpha,net.plr2,net.p), 
                       stringsAsFactors = FALSE)
  net.att
}
### Taking the Bacteria_PB as an example
result1 = net.pro(data)

### The function "net.rand" can be used to generate random networks and compare them with empirical networks
net.rand = function (assmc)
{
  require(igraph)
  amg = abs(as.matrix(assmc))
  amg[is.na(amg)] = 0
  diag(amg) = 0
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  neta.obs = net.pro(assmc)
  modul.obs = module(assmc = assmc, absolute = TRUE, methods = c("greedy"))
  mod.obs = modul.obs$module
  mod.node = modul.obs$node
  
  netp.obs = rbind(neta.obs, data.frame(NetworkIndex = c(paste0("Module.number.", mod.obs$Method), paste0("Modularity.", mod.obs$Method)), 
                                        Value = c(mod.obs$Module.number, mod.obs$Modularity), 
                                        stringsAsFactors = FALSE))
  rand.method = c("swap")
  randone <- function(ag, rand.method, swap.iter = 100, 
                      endpoint.p = 0.5) {
    require(igraph)
    if (rand.method == "swap") {
      randm = igraph::rewire(ag, with = igraph::keeping_degseq(loops = FALSE, 
                                                               niter = igraph::vcount(ag) * 10))
      randm = igraph::set_edge_attr(graph = randm, 
                                    name = "weight", value = E(ag)$weight)
    }
    else if (rand.method == "endpoints") {
      randm = igraph::rewire(ag, with = igraph::each_edge(prob = endpoint.p, 
                                                          loops = FALSE, multiple = FALSE))
    }
    else {
      stop("rand.method is not valid.")
    }
    igraph::as_adj(randm, sparse = FALSE, names = TRUE, 
                   attr = "weight")
  }
  netp.rand <- sapply(1:100, function(i) {
    message("Now i=", i, " in ", 100, ". ", date())
    ramgi = randone(ag = ag, rand.method = rand.method, 
                    swap.iter = swap.iter, endpoint.p = endpoint.p)
    netai = net.pro(ramgi)
    moduli = module(assmc = ramgi, absolute = TRUE,methods = c("greedy"))
    modi = moduli$module
    as.numeric(c(netai[,2], modi$Module.number, 
                 modi$Modularity))
  })
  z.test = function(a, mu, two.tail = FALSE) {
    zeta = (mean(a, na.rm = TRUE) - mu)/sd(a, na.rm = TRUE)
    if (two.tail) {
      p = 2 * pnorm(-abs(zeta))
    }
    else {
      p = pnorm(-abs(zeta))
    }
    list(z = zeta, p = p)
  }
  random.mean = rowMeans(netp.rand, na.rm = TRUE)
  random.sd = apply(netp.rand, 1, sd, na.rm = TRUE)
  z.value = (random.mean - as.numeric(netp.obs$Value))/random.sd
  p.ztest = pnorm(-abs(z.value))
  EPS <- sqrt(.Machine$double.eps)
  p.lower = (rowSums(netp.rand >= (matrix(as.numeric(netp.obs$Value), 
                                          nrow = nrow(netp.rand), ncol = ncol(netp.rand)) - EPS), 
                     na.rm = TRUE) + 1)/(1 + rowSums(!is.na(netp.rand)))
  p.higher = (rowSums(netp.rand <= (matrix(as.numeric(netp.obs$Value), 
                                           nrow = nrow(netp.rand), ncol = ncol(netp.rand)) + EPS), 
                      na.rm = TRUE) + 1)/(1 + rowSums(!is.na(netp.rand)))
  p.count = apply(cbind(p.lower, p.higher), 1, min)
  out = data.frame(netp.obs, random.mean = random.mean, random.sd = random.sd, 
                   z.value = z.value, p.ztest = p.ztest, p.count = p.count,stringsAsFactors = FALSE)
  colnames(out) = c("NetworkIndex", "Empirical.Results", "Random.Mean", 
                    "Random.Stdev", "Z.value", "P.Ztest", "P.count")
  out
}
### The function “module” can be used to calculate network modularity-related parameters and identify keystone nodes
module = function (assmc, methods = c("greedy", "walk", "eigen", "betweenness","infomap", "MLopt"), absolute = TRUE, zi.threshold = 2.5, 
                   Pi.threshold = 0.62) 
{
  require(igraph)
  if (absolute) {
    amg = abs(as.matrix(assmc))
  }
  else {
    amg = as.matrix(assmc)
  }
  amg[is.na(amg)] = 0
  diag(amg) = 0
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  mods = list()
  if ("greedy" %in% methods) {
    mods$Greedy.optimization = try(igraph::cluster_fast_greedy(graph = ag, 
                                                               merges = TRUE, modularity = TRUE, membership = TRUE, 
                                                               weights = E(ag)$weight))
  }
  if ("eigen" %in% methods) {
    mods$Leading.eigenvector = try(igraph::cluster_leading_eigen(ag))
  }
  if ("walk" %in% methods) {
    mods$Short.random.walk = try(igraph::cluster_walktrap(ag, 
                                                          weights = E(ag)$weight, steps = 4, merges = TRUE, 
                                                          modularity = TRUE, membership = TRUE))
  }
  if ("betweenness" %in% methods) {
    mods$Edge.betweenness = try(igraph::cluster_edge_betweenness(ag, 
                                                                 weights = E(ag)$weight, directed = FALSE, edge.betweenness = TRUE, 
                                                                 merges = TRUE, bridges = TRUE, modularity = TRUE, 
                                                                 membership = TRUE))
  }
  if ("infomap" %in% methods) {
    mods$infomap = try(igraph::cluster_infomap(ag, nb.trials = 10, 
                                               modularity = TRUE))
  }
  if ("MLopt" %in% methods) {
    mods$Multi.level.optimization = try(igraph::cluster_louvain(ag))
  }
  if (length(mods) == 0) 
    stop("All method names are not valid.")
  module.att = data.frame(Method = names(mods), Module.number = sapply(1:length(mods), 
                                                                       function(i) {
                                                                         if (inherits(mods[[i]], "try-error")) {
                                                                           out = NA
                                                                         }
                                                                         else {
                                                                           out = length(mods[[i]])
                                                                         }
                                                                         out
                                                                       }), Modularity = sapply(1:length(mods), function(i) {
                                                                         if (inherits(mods[[i]], "try-error")) {
                                                                           out = NA
                                                                         }
                                                                         else {
                                                                           out = modularity(mods[[i]])
                                                                         }
                                                                         out
                                                                       }), stringsAsFactors = FALSE)
  rownames(module.att) = c()
  zpvalue <- function(mod, amg, method.name) {
    if (inherits(mod, "try-error")) {
      out = data.frame(Method = method.name, ID = rownames(amg), 
                       ModuleID = NA, zi = NA, Pi = NA, Classification = NA, 
                       stringsAsFactors = FALSE)
    }
    else {
      nod.modid = igraph::membership(mod)
      modlen = length(mod)
      zi <- kis2 <- rep(0, length(nod.modid))
      ki = rowSums(amg > 0)
      for (i in 1:modlen) {
        idi = which(nod.modid == i)
        nodi = names(nod.modid)[idi]
        idii = match(nodi, rownames(amg))
        kib = rowSums(amg[idii, idii, drop = FALSE] > 
                        0)
        zib = (kib - mean(kib))/sd(kib)
        zib[which(is.na(zib))] = 0
        zi[idi] = zib
        kisi = rowSums(amg[, idii, drop = FALSE] > 0)
        kis2 = kis2 + (kisi^2)
      }
      Pi = 1 - (kis2/(ki^2))
      nod.class = rep("Peripheral", length(nod.modid))
      nod.class[which(zi <= zi.threshold & Pi > Pi.threshold)] = "Connectors"
      nod.class[which(zi > zi.threshold & Pi <= Pi.threshold)] = "Module.hubs"
      nod.class[which(zi > zi.threshold & Pi > Pi.threshold)] = "Network.hubs"
      out = data.frame(Method = method.name, ID = names(nod.modid), 
                       ModuleID = as.vector(nod.modid), zi = zi, Pi = Pi, 
                       Classification = nod.class, stringsAsFactors = FALSE)
    }
    rownames(out) = c()
    out
  }
  node.mod.att = Reduce(rbind, lapply(1:length(mods), function(i) {
    zpvalue(mod = mods[[i]], amg = amg, method.name = names(mods)[i])
  }))
  list(module = module.att, node = node.mod.att)
}

### Taking the Bacteria_PB as an example
result2 = net.rand(data) 
modul.obs = module(data, absolute = TRUE, methods = c("greedy"))
mod.node = modul.obs$node

##### Network stability analysis - Random removal
data.raw<-data[colSums(abs(data))>0,colSums(abs(data))>0]
### relative abundance of each species
### Taking the Bacteria_PB as an example
otu <- read.table("Bacteria-zOTUtable.csv",sep=",",header = TRUE, row.names=1)
unique(treatment$all)

treat <- treatment[treatment$all=="PB",]
spc <- match.name(rn.list = list(treat=treat),cn.list = list(otu=otu))
ra <- spc$otu
ra <- ra/colSums(ra)
spc <- match.name(rn.list = list(data.raw=data.raw,ra=ra))
ra <- spc$ra
ra.list <- colMeans(t(ra))
sum(row.names(data.raw)==names(ra.list))

### Function
rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=data.raw, rm.p.list=seq(0.1,1,by=0.1), sp.ra=ra.list, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=data.raw, rm.p.list=seq(0.1,1,by=0.1), sp.ra=ra.list, abundance.weighted=F,nperm=100)

##### Network stability analysis - Natural connectivity
### Index calculation
nc <- function(adj_matrix) {
  adj_matrix <- as.matrix(adj_matrix)
  adj_matrix[abs(adj_matrix) != 0] <- 1
  
  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)
  
  lambda_sum <- 0
  N = length(lambda)
  for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum/N, base = exp(1))
  lambda_average
}

### Simulation
nc.simu <- function(netRaw, rm.p.list,nperm=100){
  t(sapply(rm.p.list,function(x){
    nc1=sapply(1:nperm,function(i){
      id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*x))
      net.Raw=netRaw
      net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;
      net.Raw<-net.Raw[colSums(abs(net.Raw))>0,colSums(abs(net.Raw))>0]
      nc(net.Raw)
    })
    nc1.mean=mean(nc1)
    nc1.sd=sd(nc1)
    nc1.se=sd(nc1)/(nperm^0.5)
    result<-c(nc1.mean,nc1.sd,nc1.se)
    names(result)<-c("Natural.connectivity.mean","Natural.connectivity.sd",
                     "Natural.connectivity.se")
    result
  }))
}

result.nc <- nc.simu(netRaw=data.raw,rm.p.list=seq(0,0.9,by=0.1),nperm=100)

###### Niche width - Generalist -------------------------------
treatment <- read.csv("Treatment.csv",sep=",",header = TRUE, row.names=1)
treatment <- treatment[3:86,]
### Bacteria
otu <- read.table("Bacteria-zOTUtable.csv",sep=",",header = TRUE, row.names=1)
### Fungi
otu <- read.table("Fungi-zOTUtable.csv",sep=",",header = TRUE, row.names=1)

otu <- as.data.frame(otu[,1:86])
spc <- match.name(rn.list = list(treatment=treatment),cn.list = list(otu=otu))
comm <- spc$otu
comm = as.data.frame(t(comm))

### Niche width
comm1 <- comm[treatment$site=="P",]
comm1 = comm1[, colSums(comm1) > 0]
dim(comm1)
library(spaa)
niche_width <- niche.width(comm1, method = 'levins')
min(niche_width)
max(niche_width)

### Generalists or specialists 
library(EcolUtils)
set.seed(1028)
result <- spec.gen(comm1, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))

###### Structural equation model - PLS-PM-------------------------------
library(plspm)
data <- read.csv("Bacteria-Fungi-plspm.csv",sep=",",header = TRUE)
data[,3:ncol(data)] <- scale(data[,3:ncol(data)])
data$SCF <- data$SCF*(-1)
data$Bac.niche <- data$Bac.niche*(-1)
data$Fun.niche <- data$Fun.niche*(-1)
data$Bac.robustness <- data$Bac.robustness*(-1)
data$Fun.negative <- data$Fun.negative*(-1)

colnames(data)
Latent.Variable <- list(
  env = c("SOC","TC","TN","SCF"),
  Bac.diversity = "Bac.rich",
  Fun.diversity = "Fun.rich",
  Bac.function = c("Bac.rrn","Bac.niche"),
  Fun.function = c("Fun.plant","Fun.niche"),
  Bac.network = c("Bac.negative","Bac.robustness"),
  Fun.network = c("Fun.negative","Fun.robustness"),
  env.es = c("SOC.SE","TC.SE","TN.SE","SCF.SE")
)

env <- c(0,0,0,0,0,0,0,0)
Bac.diversity <- c(1,0,0,0,0,0,0,0)
Fun.diversity <- c(1,0,0,0,0,0,0,0)
Bac.function <- c(0,1,0,0,0,0,0,0)
Fun.function <- c(0,0,1,0,0,0,0,0)
Bac.network <- c(1,0,0,1,0,0,0,0)
Fun.network <- c(0,0,0,0,1,0,0,0)
env.es <- c(0,1,0,0,1,1,1,0)

path <- rbind(env,Bac.diversity,Fun.diversity,Bac.function,Fun.function,
              Bac.network,Fun.network,env.es)
colnames(path) <- rownames(path)

modes <- rep('A',8)

pls <- plspm(data,path,Latent.Variable,modes=modes)

#Goodness-of-fit > 0.7
pls$gof
#loading > 0.7
outer_model <- pls$outer_model

inner_model <- pls$inner_model
path_coefficients <- bind_rows(lapply(names(inner_model), function(name) {  data.frame(    Predictor = rownames(inner_model[[name]]),    Response = name,    Estimate = inner_model[[name]][, "Estimate"],    StdError = inner_model[[name]][, "Std. Error"],    tValue = inner_model[[name]][, "t value"],    pValue = inner_model[[name]][, "Pr(>|t|)"]  )})) %>%  filter(Predictor != "Intercept")
r_squared <- as.data.frame(pls$inner_summary) %>%  select(Type, R2) %>%  filter(Type == "Endogenous") %>%  mutate(Response = rownames(.)) %>%  select(Response, R2)
path_coefficients <- path_coefficients %>%  left_join(r_squared, by = "Response") %>%  arrange(Response, Predictor) %>%  mutate(Significance = ifelse(pValue < 0.001, "***",                               ifelse(pValue < 0.01, "**",                                      ifelse(pValue < 0.05, "*", "ns"))))
variable_effect <- as.data.frame(pls$effects)

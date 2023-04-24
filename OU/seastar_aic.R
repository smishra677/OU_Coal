


library(matlib)
library(seastaR)
library(ape)
library(dplyr)
library(gtools)
library(combinat)
library(geiger)
library(stringr)
library(ouch)


.prepare.bm.univariate=function(phy, dat, nodes=NULL, SE=NA, control=list(binary=TRUE, ultrametric=FALSE)){

  ## CONTROL OBJECT
  ct=list(binary=TRUE, ultrametric=FALSE)
  ct[names(control)]=control

  print(1)
  ## MATCHING and major problems
  td=treedata(phy, dat, sort=TRUE, warnings=FALSE)
  phy=reorder(td$phy, "postorder")
  if(ct$binary) if(!is.binary.phylo(phy)) stop("'phy' should be a binary tree")
  if(ct$ultrametric) if(!is.ultrametric(phy)) stop("'phy' should be an ultrametric tree")
  if(is.null(phy$edge.length)) stop("'phy' must include branch lengths in units of time")

  if(ncol(td$data)>1) stop("'dat' should be univariate")
  dat=td$data[,1]

  ## RESOLVE SE
  seTMP=structure(rep(NA, length(dat)), names=names(dat))

  if(is.null(SE)) SE=NA

  if(length(SE)>1){
    if(is.null(names(SE))) stop("'SE' should be a named vector")
    if(!all(names(dat)%in%names(SE))) stop("names in 'SE' must all occur in names of 'dat'")
    seTMP[names(SE[names(dat)])]=SE[names(dat)]
    SE=seTMP
  } else {
    if(is.numeric(SE)){
      seTMP[]=SE
      SE=seTMP
    } else {
      SE=seTMP
    }
  }

  if(!all(is.na(SE) | SE >= 0)) stop("'SE' values should be positive (including 0) or NA")

  ## CACHE tree
  cache=.cache.tree(phy)
  N=cache$n.tip
  n=cache$n.node
  m<-s<-g<-numeric(N+n)

  ## RESOLVE data: given trait values (m and g) and SE (s) for every node (tips and internals)
  g[1:N]=1
  m[]=NA; m[1:N]=dat
  s[1:N]=SE

  ## RESOLVE nodes
  if(!is.null(nodes)){
    nn=(N+1):(N+n)
    vec=.cache.y.nodes(m, s, g, nn, phy, nodes=nodes)
  } else {
    vec=rbind(m=m, s=s)
    attr(vec, "given")=g
    attr(vec, "adjse")=as.numeric(is.na(s))[1:N]
  }

  cache$SE=SE
  cache$dat=dat[match(phy$tip.label, names(dat))]
  cache$phy=phy

  cache$y=vec

  return(cache)
}




sigma2_AIC <- function(sigma2,trait, var_covar) {

  one <-c(1, 1, 1)
  zhat_root <- (t(one)%*%inv(var_covar)%*%one)%*%(t(one)%*%inv(var_covar)%*%trait)
  top <- exp((-1/2)*((t(trait - zhat_root*one))%*%inv(as.vector(sigma2)*var_covar)%*%(trait - zhat_root*one)))
  bottom = sqrt(((2*pi)^3)*det(as.vector(sigma2)*var_covar))
  logL = log(top/bottom)

  return(logL)
}


.ou.vcv<-function(vcv){
  fx=function(alpha){
    vcvDiag <- diag(vcv)
    diagi <- matrix(vcvDiag, nrow = length(vcvDiag), ncol = length(vcvDiag))
    diagj <- matrix(vcvDiag, nrow = length(vcvDiag), ncol = length(vcvDiag), byrow = TRUE)
    Tij = diagi + diagj - (2 * vcv)
    vcvRescaled = (1/(2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcv))
    return(vcvRescaled)
  }
  fx
}


#
ou.lik <- function (phyloe,SE=0,dat)
{
  model="OU"
  #if(is.na(SE)){SE=0}
  cache = .prepare.bm.univariate(phy, dat, SE = SE)
  cache$ordering = attributes(cache$phy)$order
  cache$N = cache$n.tip
  cache$n = cache$n.node
  cache$nn = (cache$root + 1):(cache$N + cache$n)
  cache$intorder = as.integer(cache$order[-length(cache$order)])
  cache$tiporder = as.integer(1:cache$N)
  cache$z = length(cache$len)
  ll.ou.vcv <- function(cache, alpha, sigsq, z0, se){
    N <- cache$N

    #phyvcv <-vcv.phylo(cyp_phy)

    #ouV <- .ou.vcv(C_matrix_seastaR)

    phyvcv <- vcv.phylo(cache$phy)
    ouV <- .ou.vcv(phyloe)


    V <- sigsq*ouV(alpha)
    if(!is.null(se)){
      diag(V) <- diag(V)+se^2
    } else {
      if(any(attr(cache$y, "given")==1)){
        var <- cache$y['s',][attributes(cache$y)$given==1]^2
        var[is.na(var)] <- 0
        diag(V) <- diag(V)+var
      }
    }
    ci <- solve(V)
    invV <- ci
    if(is.null(z0)){
      o <- rep(1, length(cache$dat))
      m1 <- 1/(t(o) %*% ci %*% o)
      m2 <- t(o) %*% ci %*% cache$dat
      mu <- m1 %*% m2
    } else {mu <- z0}
    logd <- determinant(V, logarithm=TRUE)
    logl <- -0.5 * (t(cache$dat-mu) %*% invV %*%
                      (cache$dat-mu)) - 0.5 * logd$modulus -
      0.5 * (N * log(2 * pi))
    attributes(logl) <- NULL
    if (is.na(logl)) logl = -Inf
    attr(logl, "ROOT.MAX") = mu
    class(logl) = c("glnL", class(logl))
    return(logl)
  }
  class(ll.ou.vcv) <- c("bm","dtlik","function")
  fx_exporter = function() {
    attb=c("alpha","sigsq")
    cache$attb <- attb
    if (any(attr(cache$y, "adjse") == 1)) {
      attb = c(attb, "SE")
    }
    lik <- function(pars, root="max", ...){
      attb=c("alpha","sigsq")
      cache$attb <- attb
      if (any(attr(cache$y, "adjse") == 1)) {
        attb = c(attb, "SE")
      }
      if (root == "max"){
        rtmx = TRUE
      } else {
        if (root %in% c("obs", "given")) {
          attb = c(attb,"z0")
        } else {
          stop("unusable 'root' type specified")
        }
      }
      if (missing(pars))
        stop(paste("The following 'pars' are expected:\n\t",
                   paste(attb, collapse = "\n\t", sep = ""), sep = ""))
      pars = .repars(pars, attb)
      names(pars) = attb
      if ("alpha" %in% attb)
        alpha = pars[["alpha"]]
      sigsq = pars[["sigsq"]]
      if ("SE" %in% attb){
        se = pars[["SE"]]
      } else se = NULL
      if ("z0" %in% attb) {
        z0 = pars[["z0"]]
      } else {z0 = NULL}
      ll = ll.ou.vcv(cache = cache, alpha=alpha, sigsq = sigsq,
                     z0=z0, se = se)
      return(ll)
    }
    attr(lik, "argn") = attb
    attr(lik, "cache") <- cache
    class(lik) = c("bm", "function")
    lik
  }
  likfx <- fx_exporter()
  return(likfx)
}
#
#






#write.csv(df, "C:\\Users\\smish\\Documents\\OU\\datasets\\Mat.csv", row.names=TRUE,col.names = FALSE)


cyp_phy <- read.tree("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt")
as.matrix(cyp_phy$edge)
cyp_phy$edge.length
cyp_phy$tip.label

species <- seastaR::parse_input_file("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt", genetrees = FALSE)






datalist<-list("sp1", "sp2","sp3")

C_matrix <- ape::vcv(species) #Standard species tree matrix

C_matrix



C_matrix_seastaR <- seastaR::get_full_matrix(species) #Standard species tree matrix
C_matrix_seastaR


test_trait <- seastaR::simulate_traits(1, C_matrix, 1)



test_trait_seastaR <- seastaR::simulate_traits(1, C_matrix_seastaR, 1)
test_trait_seastaR


phyvcv <-vcv.phylo(cyp_phy)
ouV <- .ou.vcv(C_matrix_seastaR)
ouV(0.3)

#
# df <- data.frame(test_trait_seastaR)
# #df <-t(df)
#
# cyp_dat <- df
#
# cyp_td <- treedata(cyp_phy, cyp_dat)
# states <- cyp_td$data[,"test_trait_seastaR"]
#
# tree <- cyp_td$phy
# states <- states[tree$tip.label]
#
#
#
# ou.lik(C_matrix_seastaR,cyp_dat)
#
#
#
#






# tree1 <- vcv2phylo(C_matrix_seastaR)
# ape::vcv(tree1) #Standard species tree matrix
#
# phyvcv <- vcv.phylo(cyp_phy)
# ouV <- .ou.vcv(phyvcv)
#
# ouV <- .ou.vcv(phyvcv)
#
#
#
#
#
#
# df <- data.frame(test_trait_seastaR)
# #df <-t(df)
#
# cyp_dat <- df
#
#
#
# cyp_td <- treedata(cyp_phy, cyp_dat)
# states <- cyp_td$data[,"test_trait_seastaR"]
#
# tree <- cyp_td$phy
# states <- states[tree$tip.label]
#
#
#
# fit_ou_seastaR <- fitContinuous(tree, states, model="BM")
# fit_ou_seastaR
# fit_ou_seastaR$lik
#
#
# flik_ou_seastaR <- fit_ou_seastaR$lik
#
# siq_bm_seastaR <- fit_ou_seastaR$opt$sigsq
# lik_bm_seastaR1 <- flik_ou_seastaR(c(siq_bm_seastaR ,0), root="given")
# lik_bm_seastaR1
#
#
# sigma2_coal <- seastaR::sigma2_inference(C_matrix_seastaR, test_trait_seastaR)
# sigma2_coal
# lik_coal_coal <- sigma2_AIC(sigma2_coal,test_trait_seastaR,C_matrix_seastaR)
# lik_coal_coal
#
# ee_k_ape<-NULL
# ee_k_ou_ape<-NULL
#
# ee_k_seastaR_bm<-NULL
# ee_k_seastaR_ou<-NULL
#
# ee_k_seastaR<-NULL
#
# ee_k_seastaR_seastaR<-NULL
#
# for (kk in 3:length(datalist)){
#
#   dli<-datalist[1:kk]
#
#   k<-NULL
#   for (i in dli){
#
#     k<-c(k,i)
#
#   }
#
#
#   ee_k_ape<-str_c(c(ee_k_ape,str_c(c(k),collapse = "."),collapse=""))
#   ee_k_ou_ape<-str_c(c(ee_k_ou_ape,str_c(c(k),collapse = "."),collapse=""))
#
#   ee_k_seastaR_bm<-str_c(c(ee_k_seastaR_bm,str_c(c(k),collapse = "."),collapse=""))
#   ee_k_seastaR_ou<-str_c(c(ee_k_seastaR_ou,str_c(c(k),collapse = "."),collapse=""))
#
#   ee_k_seastaR<-str_c(c(ee_k_seastaR,str_c(c(k),collapse = "."),collapse=""))
#   ee_k_seastaR_seastaR<-str_c(c(ee_k_seastaR_seastaR,str_c(c(k),collapse = "."),collapse=""))
#
#   print(ee_k_ape)
#
#   for (i in 1:10)
#   {
#     test_trait <- seastaR::simulate_traits(1, C_matrix, 1)
#     test_trait_seastaR <- seastaR::simulate_traits(1, C_matrix_seastaR, 1)
#
#
#
#     sigma2_ <- seastaR::sigma2_inference(C_matrix, test_trait)
#     sigma2_coal <- seastaR::sigma2_inference(C_matrix_seastaR, test_trait_seastaR)
#
#     lik_coal <- sigma2_AIC(sigma2_,test_trait,C_matrix)
#     lik_coal_coal <- sigma2_AIC(sigma2_coal,test_trait_seastaR,C_matrix_seastaR)
#
#     aic_coal <- (2*(2-lik_coal))
#     aic_coal_coal <- (2*(2-lik_coal_coal))
#
#     df <- data.frame(test_trait)
#     #df <-t(df)
#
#     cyp_dat <- df
#     cyp_td <- treedata(cyp_phy, cyp_dat)
#     states <- cyp_td$data[,"test_trait"]
#
#     tree <- cyp_td$phy
#     states <- states[tree$tip.label]
#
#
#
#     fit_bm <- fitContinuous(tree, states[k], model="BM")
#     fit_ou <- fitContinuous(tree, states[k], model="OU")
#
#     flik_bm <- fit_bm$lik
#     flik_ou <- fit_ou$lik
#
#     aic_bm <- fit_bm$opt$aic
#     lik_bm <- fit_bm$opt$lnL
#     siq_bm <- fit_bm$opt$sigsq
#     lik_bm1<- flik_bm(c(siq_bm ,0), root="given")
#
#
#     aic_ou <- fit_ou$opt$aic
#     lik_ou <- fit_ou$opt$lnL
#     siq_ou <- fit_ou$opt$sigsq
#     lik_ou1<- flik_ou(c(fit_ou$opt$alpha,siq_ou ,0), root="given")
#
#
#
#     df1 <- data.frame(test_trait_seastaR)
#     #df <-t(df)
#
#     cyp_dat1 <- df1
#     cyp_td1 <- treedata(cyp_phy, cyp_dat1)
#     states1 <- cyp_td1$data[,"test_trait_seastaR"]
#
#     tree <- cyp_td$phy
#     states1 <- states1[tree$tip.label]
#
#
#
#
#     fit_bm_seastaR <- fitContinuous(tree, states1[k], model="BM")
#     fit_ou_seastaR <- fitContinuous(tree, states1[k], model="OU")
#
#     flik_bm_seastaR <- fit_bm_seastaR$lik
#     flik_ou_seastaR <- fit_ou_seastaR$lik
#
#
#     aic_bm_seastaR <- fit_bm_seastaR$opt$aic
#     lik_bm_seastaR <- fit_bm_seastaR$opt$lnL
#     siq_bm_seastaR <- fit_bm_seastaR$opt$sigsq
#     lik_bm_seastaR1 <- flik_bm(c(siq_bm_seastaR ,0), root="given")
#
#
#
#     aic_ou_seastaR <- fit_ou_seastaR$opt$aic
#     lik_ou_seastaR <- fit_ou_seastaR$opt$lnL
#     siq_ou_seastaR <- fit_ou_seastaR$opt$sigsq
#     lik_ou_seastaR1<- flik_ou(c(fit_ou_seastaR$opt$alpha,siq_ou_seastaR ,0), root="given")
#
#
#
#
#     tmp1_bm1 <- str_c(c(",(",str_c(aic_bm),";",str_c(siq_bm),";",str_c(lik_bm),";",str_c(lik_bm1),")"),collapse="")
#     ee_k_ape<- str_c(c(ee_k_ape,tmp1_bm1),collapse="")
#
#
#
#     tmp1_bm2 <- str_c(c(",(",str_c(aic_bm_seastaR),";",str_c(siq_bm_seastaR),";",str_c(lik_bm_seastaR),";",str_c(lik_bm_seastaR1),")"),collapse="")
#     ee_k_seastaR_bm<- str_c(c(ee_k_seastaR_bm,tmp1_bm2),collapse="")
#
#
#
#     tmp1_ou3 <- str_c(c(",(",str_c(aic_ou),";",str_c(siq_ou),";",str_c(lik_ou),";",str_c(lik_ou1),")"),collapse="")
#     ee_k_ou_ape<- str_c(c(ee_k_ou_ape,tmp1_ou3),collapse="")
#
#     tmp1_ou4 <- str_c(c(",(",str_c(aic_ou_seastaR),";",str_c(siq_ou_seastaR),";",str_c(lik_ou_seastaR),";",str_c(lik_ou_seastaR1),")"),collapse="")
#     ee_k_seastaR_ou<- str_c(c(ee_k_seastaR_ou,tmp1_ou4),collapse="")
#
#
#
#     tmp1_seastaR5 <- str_c(c(",(",str_c(aic_coal),";",str_c(sigma2_),";",str_c(lik_coal),";",str_c(lik_coal),")"),collapse="")
#     ee_k_seastaR<- str_c(c(ee_k_seastaR,tmp1_seastaR5),collapse="")
#
#     tmp1_seastaR6 <- str_c(c(",(",str_c(aic_coal_coal),";",str_c(sigma2_coal),";",str_c(lik_coal_coal),";",str_c(lik_coal_coal),")"),collapse="")
#     ee_k_seastaR_seastaR<- str_c(c(ee_k_seastaR_seastaR,tmp1_seastaR6),collapse="")
#
#
#
#
#
#   }
#   ee_k_ape<- str_c(c(ee_k_ape,"\n"),collapse="")
#   ee_k_ou_ape<- str_c(c(ee_k_ou_ape,"\n"),collapse="")
#
#   ee_k_seastaR_bm<- str_c(c(ee_k_seastaR_bm,"\n"),collapse="")
#   ee_k_seastaR_ou<- str_c(c(ee_k_seastaR_ou,"\n"),collapse="")
#
#   ee_k_seastaR<- str_c(c(ee_k_seastaR,"\n"),collapse="")
#   ee_k_seastaR_seastaR<- str_c(c(ee_k_seastaR_seastaR,"\n"),collapse="")
#
#
# }
#
#
#
# write.csv(substring(ee_k_ape,1,nchar(ee_k_ape)), "C:\\Users\\smish\\Documents\\OU\\databmbm.csv", row.names = FALSE)
# write.csv(substring(ee_k_ou_ape,1,nchar(ee_k_ou_ape)), "C:\\Users\\smish\\Documents\\OU\\dataoubm.csv", row.names = FALSE)
# write.csv(substring(ee_k_seastaR_bm,1,nchar(ee_k_seastaR_bm)), "C:\\Users\\smish\\Documents\\OU\\databmcoalbm.csv", row.names = FALSE)
#
# write.csv(substring(ee_k_seastaR_ou,1,nchar(ee_k_seastaR_ou)), "C:\\Users\\smish\\Documents\\OU\\databmcoalou.csv", row.names = FALSE)
# write.csv(substring(ee_k_seastaR,1,nchar(ee_k_seastaR)), "C:\\Users\\smish\\Documents\\OU\\databmbmcoal.csv", row.names = FALSE)
# write.csv(substring(ee_k_seastaR_seastaR,1,nchar(ee_k_seastaR_seastaR)), "C:\\Users\\smish\\Documents\\OU\\dataseabmcoalbmcoal.csv", row.names = FALSE)
#
#


library(matlib)
library(seastaR)
library(ape)
library(dplyr)
library(gtools)
library(combinat)
library(geiger)
library(stringr)
library(ouch)
library(subplex)
library(mvtnorm)
library(tidyverse)
library(plotly)
library(MASS)


sigma2_likelihood1 <- function(sigma2, trait, var_covar) {

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
    if(alpha==0)
    {
      vcvRescaled= vcv
    }
    else{
    vcvRescaled = (1/(2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcv))
    }

    return(vcvRescaled)
  }
  fx
}


.ou.vcv_mean<-function(vcv){
  fx=function(alpha,z0,mu){
    vcvDiag <- diag(vcv)
    vcvRescaled = z0* exp(-alpha * vcvDiag) + mu*(1 - exp(-2 * -alpha * vcvDiag))


    return(vcvRescaled)
  }
  fx
}






sigma2_likelihood <- function(param,trait,matr) {
  z0 <- param[1];alpha <-param[2];
  if(alpha!=0)
  {
  ouV <- .ou.vcv(matr)
  mu <- sum(trait)/300;
  var_covar <-ouV(alpha)
  ouM <- .ou.vcv_mean(var_covar)

  mean<-ouM(alpha,z0,mu)

  }
  else
  {
    var_covar<-matr
    mu <- sum(trait)/30000;
    ouM <- .ou.vcv_mean(var_covar)

    mean<-ouM(alpha,z0,mu)
  }
  sigsum<-0
  sigsum1<-0
  vale<-0
  for(i in 1:10000)
  {









  if ( alpha>0 && alpha<15) {
    # one <-c(1, 1, 1)
    # zhat_root <- (t(one)%*%inv(var_covar)%*%one)%*%(t(one)%*%inv(var_covar)%*%trait)
    # print(c(rep(zhat_root,3)))
    # print(det(as.vector(sigma2)*var_covar))
    # top <- exp((-1/2)*((t(trait - zhat_root*one))%*%inv(as.vector(sigma2)*var_covar)%*%(trait - zhat_root*one)))
    # bottom = sqrt(((2*pi)^3)*det(as.vector(sigma2)*var_covar))
    # dmvnorm(trait,rep(zhat_root,3),det(as.vector(sigma2)*var_covar), log=TRUE)
    sigma2<-seastaR::sigma2_inference(var_covar,  trait[i,])
    print(sigma2)
    true_sigma <- seastaR::sigma2_inference(matr,  trait[i,])
    print(true_sigma)
    sigsum1 <-sigsum1+(true_sigma/sigma2)
    if (sigma2<=100*true_sigma || sigma2<=(1/100)*true_sigma)
    {
      sigsum <-sigma2+sigsum
    vale<- vale+dmvnorm(trait[i,],mean,(as.vector(sigma2)*var_covar), log=TRUE)
    }
    else
    {
      10000
    }
}
  else {10000 }
  }
  print('###################')
  print(sigsum/(29999))
  print(sigsum1/29999)
  print('@@@@@@@@@@@@@@@@@@@')
  vale/10000


}


  ou <- function(param,data){
    alpha <- param[1]; sigsq <- param[2];
    ouV <- .ou.vcv(C_matrix_seastaR)
    V <- ouV(alpha)

    if ( alpha>=0 & sigsq>0 ) {
      #dmvnorm(data,c(1, 1, 1),sigsq*V, log=FALSE)
      sigma2_likelihood(sigsq,data,V)
    }
    else { 10000 }
  }





  packageVersion("geiger")
  cyp_phy <- read.tree("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt")




  species <- seastaR::parse_input_file("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt", genetrees = FALSE)


  datalist<-list("sp1", "sp2","sp3")
  C_matrix <- seastaR::get_full_matrix(species) #Standard species tree matrix
  C_matrix1<-C_matrix

  #N<-10000
  alpha<-10
  n<- length(datalist)




  for (i in 1:length(C_matrix))
  {
    if(i%%(n+1)!=1)

    {

      C_matrix[i]=exp(-2*alpha*(C_matrix[(floor(i/(n+1))*(n+1))+1]-C_matrix[i]))*(1-exp(-2*alpha*(C_matrix[i])))/(2*alpha)
    }


  }



  for (i in 1:length(C_matrix))
  {
    if(i%%(n+1)==1)
    {
      C_matrix[i]= (1-exp(-2*alpha*C_matrix[i]))/(2*alpha)
    }
  }



  C_matrix
  C_matrix1


  test_trait1 <- seastaR::simulate_traits(10000, C_matrix, 10)
  test_trait1


  test_trait <- seastaR::simulate_traits(10000, C_matrix1, 10)
  test_trait
#
# test_trait3 <- seastaR::simulate_traits(100, C_matrix, 5)
# test_trait3


d1 <- data.frame(test_trait1)
names(d1)
d1$flag <- rep(1,10000)
d2 <- data.frame(test_trait)
names(d2)
d2$flag <- rep(0,10000)
df3 <- rbind(d1,d2)
df3

p1<- plot_ly(df3, x = ~ sp1, y = ~ sp2, z = ~ sp3,split= ~ flag, label= ~flag,
              marker = list(showscale = TRUE)) %>%

  add_markers()


p1
  # p1<- plot_ly(d2, x = ~ sp1, y = ~ sp2, z = ~ sp3,
  #              marker = list(color='red', showscale = TRUE)) %>%
  #   add_markers()
  #
  # p1
  #
  #

  C_matrix




  optim(c(0,0),sigma2_likelihood,gr=NULL,test_trait,C_matrix,
        method = c("BFGS"),
        control = list(maxit = 1000000,abstol=1e-7), hessian = TRUE)
  test_trait

  C_matrix1



  fit_OU_seastaR <- fitContinuous(tree, states, model="OU")
  print(fit_OU_seastaR)



  # flik_fit_OU_1 <- fit_OU_seastaR$lik
  # lik_ou_seastaR2 <- flik_fit_OU_1(c( 0.07230993, 1.00198148,0), root="given")
  # lik_ou_seastaR2
  #
  # sigma2_likelihood1()
  #
  #   df <- data.frame(test_trait)
  #   cyp_dat <- df
  #
  #   cyp_td <- treedata(cyp_phy, cyp_dat)
  #   states <- cyp_td$data[,"test_trait"]
  #
  #
  #
  #   tree <- cyp_td$phy
  #   states <- states[tree$tip.label]
  #
  #   fit_OU_seastaR <- fitContinuous(tree, states, model="OU")
  #   print(fit_OU_seastaR)
  #
  #
  #   warnings()
  #


library(matlib)
library(seastaR)
library(ape)
library(dplyr)
library(gtools)
library(combinat)
library(geiger)
library(stringr)
library(ouch)



sigma2_AIC <- function(sigma2,trait, var_covar) {

  one <-c(1, 1, 1)
  zhat_root <- (t(one)%*%inv(var_covar)%*%one)%*%(t(one)%*%inv(var_covar)%*%trait)
  top <- exp((-1/2)*((t(trait - zhat_root*one))%*%inv(as.vector(sigma2)*var_covar)%*%(trait - zhat_root*one)))
  bottom = sqrt(((2*pi)^3)*det(as.vector(sigma2)*var_covar))
  logL = log(top/bottom)

  return(2*(2-logL))
}





#write.csv(df, "C:\\Users\\smish\\Documents\\OU\\datasets\\Mat.csv", row.names=TRUE,col.names = FALSE)


cyp_phy <- read.tree("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt")


species <- seastaR::parse_input_file("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt", genetrees = FALSE)





C_matrix_seastaR <- seastaR::get_full_matrix(species)
C_matrix_seastaR


test_trait_seastaR <- seastaR::simulate_traits(1, C_matrix_seastaR, 1)
test_trait_seastaR

sigma2_coal <- seastaR::sigma2_inference(C_matrix_seastaR, test_trait_seastaR)
print(sigma2_coal)

aic_coal_coal <- sigma2_AIC(sigma2_coal,test_trait_seastaR,C_matrix_seastaR)
print(aic_coal_coal)




df1 <- data.frame(test_trait_seastaR)


cyp_dat1 <- df1
cyp_td1 <- treedata(cyp_phy, cyp_dat1)
states1 <- cyp_td1$data[,"test_trait_seastaR"]

tree <- cyp_td$phy
states1 <- states1[tree$tip.label]



fit_bm_seastaR <- fitContinuous(tree, states1, model="BM")
print(fit_bm_seastaR)

fit_ou_seastaR <- fitContinuous(tree, states1, model="OU")


aic_bm_seastaR <- fit_bm_seastaR$opt$aic
siq_bm_seastaR <- fit_bm_seastaR$opt$sigsq

aic_ou_seastaR <- fit_ou_seastaR$opt$aic
siq_ou_seastaR <- fit_ou_seastaR$opt$sigsq








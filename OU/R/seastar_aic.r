
library(matlib)
library(seastaR)
library(ape)
library(dplyr)
library(gtools)
library(combinat)
library(geiger)
library(stringr)
library(ouch)



sigma2_likelihood <- function(sigma2, trait, var_covar) {

  one <-c(1, 1, 1)
  zhat_root <- (t(one)%*%inv(var_covar)%*%one)%*%(t(one)%*%inv(var_covar)%*%trait)
  top <- exp((-1/2)*((t(trait - zhat_root*one))%*%inv(as.vector(sigma2)*var_covar)%*%(trait - zhat_root*one)))
  bottom = sqrt(((2*pi)^3)*det(as.vector(sigma2)*var_covar))
  logL = log(top/bottom)

  return(logL)
}


#write.csv(df, "C:\\Users\\smish\\Documents\\OU\\datasets\\Mat.csv", row.names=TRUE,col.names = FALSE)


cyp_phy <- read.tree("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt")


species <- seastaR::parse_input_file("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt", genetrees = FALSE)



datalist<-list("sp1", "sp2","sp3")
C_matrix <- seastaR::get_full_matrix(species) #Standard species tree matrix



ee_k<-NULL
ee_k_ou<-NULL

for (kk in 3:length(datalist)){

  dli<-datalist[1:kk]

  k<-NULL
  for (i in dli){

    k<-c(k,i)

  }

  ee_k<-str_c(c(ee_k,str_c(c(k),collapse = "."),collapse=""))
  ee_k_ou<-str_c(c(ee_k,str_c(c(k),collapse = "."),collapse=""))

  for (i in 1:10)
  {
    test_trait <- seastaR::simulate_traits(1, C_matrix, 1)


    sptree_rate <- seastaR::sigma2_inference(C_matrix, test_trait)
    lik_star <- sigma2_likelihood(sptree_rate,test_trait,C_matrix)


    df <- data.frame(test_trait)
    #df <-t(df)

    cyp_dat <- df

    #hyphy
    cyp_td <- treedata(cyp_phy, cyp_dat)
    states <- cyp_td$data[,"test_trait"]

    tree <- cyp_td$phy
    states <- states[tree$tip.label]




    fit_bm <- fitContinuous(tree, states[k], model="BM")
    fit_ou <- fitContinuous(tree, states[k], model="OU")


    aic_bm <- fit_bm$opt$aic
    siq_bm <- fit_bm$opt$sigsq
    tmp1_bm <- str_c(c(",(",str_c(aic_bm),"-",str_c(siq_bm),")"),collapse="")
    ee_k<- str_c(c(ee_k,tmp1_bm),collapse="")


    aic_ou <- fit_ou$opt$aic
    siq_ou <- fit_ou$opt$sigsq
    tmp1_ou <- str_c(c(",(",str_c(aic_ou),"-",str_c(siq_ou),")"),collapse="")
    ee_k_ou<- str_c(c(ee_k_ou,tmp1_ou),collapse="")
  }
  ee_k<- str_c(c(ee_k,"\n"),collapse="")
  ee_k_ou<- str_c(c(ee_k_ou,"\n"),collapse="")
}
write.csv(ee_k_ou, "C:\\Users\\smish\\Documents\\OU\\dataou.csv", row.names = FALSE)
write.csv(ee_k, "C:\\Users\\smish\\Documents\\OU\\databm.csv", row.names = FALSE)





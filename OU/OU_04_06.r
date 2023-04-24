
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




#write.csv(df, "C:\\Users\\smish\\Documents\\OU\\datasets\\Mat.csv", row.names=TRUE,col.names = FALSE)


cyp_phy <- read.tree("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt")


species <- seastaR::parse_input_file("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt", genetrees = FALSE)



datalist<-list("sp1", "sp2","sp3")

#C_matrix_seastaR <- seastaR::get_full_matrix(species) #Standard species tree matrix
# C_matrix_seastaR
C_matrix_seastaR <- ape::vcv(species)

pdf <- function(param,data) {
  drift <- param[1];
  if ( sel>=0 & drift>0 ) {
    dmvnorm(data,0,drift*C_matrix_seastaR, log=FALSE)
  }
  else { 10000 }
}


pdf1 <- function(param,data) {
  sel <- param[1]; drift <- param[2]; len <- param[3];
  print((drift/(2*sel))*(1-exp(-2*len*sel)))
  print(C_matrix_seastaR[1])
  print('#########')
  print(sel)
  print(drift)
  print(len)
  print('@@@@@@@@@')
  if (sel<1 & sel>1 & drift<20 & drift>0 & len>0){
    print(1)
  min((drift/(2*sel))*(1-exp(-2*len*sel))-C_matrix_seastaR[1])
}
}


#*(1-exp(-meanslist*sel))

meanslist <- diag(C_matrix_seastaR);

test_trait <- seastaR::simulate_traits(1, C_matrix_seastaR, 2)

subplex(c(1,0,0),pdf1)


optim(c(0),pdf,gr=NULL,test_trait,
      method = c("Brent"),
      lower = -Inf, upper = Inf,
      control = list(maxit = 10000,REPORT=109), hessian = FALSE)
optimHess(c(1,1), pdf, gr = NULL,test_trait)



test_trait <- seastaR::simulate_traits(1, C_matrix_seastaR, 2)

optim(c(1,1,1),pdf1,gr=NULL,
      method = c("Nelder-Mead"),
      lower = -Inf, upper = Inf,
      control = list(maxit = 10000,REPORT=109), hessian = FALSE)



ee_k_ape<-NULL
ee_k_ou_ape<-NULL
ee_k_ou_ape1<-NULL
for (kk in 3:length(datalist)){

  dli<-datalist[1:kk]

  k<-NULL
  for (i in dli){

    k<-c(k,i)

  }


  ee_k_ape<-str_c(c(ee_k_ape,str_c(c(k),collapse = "."),collapse=""))
  ee_k_ou_ape<-str_c(c(ee_k_ou_ape,str_c(c(k),collapse = "."),collapse=""))
  ee_k_ou_ape1<-str_c(c(ee_k_ou_ape1,str_c(c(k),collapse = "."),collapse=""))



  for (i in 1:1)
  {




    print(C_matrix_seastaR)

    test_trait <- seastaR::simulate_traits(1, C_matrix_seastaR, 1)
    print(test_trait)
     sptree_rate <- seastaR::sigma2_inference(C_matrix_seastaR, test_trait)



    like1 <-sigma2_likelihood(sptree_rate,test_trait,C_matrix_seastaR)

    df <- data.frame(test_trait)
    cyp_dat <- df

    cyp_td <- treedata(cyp_phy, cyp_dat)
    states <- cyp_td$data[,"test_trait"]



    tree <- cyp_td$phy
    states <- states[tree$tip.label]


    fit_ou <- fitContinuous1(C_matrix_seastaR,tree, states, model="OU")

    flik_fit_ou <- fit_ou$lik


    lik_ou_seastaR1 <- flik_fit_ou(c(0,fit_ou$opt$sigsq,0), root="given")
    lik_ou_seastaR1



    fit_OU_seastaR <- fitContinuous(tree, states, model="OU")
    print(fit_OU_seastaR)

    fit_BM_seastaR <- fitContinuous(tree, states, model="BM")
    print(fit_BM_seastaR)



    flik_fit_OU_1 <- fit_OU_seastaR$lik
    lik_ou_seastaR2 <- flik_fit_OU_1(c(0,fit_OU_seastaR$opt$sigsq,0 ), root="given")
    lik_ou_seastaR2




    tmp1_bm1 <- str_c(c(",(",str_c(sptree_rate),";",str_c(like1),")"),collapse="")
    ee_k_ape<- str_c(c(ee_k_ape,tmp1_bm1),collapse="")

    tmp1_bm2 <- str_c(c(",(",str_c(fit_ou$opt$sigsq),";",str_c(lik_ou_seastaR1),")"),collapse="")
    ee_k_ou_ape<- str_c(c(ee_k_ou_ape,tmp1_bm2),collapse="")

    tmp1_bm3 <- str_c(c(",(",str_c(fit_OU_seastaR$opt$sigsq),";",str_c(lik_ou_seastaR2),")"),collapse="")
    ee_k_ou_ape1<- str_c(c(ee_k_ou_ape1,tmp1_bm3),collapse="")





  }
  ee_k_ape<- str_c(c(ee_k_ape,"\n"),collapse="")
  ee_k_ou_ape<- str_c(c(ee_k_ou_ape,"\n"),collapse="")
  ee_k_ou_ape1<-str_c(c(ee_k_ou_ape1,"\n"),collapse="")
}



# write.csv(substring(ee_k_ape,1,nchar(ee_k_ape)), "C:\\Users\\smish\\Documents\\OU\\databmbm.csv", row.names = FALSE)
# write.csv(substring(ee_k_ou_ape,1,nchar(ee_k_ou_ape)), "C:\\Users\\smish\\Documents\\OU\\databmou.csv", row.names = FALSE)
# write.csv(substring(ee_k_ou_ape1,1,nchar(ee_k_ou_ape1)), "C:\\Users\\smish\\Documents\\OU\\databmou1.csv", row.names = FALSE)
#
#
#
#
#






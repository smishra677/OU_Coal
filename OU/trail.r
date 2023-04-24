

library(seastaR)
library(ape)
library(dplyr)
library(gtools)
library(combinat)
library(geiger)
library(stringr)

#write.csv(df, "C:\\Users\\smish\\Documents\\OU\\datasets\\Mat.csv", row.names=TRUE,col.names = FALSE)


packageVersion("geiger")
cyp_phy <- read.tree("C:\\Users\\smish\\Documents\\OU\\datasets\\cyprinodon.tre")

print(cyp_phy$tip.label)


genetree_example <- seastaR::parse_input_file("C:\\Users\\smish\\Documents\\OU\\tests\\/seastaR_sptree_example.txt", genetrees = FALSE)



datalist<-list("sp1", "sp2","sp3","sp4", "sp5","sp6","sp7", "sp8","sp9","sp10", "sp11","sp12","sp13", "sp14","sp15","sp16", "sp17","sp18", "sp19","sp20")
C_matrix <- ape::vcv(sptree_example) #Standard species tree matrix 
Cstar_matrix <- seastaR::get_full_matrix(sptree_example) #Gene tree matrix 




dli<-combn(datalist[1:8], 3)



Matrix_x <- matrix(unlist(dli), ncol = 3, byrow = FALSE)


k<-c(Matrix_x[1,1],Matrix_x[1,2])



for (i in 1:length(Matrix_x[,1])){
  k<-NULL
  for (j in 1:length(Matrix_x[1,]))
  {
    k<-c(k,Matrix_x[i,j])
    
  }
  
  cstar_matrix<-Cstar_matrix[k,k]
  test_trait <- seastaR::simulate_traits(1, cstar_matrix, 1)
  df <- data.frame(test_trait)
  df <-t(df)
  
  cyp_dat <- df
  
  
  
  cyp_td <- treedata(cyp_phy[4][k], cyp_dat)
  print(cyp_td)
  states  <- cyp_td$data
  
  tree <- cyp_td$phy
  #print(tree)
  
  
  states <- states[tree$tip.label]
  #print(states)
  
  fit_bm <- fitContinuous(tree, states, model="BM")
  aic_bm <- fit_bm$opt$aic
  lik <- fit_bm$opt$sigsq 
  
  
  
  
}





df <- data.frame(test_trait)
df <-t(df)
datalist<-list("sp1", "sp2","sp3","sp4", "sp5","sp6","sp7", "sp8","sp9","sp10", "sp11","sp12","sp13", "sp14","sp15","sp16", "sp17","sp18", "sp19","sp20")
expand.grid(0:1, 0:1, 0:1)
rownames(df) <- c()unlist(l), nrow=length(l), byrow=TRUE
print(df)

write.csv(df, "C:\\Users\\smish\\Documents\\OU\\datasets\\Mat.csv", row.names=FALSE)



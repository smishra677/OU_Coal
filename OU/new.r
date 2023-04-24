


library(matlib)
library(seastaR)
library(ape)
library(dplyr)
library(gtools)
library(combinat)
library(geiger)
library(stringr)
library(ouch)



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

typeof(C_matrix_seastaR)

df1 <- data.frame(C_matrix_seastaR)

df1[nrow(df1) + 1,] = c(0,0,0)

df1[length(df1) + 1] = c(0,0 ,0,df1[1,1])
df1
colnames(df1)[colnames(df1) == 'V4'] <- '4'

typeof(C_matrix_seastaR)
typeof(data.matrix(df1))

data1 <- data.matrix(df1)
data1

tree1 <- vcv2phylo(data1)



 
cyp_td <- treedata(tree1, cyp_dat)
states <- cyp_td$data[,"test_trait_seastaR"]
 
tree <- cyp_td$phy
states <- states[tree$tip.label]
 
 
 
fit_ou_seastaR <- fitContinuous(tree, states, model="BM")
fit_ou_seastaR






# 
# phyvcv <-vcv.phylo(cyp_phy)
# ouV <- .ou.vcv(C_matrix_seastaR)
# ouV(0.3)

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

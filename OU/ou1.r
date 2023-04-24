

library(geiger)
library(stringr)

packageVersion("geiger")
cyp_phy <- read.tree("C:\\Users\\smish\\Documents\\OU\\datasets\\cyprinodon.tre")


cyp_dat <- read.csv("C:\\Users\\smish\\Documents\\OU\\datasets\\Mat.csv", row.names=1)

head(cyp_dat)


cyp_td <- treedata(cyp_phy, cyp_dat)
cyp_td

ee<-""

for (x in 1:20) {
  tmp <- cbind("V", toString(x))
  k<-str_c(tmp, collapse = "")

  states  <- cyp_td$data[,k]
  tree <- cyp_td$phy
  tree
  
  states <- states[tree$tip.label]
  
  states
  
  fit_bm <- fitContinuous(tree, states, model="BM")
  aic_bm <- fit_bm$opt$aic
  lik <- fit_bm$opt$sigsq 
  temp2 <- cbind(toString(lik),"\n")
  tmp2 <-str_c(temp2,collapse = "")
  
  
  tmp1 <- cbind(toString(aic_bm),tmp2)
  kk<-str_c(tmp1, collapse = ",")
  print(kk)
  ee <- str_c(ee,kk)
  print(ee)


}

print(ee)

write.csv(ee, "C:\\Users\\smish\\Documents\\OU\\data.csv", row.names = TRUE)

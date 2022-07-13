
#----------------------------------------------------------------------# 
# SIMULAION (small n, large p)
#----------------------------------------------------------------------# 
rm(list = ls())
library(dbplyr)
library(dplyr)
library(tidyverse)
library(WGCNA)
library(ggplot2)
library(igraph)
library(MASS)
library(lqmm)
library(mvtnorm)
library(igraph)


Danaher_pos_def <- function(m, cc = 3.4, dd = 0.95){
  AA <- m*cc + t(m*cc)
  AA[AA>0] <- 1
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA <- AA/( as.matrix(rep(1,nrow(m))) ) %*% t(1.4*rowSums(abs(AA)))
  AA <- (AA+t(AA))/2
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA
}

pos_def_func = function(X){
  e_val = eigen(X)$values
  e_vec = eigen(X)$vectors
  neg_e = e_val[e_val < 0]
  if (length(neg_e) > 0) {e_val = e_val + max(abs(neg_e)) + 0.1 }
  else e_val = e_val
  Xnew = e_vec %*% diag(e_val) %*% t(e_vec)#solve(e_vec)
  return(Xnew)
}



pve = function(X, eigen_vec){
  num1 = X %*% eigen_vec
  deno1 = sum(X^2)
  res = sum(num1^2)/deno1
  return(res)
}




set.seed(4591)  #(4591) #(3041) for prob = 0.75

#Erdos-Renyi model
pc_id = 1



prob = 0.1#close to graph density
p =  70# nodes
n = 50 # sample size
n_iter = 100 # iterations
g1 = sample_gnp(p, prob, directed = FALSE, loops = FALSE)
#g1 = make_star(p, mode = "undirected")
graph.density(g1)
E(g1)$weight = runif(length(E(g1)), 0.1, 0.8)
#plot(g1)
#hist(degree(g1), col="lightblue", xlab="Degree", ylab="Frequency", main="")


# adjacency matrix
A = as.matrix(get.adjacency(g1, attr = "weight"))
A = as.matrix(A) + diag(p)

B_true = Danaher_pos_def(A)

if(is.positive.definite(B_true))
{C_true = solve(B_true)

}else {B_true = pos_def_func(A) #Danaher_pos_def(A)
B_true = round(cov2cor(B_true), 4) ###
C_true = solve(B_true)}



Ctest = eigen(C_true)
Ctest_value = Ctest$values
pct = .7
Ctest_value[1] = (pct/(1-pct))*sum(Ctest_value[-1])
Cnew = Ctest$vectors %*% diag(Ctest_value) %*% t(Ctest$vectors)


graph_true = graph_from_adjacency_matrix(B_true, mode = "undirected", weighted = TRUE, diag = FALSE) #B_true


D0 = rowSums(as.matrix(A))
D1 = rowSums(as.matrix(B_true))
L0 = diag(D1) - B_true
#L1 = scale(L0, center = T, scale = T)
TOM_A = TOMsimilarity(B_true, verbose = FALSE)
L2 <- B_true %>% 
  graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE,diag = FALSE) %>% 
  graph.laplacian(normalized = F) %>% as.matrix

alpha = graph.density(graph_true)

#"truth"
#Multivariate normal distribution
mu1 = rep(0,p)

sd_e1 = c(6,3,1)
n_noise = length(sd_e1)

temp0_cor = temp1_cor = temp4_cor = temp0_var = temp1_var = temp4_var = matrix(0, n_iter, n_noise)



for (i in 1:n_iter){
  
  for (e in 1:n_noise){
    
    #print(e)
    sd_e = sd_e1[e]
    
    X0 = mvrnorm(n, mu1, Cnew) 
    #X0 = scale(X0, center = T, scale = T)
    X1 =  X0 + matrix(rnorm(n*p, 0, sd_e), n, p) #observed
    
    X0 = scale(X0, center = T, scale = T)
    X1 = scale(X1, center = T, scale = T)
    
    b1 = c(rnorm(1), D1) 
    e1 = rnorm(n, 0, sd_e) #1) 
    Y0 =   b1[1] + X0 %*% as.matrix(b1[-1]) + e1
    Y0 = scale(Y0, center = T)#, scale = T)
    
    temp0 = prcomp(X0, center = T, scale. = T)
    temp0_var[i,e] = pve(X0, temp0$rotation[,pc_id])
    temp0_cor[i,e] = cor(X0 %*% temp0$rotation[,pc_id], Y0)
    
    
    temp1 = prcomp(X1, center = T, scale. = T)
    temp1_var[i,e] = temp1$sdev[pc_id]^2/sum(temp1$sdev^2) 
    temp1_cor[i,e] = cor(temp1$x[,pc_id], Y0) 
    
    
    temp4  = (1-alpha)*(X1 %*% L2) + alpha*(X1 %*% TOM_A)
    temp4 = prcomp(temp4, center = T, scale. = T)
    temp4_var[i,e] = temp4$sdev[pc_id]^2/sum(temp4$sdev^2)  
    temp4_cor[i,e] = cor(temp4$x[,pc_id], Y0) 
    
    
  }  
}  

colMeans(temp0_cor)
colMeans(temp1_cor)
colMeans(temp4_cor)




#----------------------------------------------------------------------------#
# mean corr plot
# COMBINED PLOT 
#----------------------------------------------------------------------------#
# mean ratio to truth


df = data.frame(Mean_corr = c(mean_all), 
                sd_corr = c(var_all),
                noise = factor(c(rep("high", 12), rep("med", 12), rep("low", 12)), levels = c("high", "med", "low")),
                size = c(rep(c(n,n_all),(n_version - 1))),
                method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
ggplot(df, aes(x = size, y = Mean_corr, ymin = Mean_corr - sd_corr, ymax = Mean_corr + sd_corr, colour = method)) + ylab("Ratio") + ylim(0.5,1.01) +
  geom_point()+ geom_line(aes(linetype = noise)) + scale_linetype_manual(values= c("high" = 1, "med" = 5, "low" = 3)) + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(axis.text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# pc1 explained

df = data.frame(PC1_exp = c(pc_all), 
                sd_corr = c(pc_sd),
                noise = factor(c(rep("high", 12), rep("med", 12), rep("low", 12)), levels = c("high", "med", "low")),
                size = c(rep(c(n,n_all),(n_version - 1))),
                method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))

ggplot(df, aes(x = size, y = PC1_exp, ymin = PC1_exp - sd_corr, ymax = PC1_exp + sd_corr, colour = method)) + ylab("Ratio")+  ylim(0.1,1.1) +
  geom_point()+ geom_line(aes(linetype = noise)) + scale_linetype_manual(values= c("high" = 1, "med" = 5, "low" = 3))  + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(axis.text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 




#----------------------------------------------------------------------# 
# REAL NETWORKS
#----------------------------------------------------------------------# 

rm(list = ls())
library(dbplyr)
library(dplyr)
library(tidyverse)
library(WGCNA)
library(ggplot2)
library(igraph)
library(MASS)
library(lqmm)
library(pls)

pve = function(X, eigen_vec){
  num1 = X %*% eigen_vec
  deno1 = sum(X^2)
  res = sum(num1^2)/deno1
  return(res)
}


#----------------------------------------------------------------------# 
# Adjusted
#----------------------------------------------------------------------# 
# FEV1
setwd("~/Documents/Multiomic project/Subnetwork summarization/MastejEmily/Final FEV1 Adjusted Network Results")
# Emphysema
#setwd("~/Documents/Multiomic project/Subnetwork summarization/MastejEmily/Final pctEmph Adjusted Network Results")

load('final_workspace.RData')
load("SmCCNetWeights.RData")


#----------------------------------------------------------------------# 
# Unadjusted
#----------------------------------------------------------------------# 

#setwd("~/Documents/Multiomic project/Subnetwork summarization/MastejEmily_Unadjusted Net")
#load("trimmed_FEV1_workspace.RData")

#setwd("~/Documents/Multiomic project/Subnetwork summarization/MastejEmily_Unadjusted Net")
#load("trimmed_pctEmph_workspace.RData")

pc_id = 1
idx <- as.logical(moduleColors_trimmed)
net <- as.matrix(X[,idx])
n_trimmed = dim(net)[2]
# this already has both metabolite and protein names
node_name = AbarLabel[idx]


# adjusted module 1 trimmed

CDir = "~/Documents/Multiomic project/Subnetwork summarization/MastejEmily/"
pro_name = read.csv(paste0(CDir,"COPDGene_SOMA1_3_MetaData_Sep18.csv"))

n_trimmed = dim(net)[2]
n_met = 7 #FEV1: 7 metabolites; emph: 10 metabolites
pro_label = AbarLabel[idx][(n_met+1):n_trimmed]
pro_name = pro_name %>% dplyr::select(pro_label)
pro_name = pro_name[3,]
met_name = AbarLabel[idx][1:n_met]
node_name = unlist(c(met_name, pro_name))




edgeCut  = EdgeCut #FEV1 adjusted: EdgeCut - alpha = 0.51, 0.015 - alpha = 0.25; FEV1 unadj: EdgeCut - alpha = 0.09, 0.65 - alpha = 0.048
subA = as.matrix(Abar[idx,idx]) #AbarSigned
subA = apply(subA, 1, function(x) (ifelse(x < edgeCut, 0, x)))



#id_rm = which(rowSums(subA) == 0)
#if (length(id_rm) > 0){
#subA = subA[-id_rm, -id_rm]
#net = net[,-id_rm]
#} else {subA = subA}



subA = subA + diag(nrow(subA))
test_g = graph_from_adjacency_matrix(subA, mode = "undirected", weighted = TRUE, diag = FALSE)
TOM_A = TOMsimilarity(subA, verbose = FALSE)

alpha_star = graph.density(test_g)
alpha_star


# regularized laplacian
L2 <- subA %>% 
  graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE,diag = FALSE) %>% 
  graph.laplacian(normalized = F) %>% as.matrix


#tau = round(alpha, 2) 
#L2 <- solve(diag(nrow(L2)) + tau*L2)

net = scale(net, center = TRUE, scale = TRUE)
Y = scale(Y, center = TRUE, scale = TRUE)




temp0 = prcomp(net)#, center = T, scale. = T)
temp0_var = pve(net, temp0$rotation[,pc_id])  #temp0$sdev[pc_id]^2/sum(temp0$sdev^2)
temp0_cor = cor(net %*% temp0$rotation[,pc_id], Y)




alpha_seq = c(0, alpha_star) #c(seq(0,1,0.25), alpha_star)
n_alpha = length(alpha_seq)
temp4_cor = temp4_var = c()


for (a in 1: n_alpha){
alpha = alpha_seq[a]

temp4 = (1-alpha)*(net %*% L2) + alpha*(net %*% TOM_A)
temp4 = prcomp(temp4)#, center = T, scale. = T)

#if(alpha == 0|alpha == 1){
#  temp4 = (1-alpha)*(net %*% L2) + alpha*(net %*% TOM_A)
#  temp4 = scale(temp4, center = TRUE, scale = TRUE)
#  temp4 = prcomp(temp4)#, center = T, scale. = T)
#}else{
#  temp4 =(1-alpha) * scale((net %*% L2), center = T, scale = T) + alpha * scale((net %*% TOM_A), center = T, scale = T)
#  temp4 = prcomp(temp4)}


temp4_var[a] =  pve(net, temp4$rotation[,pc_id]) #temp4$sdev[pc_id]^2/sum(temp4$sdev^2) 
temp4_cor[a] = cor(temp4$x[,pc_id], Y) 
}


df1 = data.frame(x = seq(1:n_trimmed), l1 = abs(temp0$rotation[,pc_id]), node_name = node_name)
ggplot(df1, aes(x = x, y = l1, label = node_name)) + geom_point(size = 0.1) + geom_text(label = ifelse(df1$l1 > 0.1, as.character(node_name), '')) + xlab("Node index") + #ylim(0,0.8)+
 ylab("Loadings") +  theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


df1 = data.frame(x = seq(1:n_trimmed), l1 = abs(temp4$rotation[,pc_id]), node_name = node_name)
ggplot(df1, aes(x = x, y = l1, label = node_name)) + geom_point(size = 0.1) + geom_text(label = ifelse(df1$l1 > 0.1, as.character(node_name), '')) + xlab("Node index") + ylim(0,0.8)+
 ylab("Loadings") +  theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#apply(temp4$rotation, 2, function(x) pve(net, x))


true_est = c(temp0_cor,  temp4_cor) #c(temp0_cor,  temp2_cor, temp3_cor, temp4_cor) 
true_var = c(temp0_var,  temp4_var) #c(temp0_var, temp2_var, temp3_var, temp4_var) 

true_est
true_var


# subsampling 

set.seed(8097)
n_all = c(500, 300, 200, 100, 50)#, 25, 10)
n_iter = 1000
n_version = n_alpha + 1
PC1_cor = PC1_exp = R2_all = array(0, dim = c(n_iter, n_version, length(n_all)))



  for (v in 1:length(n_all)){
      n1 = n_all[v]
    for (i in 1: n_iter){

    idx = sample(n, n1, replace = FALSE) 
    X_i = net[idx,]
    Y_i = Y[idx,]
    X_i = scale(X_i, center = T)
    Y_i = scale(Y_i, center = T)
    

  
    #observed
    pca0 = prcomp(X_i)#, center = T, scale. = T)  
    pca0_exp = pve(X_i, pca0$rotation[,pc_id])  #pca0$sdev[pc_id]^2/sum(pca0$sdev^2) 
    Z0_score = pca0$x[,pc_id] #X_i %*% pca0$rotation[,pc_id]
    #mod0 = lm(Y_i ~ Z0_score)  
    #R2_0 = summary(mod0)$adj.r.squared
    #R2_all[i,j,v] = R2_0
    PC1_exp[i,1,v] = pca0_exp
    PC1_cor[i,1,v] = cor(Y_i, Z0_score) 
    
    for (j in 2: (n_alpha+1)){
      alpha = alpha_seq[j-1]
      
      Z6 =  (1-alpha)* X_i %*% L2 + alpha* X_i %*% TOM_A
      pca6 = prcomp(Z6)#, center = T, scale. = T)
      
      #if(alpha == 0|alpha == 1){
      #  Z6 =  (1-alpha)* X_i %*% L2 + alpha* X_i %*% TOM_A
      #  pca6 = prcomp(Z6, center = T, scale. = T)
      #}else{
      #  Z6 = (1-alpha)*scale( X_i %*% L2, center = T, scale = T) + alpha*scale(X_i %*% TOM_A, center = T, scale = T)
      #  pca6 = prcomp(Z6)}
      
      
      pca6_exp = pve(X_i, pca6$rotation[,pc_id]) #pca6$sdev[pc_id]^2/sum(pca6$sdev^2)  
      Z6_score =  pca6$x[,pc_id]
      #mod6 = lm(Y_i ~ Z6_score) 
      #R2_6 = summary(mod6)$adj.r.squared
      #R2_all[i, j,v] = R2_6
      PC1_exp[i,j,v] = pca6_exp
      PC1_cor[i,j,v] = cor(Y_i, Z6_score) 
    }
  }
}

mu_cor = sd_cor = mu_pc1 = sd_pc1 = matrix(0, length(n_all), n_version)
for (v in 1:length(n_all)){
  mu_cor[v,] = round(apply(PC1_cor[,,v], 2, function(x) mean(abs(x))), 3)
  sd_cor[v,] = round(apply(PC1_cor[,,v], 2, function(x) sd(abs(x))), 3)
  mu_pc1[v,] = round(colMeans(PC1_exp[,,v]), 3)
  sd_pc1[v,] = round(apply(PC1_exp[,,v], 2, sd),3)
}



#----------------------------------------------------------------------------#
# mean corr plot
#----------------------------------------------------------------------------#
method_name = c()
for (a in 1:(n_alpha-1)){
  method_name[a] = paste0("alpha_", round(alpha_seq[a],2))
}

method_name = c(method_name, paste0("alpha_star"))




df = data.frame(Mean_corr = c(rbind(abs(true_est), mu_cor)), 
                sd_corr = c(rbind(c(rep(0,n_version)), sd_cor/sqrt(n_iter))),
                size = c(rep(c(n,n_all),n_version)),
                #method = c(rep("Xsub", length(n_all)+1), rep("X_L", length(n_all)+1), rep("X_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                #method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                method = c(rep("X only", length(n_all)+1), rep(method_name, each = length(n_all)+1)))

ggplot(df, aes(x = size, y = Mean_corr, ymin = Mean_corr - sd_corr, ymax = Mean_corr + sd_corr, colour = method)) + ylab("Mean Correlation")+  xlab("Sample size")+ ylim(0,0.5) +
  geom_point()+ geom_line() + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(axis.text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())# + theme(legend.position = "none")



#----------------------------------------------------------------------------#
# mean PC1 var plot 
#----------------------------------------------------------------------------#
df = data.frame(PC1_exp = c(rbind(abs(true_var), mu_pc1)), 
                sd_pc1 = c(rbind(c(rep(0,n_version)), sd_pc1/sqrt(n_iter))),
                size = c(rep(c(n,n_all),n_version)),
                #method = c(rep("Xsub", length(n_all)+1), rep("X_L", length(n_all)+1), rep("X_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                #method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                method = c(rep("X only", length(n_all)+1), rep(method_name, each = length(n_all)+1)))

ggplot(df, aes(x = size, y = PC1_exp, ymin = PC1_exp - sd_pc1, ymax = PC1_exp + sd_pc1, colour = method)) + ylab("Mean PC1 variance explained")+  xlab("Sample size") + ylim(0, 1) +
  geom_point()+ geom_line() + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(axis.text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


 #----------------------------------------------------------------------------#
# PC1 vs Corr
#----------------------------------------------------------------------------#
df2 = data.frame(pc = c(rbind(abs(true_var), mu_pc1 )), cor = c(rbind(abs(true_est), mu_cor)),
                 #method = c( rep("Xsub", length(n_all)+1), rep("X_L0_TOM", length(n_all)+1), rep("X_Lnorm_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)),
                 method = c( rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)),
                 size = c(rep(c(n,n_all), n_version)))


ggplot(df2, aes(x = pc, y = cor, colour = method)) + geom_point(size = 0.005*df2$size) + #ylim(0.25, 0.30) + #ylim(0.96, 1)+
  theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +xlab("PC1 Variance exp") + ylab("Correlation") 



plot(abs(temp4$rotation[,1]))
text(abs(temp4$rotation[,1]), labels = colnames(subA), cex = 0.6)


setwd("/Users/vut3/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation ")
save.image(file = "real_mp_FEV1unadj_a0048.RData") 


#----------------------------------------------------------------------# 
# PROTEIN NETWORKS
#----------------------------------------------------------------------# 

rm(list = ls())
library(dbplyr)
library(dplyr)
library(tidyverse)
library(WGCNA)
library(ggplot2)
library(igraph)
library(MASS)
library(lqmm)
library(pls)

pve = function(X, eigen_vec){
  num1 = X %*% eigen_vec
  deno1 = sum(X^2)
  res = sum(num1^2)/deno1
  return(res)
}

#----------------------------------------------------------------------# 
# Previous version
#----------------------------------------------------------------------# 
#setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation /Sim_RealData")
#load("FEV1_AA1a_adjusted.RData")   
#load("FEV1_NHW1a_adjusted.RData")

# new version (Dec 2021)
setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation /Sim_RealData/FEV1_NHW_Unadjusted")
load("FEV1_NHWnet16_1AA.Rdata") # scores and M
load("FEV1_NHW_Unadjusted_Data.Rdata") # X and Y
subA = M
test_g = graph_from_adjacency_matrix(subA, mode = "undirected", weighted = TRUE, diag = FALSE)



#----------------------------------------------------------------------# 
# FEV1
#----------------------------------------------------------------------# 
#setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation /ProteinNetwork_July/fev1_aa_unadjusted")
#load("FEV1_AA_Unadjusted_Data.Rdata")
#load("fev1_aa_unadjustednet55_1AA.Rdata")


#----------------------------------------------------------------------# 
# Smoking status
#----------------------------------------------------------------------# 
# AA
#setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation /ProteinNetwork_July/smoking_aa_unadjusted")
#load("Data.Rdata")
#load("Smoking_AAnet35_1AA.Rdata")

#NHW
#setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation /ProteinNetwork_July/smoking_nhw_unadjusted")
#load("Data.Rdata")
#load("Smoking_NHWnet35_1.Rdata")




test_g = graph_from_adjacency_matrix(subA, mode = "undirected", weighted = TRUE, diag = FALSE)
TOM_A = TOMsimilarity(subA, verbose = FALSE)

alpha = graph.density(test_g)
# regularized laplacian
L2 <- subA %>% 
  graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE,diag = FALSE) %>% 
  graph.laplacian(normalized = F) %>% as.matrix


#tau = round(alpha, 2) 
#L2 <- solve(diag(nrow(L2)) + tau*L2)

M_name = colnames(M)
X1_name = colnames(X1)
idx = sapply(M_name, function(x) which(X1_name == x))
idx = unlist(idx)
net = X1[,idx]
p = dim(net)[2]
n = dim(net)[1]
pc_id = 1


net = scale(net, center = TRUE, scale = TRUE)
Y = scale(Y, center = TRUE, scale = TRUE)




temp0 = prcomp(net, center = T, scale. = T)
temp0_var = pve(net, temp0$rotation[,pc_id]) #temp0$sdev[pc_id]^2/sum(temp0$sdev^2)
temp0_cor = cor(net %*% temp0$rotation[,pc_id], Y)



temp4  = (net %*% L2) #(1-alpha)*(net %*% L2) + alpha*(net %*% TOM_A)
temp4 = prcomp(temp4, center = T, scale. = T)
temp4_var = temp4$sdev[pc_id]^2/sum(temp4$sdev^2)   #pve(net, temp4$rotation[,pc_id]) 
temp4_cor = cor(temp4$x[,pc_id], Y) #cor(net %*% temp4$rotation[,pc_id], Y) 

plot(abs(temp4$rotation[,pc_id]))
text(abs(temp4$rotation[,pc_id]), labels = colnames(subA))



#apply(temp4$rotation, 2, function(x) pve(net, x))


true_est = c(temp0_cor,  temp4_cor)
true_var = c(temp0_var,  temp4_var)

true_est
true_var


# subsampling 
set.seed(8097)
n_all = c(1000, 500, 300, 200, 100, 50)
n_iter = 1000
n_version = 2
PC1_cor = PC1_exp = R2_all = array(0, dim = c(n_iter, n_version, length(n_all)))

for (v in 1:length(n_all)){
  n1 = n_all[v]
  for (i in 1: n_iter){
    j = 1
    idx = sample(n, n1, replace = FALSE) 
    
    X_i = net[idx,]
    Y_i = Y[idx,]
    X_i = scale(X_i, center = T)#, scale = T)
    Y_i = scale(Y_i, center = T)
    
    
    
    #observed
    pca0 = prcomp(X_i, center = T, scale. = T)  
    pca0_exp = pca0$sdev[pc_id]^2/sum(pca0$sdev^2) #pve(X_i, pca0$rotation[,pc_id]) 
    Z0_score = pca0$x[,pc_id] #X_i %*% pca0$rotation[,pc_id]
    mod0 = lm(Y_i ~ Z0_score)  
    R2_0 = summary(mod0)$adj.r.squared
    R2_all[i,j,v] = R2_0
    PC1_exp[i,j,v] = pca0_exp[pc_id]
    PC1_cor[i,j,v] = cor(Y_i, Z0_score) 
    
    j = j+1
    
    Z6 =  X_i %*% L2  #(1-alpha)* X_i %*% L2 + alpha* X_i %*% TOM_A
    pca6 = prcomp(Z6, center = T, scale. = T) 
    pca6_exp = pca6$sdev[pc_id]^2/sum(pca6$sdev^2) #pve(X_i, pca6$rotation[,pc_id])  
    Z6_score =  pca6$x[,pc_id] #X_i %*% pca6$rotation[,pc_id] 
    mod6 = lm(Y_i ~ Z6_score) # lm(Y0 ~ Z6_score)
    R2_6 = summary(mod6)$adj.r.squared
    
    R2_all[i, j,v] = R2_6
    PC1_exp[i,j,v] = pca6_exp[pc_id]
    PC1_cor[i,j,v] = cor(Y_i, Z6_score) 
    
  }
}

mu_cor = sd_cor = mu_pc1 = sd_pc1 = matrix(0, length(n_all), n_version)
for (v in 1:length(n_all)){
  mu_cor[v,] = round(apply(PC1_cor[,,v], 2, function(x) mean(abs(x))), 3)
  sd_cor[v,] = round(apply(PC1_cor[,,v], 2, function(x) sd(abs(x))), 3)
  mu_pc1[v,] = round(colMeans(PC1_exp[,,v]), 3)
  sd_pc1[v,] = round(apply(PC1_exp[,,v], 2, sd),3)
}



#----------------------------------------------------------------------------#
# mean corr plot
#----------------------------------------------------------------------------#
df = data.frame(Mean_corr = c(rbind(abs(true_est), mu_cor)), 
                sd_corr = c(rbind(c(rep(0,n_version)), sd_cor/sqrt(n_iter))),
                size = c(rep(c(n,n_all),n_version)),
                #method = c(rep("Xsub", length(n_all)+1), rep("X_L0_TOM", length(n_all)+1), rep("X_Lnorm_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))

ggplot(df, aes(x = size, y = Mean_corr, ymin = Mean_corr - sd_corr, ymax = Mean_corr + sd_corr, colour = method)) + ylab("Mean Correlation")+ ylim(0.25,0.50) +
  geom_point()+ geom_line() + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 



#----------------------------------------------------------------------------#
# mean PC1 var plot 
#----------------------------------------------------------------------------#
df = data.frame(PC1_exp = c(rbind(abs(true_var), mu_pc1)), 
                sd_pc1 = c(rbind(c(rep(0,n_version)), sd_pc1/sqrt(n_iter))),
                size = c(rep(c(n,n_all),n_version)),
                #method = c(rep("Xsub", length(n_all)+1), rep("X_L0_TOM", length(n_all)+1), rep("X_Lnorm_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))

ggplot(df, aes(x = size, y = PC1_exp, ymin = PC1_exp - sd_pc1, ymax = PC1_exp + sd_pc1, colour = method)) + ylab("Mean PC1 variance explained")+ #ylim(0.15, 0.40) +
  geom_point()+ geom_line() + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


#----------------------------------------------------------------------------#
# PC1 vs Corr
#----------------------------------------------------------------------------#
df2 = data.frame(pc = c(rbind(abs(true_var), mu_pc1 )), cor = c(rbind(abs(true_est), mu_cor)),
                 #method = c( rep("Xsub", length(n_all)+1), rep("X_L0_TOM", length(n_all)+1), rep("X_Lnorm_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)),
                 method = c( rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)),
                 size = c(rep(c(n,n_all), n_version)))


ggplot(df2, aes(x = pc, y = cor, colour = method)) + geom_point(size = 0.005*df2$size) + ylim(0.25, 0.50) + #ylim(0.96, 1)+
  theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +xlab("PC1 Variance exp") + ylab("Correlation") 



plot(abs(temp4$rotation[,1]))
text(abs(temp4$rotation[,1]), labels = colnames(subA), cex = 0.6)


E(test_g)$width = E(test_g)$weight*5
V(test_g)$size = degree(test_g)
coords = layout_in_circle(test_g)
plot(test_g, edge.color = "black", layout = coords, vertex.label.cex = 0.75)



#----------------------------------------------------------------------# 
# Smoking status
#----------------------------------------------------------------------# 
rm(list = ls())
library(dbplyr)
library(dplyr)
library(tidyverse)
library(WGCNA)
library(ggplot2)
library(igraph)
library(MASS)
library(lqmm)
library(pls)

pve = function(X, eigen_vec){
  num1 = X %*% eigen_vec
  deno1 = sum(X^2)
  res = sum(num1^2)/deno1
  return(res)
}

#----------------------------------------------------------------------# 
# Smoking status
#----------------------------------------------------------------------# 
# AA
setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation /ProteinNetwork_July/smoking_aa_unadjusted")
load("Data.Rdata")
load("Smoking_AAnet35_1AA.Rdata")

#NHW
setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation /ProteinNetwork_July/smoking_nhw_unadjusted")
load("Data.Rdata")
load("Smoking_NHWnet35_1.Rdata")


M_name = colnames(M)
X1_name = colnames(X1)
idx = sapply(M_name, function(x) which(X1_name == x))
idx = unlist(idx)
net = X1[,idx]
p = dim(net)[2]
n = dim(net)[1]
pc_id = 1

test_g = graph_from_adjacency_matrix(M, mode = "undirected", weighted = TRUE, diag = FALSE)
alpha = graph.density(test_g)
#plot(test_g)
TOM_A = TOMsimilarity(M, verbose = FALSE)



net = scale(net, center = TRUE, scale = TRUE)
#Y = scale(Y, center = TRUE, scale = TRUE)




temp0 = prcomp(net)
temp0_var = temp0$sdev[pc_id]^2/sum(temp0$sdev^2)
temp0_ld = temp0$rotation[,pc_id]
plot(abs(temp0_ld))
text(abs(temp0_ld), labels = colnames(M), cex = 0.5)
org_score = temp0$x[,pc_id]



L2 <- M %>% 
  graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE,diag = FALSE) %>% 
  graph.laplacian(normalized = F) %>% as.matrix


temp4  = (1-alpha)*(net %*% L2) + alpha*(net %*% TOM_A)
temp4 = prcomp(temp4)#, center = T, scale. = T)
temp4_var = temp4$sdev[pc_id]^2/sum(temp4$sdev^2)
temp4_ld = temp4$rotation[,pc_id]
new_score = temp4$x[,pc_id]

plot(abs(temp4_ld))
text(abs(temp4$rotation[,pc_id]), labels = colnames(M), cex = 0.5)

plot(abs(loading$score))
text(abs(loading$score), labels = colnames(M), cex = 0.5)




#==========================================================##==========================================================#
# miRNA networks from Yong
# Reduced sample size
#==========================================================##==========================================================#

rm(list = ls())
library(dbplyr)
library(dplyr)
library(tidyverse)
library(WGCNA)
library(ggplot2)
library(igraph)
library(MASS)
library(lqmm)
library(mvtnorm)
library(igraph)


Danaher_pos_def <- function(m, cc = 3.4, dd = 0.95){
  AA <- m*cc + t(m*cc)
  AA[AA>0] <- 1
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA <- AA/( as.matrix(rep(1,nrow(m))) ) %*% t(1.4*rowSums(abs(AA)))
  AA <- (AA+t(AA))/2
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA
}

pos_def_func = function(X){
  e_val = eigen(X)$values
  e_vec = eigen(X)$vectors
  neg_e = e_val[e_val < 0]
  if (length(neg_e) > 0) {e_val = e_val + max(abs(neg_e)) + 0.1 }
  else e_val = e_val
  Xnew = e_vec %*% diag(e_val) %*% t(e_vec)#solve(e_vec)
  return(Xnew)
}



pve = function(X, eigen_vec){
  num1 = X %*% eigen_vec
  deno1 = sum(X^2)
  res = sum(num1^2)/deno1
  return(res)
}


setwd("~/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation ")
load("emphysema_modules.RData")
#load("FEV1pp_modules.RData")

pc_id = 1
# pctEmph
idx = mirGeneModule_emphy_module3[[1]] #mirGeneModule_emphy_module2[[1]]
# FEV1pp
#idx = mirGeneModule_FEV1pp_module3[[1]] 

X1 = X[, idx]
Y1 = Y
n = nrow(X1)
p = length(idx)


subA = as.matrix(Abar[idx,idx]) 
# distribution of edge weight

#ttt = data.frame(edge_weight = as.numeric(subA[lower.tri(subA, diag = FALSE)]))
#ggplot(ttt, aes(edge_weight)) + geom_histogram() 

#subA = apply(subA, 1, function(x) (ifelse(x > 0, runif(1),  x)))
thres = 0.00024 # module 1: thres = 0.0095 to get alpha = 0.5, 0.0112 to get alpha = 0.25; Module 3: thres = 0.00019 alpha = 0.52, thres = 0.00024 alpha = 0.248
subA = apply(subA, 1, function(x) (ifelse(x < thres, 0, x)))

#diag(subA) = 0
#D = rowSums(as.matrix(subA))
#id_rm = which(D == 0)

#if (length(id_rm) > 0){
#subA = subA[-id_rm, -id_rm]
#X1 = X1[,-id_rm]
#A_name = AbarLabel[idx][-id_rm]
#} else {subA = as.matrix(Abar[idx,idx])
#        X1 = X[,idx]
#        A_name = AbarLabel[idx]}

#p1 = dim(X1)[2]



#id_rm = which(rowSums(subA) == 0)
#if (length(id_rm) > 0){
#  subA = subA[-id_rm, -id_rm]
#  X1 = X1[,-id_rm]
#} else {subA = subA}



subA = subA + diag(nrow(subA))
test_g = graph_from_adjacency_matrix(subA, mode = "undirected", weighted = TRUE, diag = FALSE)
TOM_A = TOMsimilarity(subA, verbose = FALSE)

alpha_star = graph.density(test_g)
alpha_star



# regularized laplacian
L2 <- subA %>% 
  graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE,diag = FALSE) %>% 
  graph.laplacian(normalized = F) %>% as.matrix


#tau = round(alpha, 2) 
#L2 <- solve(diag(nrow(L2)) + tau*L2)

net = scale(X1, center = TRUE, scale = TRUE)
Y1 =  scale(Y1, center = TRUE, scale = TRUE)




temp0 = prcomp(net, center = T, scale. = T)
temp0_var = pve(net, temp0$rotation[,pc_id])  #temp0$sdev[pc_id]^2/sum(temp0$sdev^2) 
temp0_cor = cor(net %*% temp0$rotation[,pc_id], Y1)




alpha_seq = c(0, alpha_star) #c(seq(0,1,0.25), alpha_star)
n_alpha = length(alpha_seq)
temp4_cor = temp4_var = c()


for (a in 1: n_alpha){
  alpha = alpha_seq[a]
 
  #temp4 = (1-alpha)*(net %*% L2) + alpha*(net %*% TOM_A)
  #temp4 = prcomp(temp4)#, center = T, scale. = T)
  
  if(alpha == 0|alpha == 1){
    temp4 = (1-alpha)*(net %*% L2) + alpha*(net %*% TOM_A)
    temp4 = prcomp(temp4, center = T, scale. = T)
  }else{
    temp4 = scale((1-alpha)*(net %*% L2), center = T, scale = T) + scale(alpha*(net %*% TOM_A), center = T, scale = T)
    temp4 = prcomp(temp4)}
  

  temp4_var[a] = pve(net, temp4$rotation[,pc_id])  #temp4$sdev[pc_id]^2/sum(temp4$sdev^2) 
  temp4_cor[a] = cor(temp4$x[,pc_id], Y1) 
}


true_est =  c(temp0_cor,  temp4_cor) #c(temp0_cor,  temp2_cor, temp3_cor, temp4_cor) 
true_var =  c(temp0_var,  temp4_var) #c(temp0_var, temp2_var, temp3_var, temp4_var)   

#true_est
#true_var



#df1 = data.frame(x = c(x = 1:p),l1 = abs(temp0$rotation[,pc_id]), node_name = AbarLabel[idx])
#ggplot(df1, aes(x = x, y = l1, label = node_name)) + geom_point(size = 0.1)  + xlab("Node index") +  geom_text(label = ifelse(df1$l1 > 0.1, as.character(df1$node_name), '')) + #ylim(0,0.2) +
#  ylab("Loadings") +  theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


#df1 = data.frame(x = c(x= 1:length(idx)),l1 = abs(temp2$rotation[,pc_id]), node_name = AbarLabel[idx])
#ggplot(df1, aes(x = x, y = l1, label = node_name)) + geom_point(size = 0.1)  + xlab("Node index") +  geom_text(label = ifelse(df1$l1 > 0.1, as.character(df1$node_name), '')) + #ylim(0,0.7) +
#  ylab("Loadings") +  theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


#df1 = data.frame(x = c(x = 1:p1),l1 = abs(temp4$rotation[,pc_id]), node_name = A_name)
#ggplot(df1, aes(x = x, y = l1, label = node_name)) + geom_point(size = 0.1)  + xlab("Node index") +  geom_text(label = ifelse(df1$l1 > 0.1, as.character(df1$node_name), '')) + #ylim(0,0.2) +
#  ylab("Loadings") +  theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# subsampling 
set.seed(8097)
n_all = c(350, 300, 200, 100, 50)
n_iter = 1000
n_version = n_alpha + 1
PC1_cor = PC1_exp = R2_all = array(0, dim = c(n_iter, n_version, length(n_all)))

for (v in 1:length(n_all)){
  n1 = n_all[v]
  for (i in 1: n_iter){
    j = 1
    idx = sample(n, n1, replace = FALSE) 
    
    X_i = net[idx,]
    Y_i = Y1[idx,]
    X_i = scale(X_i, center = T)#, scale = T)
    Y_i = scale(Y_i, center = T)
    
    
    
    #observed
    pca0 = prcomp(X_i, center = T, scale. = T)  
    pca0_exp = pve(X_i, pca0$rotation[,pc_id])  #pca0$sdev[pc_id]^2/sum(pca0$sdev^2) 
    Z0_score = pca0$x[,pc_id] #X_i %*% pca0$rotation[,pc_id]
    #mod0 = lm(Y_i ~ Z0_score)  
    #R2_0 = summary(mod0)$adj.r.squared
    #R2_all[i,j,v] = R2_0
    PC1_exp[i,1,v] = pca0_exp[pc_id]
    PC1_cor[i,1,v] = cor(Y_i, Z0_score) 
    
  
    for (j in 2: (n_alpha+1)){
      alpha = alpha_seq[j-1]
      
      #Z6 =  (1-alpha)* X_i %*% L2 + alpha* X_i %*% TOM_A
      #pca6 = prcomp(Z6)#, center = T, scale. = T)  
      
      if(alpha == 0|alpha == 1){
        Z6 =  (1-alpha)* X_i %*% L2 + alpha* X_i %*% TOM_A
        pca6 = prcomp(Z6, center = T, scale. = T)
      }else{
        Z6 = scale((1-alpha)* X_i %*% L2, center = T, scale = T) + scale( alpha* X_i %*% TOM_A, center = T, scale = T)
        pca6 = prcomp(Z6)}
      

      
      pca6_exp = pve(X_i, pca6$rotation[,pc_id]) #pca6$sdev[pc_id]^2/sum(pca6$sdev^2)  
      Z6_score =  pca6$x[,pc_id]
      #mod6 = lm(Y_i ~ Z6_score) 
      #R2_6 = summary(mod6)$adj.r.squared
      
      #R2_all[i, j,v] = R2_6
      PC1_exp[i,j,v] = pca6_exp
      PC1_cor[i,j,v] = cor(Y_i, Z6_score) 
      
    }
    
  }
}


mu_cor = sd_cor = mu_pc1 = sd_pc1 = matrix(0, length(n_all), n_version)
for (v in 1:length(n_all)){
  mu_cor[v,] = round(apply(PC1_cor[,,v], 2, function(x) mean(abs(x))), 3)
  sd_cor[v,] = round(apply(PC1_cor[,,v], 2, function(x) sd(abs(x))), 3)
  mu_pc1[v,] = round(colMeans(PC1_exp[,,v]), 3)
  sd_pc1[v,] = round(apply(PC1_exp[,,v], 2, sd),3)
}




method_name = c()
for (a in 1:(n_alpha-1)){
  method_name[a] = paste0("alpha_", round(alpha_seq[a],2))
}

method_name = c(method_name, paste0("alpha_star"))

# CORRELATION PLOT
# MULTIPLE ALPHA
df = data.frame(Mean_corr = c(rbind(abs(true_est), mu_cor)), 
                sd_corr = c(rbind(c(rep(0,n_version)), sd_cor/sqrt(n_iter))),
                size = c(rep(c(n,n_all),n_version)),
                #method = c(rep("Xsub", length(n_all)+1), rep("X_L", length(n_all)+1), rep("X_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                #method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                method = c(rep("X only", length(n_all)+1), rep(method_name, each = length(n_all)+1)))

ggplot(df, aes(x = size, y = Mean_corr, ymin = Mean_corr - sd_corr, ymax = Mean_corr + sd_corr, colour = method)) + ylab("Mean Correlation")+  xlab("Sample size")+ ylim(0, 0.5) +
  geom_point()+ geom_line() + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(axis.text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")




#ONE ALPHA
#df = data.frame(Mean_corr = c(rbind(abs(true_est), mu_cor)), 
#                sd_corr = c(rbind(c(rep(0,n_version)), sd_cor/sqrt(n_iter))),
#                size = c(rep(c(n,n_all),n_version)),
                #method = c(rep("Xsub", length(n_all)+1), rep("X_L", length(n_all)+1), rep("X_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                #method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
 #               method = c(rep("Xsub", length(n_all)+1), rep(method_name, each = length(n_all)+1)))

#ggplot(df, aes(x = size, y = Mean_corr, ymin = Mean_corr - sd_corr, ymax = Mean_corr + sd_corr, colour = method)) + ylab("Mean Correlation")+  xlab("Sample size")+ ylim(0,0.40) +
#  geom_point()+ geom_line() + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
#  theme(axis.text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


#----------------------------------------------------------------------------#
# mean PC1 var plot 
#----------------------------------------------------------------------------#
df = data.frame(PC1_exp = c(rbind(abs(true_var), mu_pc1)), 
                sd_pc1 = c(rbind(c(rep(0,n_version)), sd_pc1/sqrt(n_iter))),
                size = c(rep(c(n,n_all),n_version)),
                #method = c(rep("Xsub", length(n_all)+1), rep("X_L0_TOM", length(n_all)+1), rep("X_Lnorm_TOM", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                #method = c(rep("Xsub", length(n_all)+1), rep("X_weighted", length(n_all)+1)))
                method = c(rep("X only", length(n_all)+1), rep(method_name, each = length(n_all)+1)))

ggplot(df, aes(x = size, y = PC1_exp, ymin = PC1_exp - sd_pc1, ymax = PC1_exp + sd_pc1, colour = method)) + ylab("Mean PC1 variance explained")+  xlab("Sample size")+ ylim(0, 1) +
  geom_point()+ geom_line() + geom_errorbar() + scale_x_reverse()+ theme(plot.title = element_text(hjust = 0.5))+ theme_gray() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  theme(legend.position = "none")


setwd("/Users/vut3/Documents/Multiomic project/Subnetwork summarization/subnetwork simulation ")
save.image(file = "real_mm_p86_module3_a1_ver2.RData") #p=60, prob = 0.65 for alpha_0 = 0.6
#save.image(file = "real_mm_p86_module3_a025_new.RData")

#E(test_g)$width = E(test_g)$weight*5
#V(test_g)$size = degree(test_g)
#coords = layout_in_circle(test_g)
#plot(test_g, edge.color = "black", layout = coords, vertex.label.cex = 0.75)



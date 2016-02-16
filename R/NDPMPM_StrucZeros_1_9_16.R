#Fit the model
#library(gtools)
library(DirichletReg)
library(matrixStats)
library(coda)

X_house = read.csv("X_house.csv",header=TRUE)
X_indiv = read.csv("X_indiv.csv",header=TRUE)
pi_true = read.csv("pi_true.csv",header=TRUE)
omega_true = read.csv("omega_true.csv",header=TRUE)
lambda_true = read.csv("lambda_true.csv",header=TRUE)
phi_true = read.csv("phi_true.csv",header=TRUE)
G_true = read.csv("G.csv",header=TRUE)
M_true = read.csv("M.csv",header=TRUE)
Data_house = apply(X_house,2,as.factor)
Data_house = data.frame(Data_house)
Data_indiv = apply(X_indiv,2,as.factor)
Data_indiv = data.frame(Data_indiv)
N = dim(Data_indiv)[1]
n = dim(Data_house)[1]
n_i = Data_house[,1]
p = dim(Data_indiv)[2]
q = dim(Data_house)[2]

#Initialize chain
d_k_house = d_k_indiv = ini_marg_house = ini_marg_indiv = NULL
for(k in 1:q){
  d_k_house = cbind(d_k_house,nlevels(Data_house[,k]))
  ini_marg_k = as.data.frame(table(Data_house[,k]))$Freq/n
  ini_marg_k = matrix(ini_marg_k,ncol=1)
  ini_marg_house = rbind(ini_marg_house,ini_marg_k)
}

for(k in 1:p){
  d_k_indiv = cbind(d_k_indiv,nlevels(Data_indiv[,k]))
  ini_marg_k = as.data.frame(table(Data_indiv[,k]))$Freq/N
  ini_marg_k = matrix(ini_marg_k,ncol=1)
  ini_marg_indiv = rbind(ini_marg_indiv,ini_marg_k)
}

FF = 20
SS = 15
alpha = beta = 1
a_kdk = 1
a_alpha = b_alpha = a_beta = b_beta = 0.25
pr_G_post = matrix(0,ncol=FF,nrow=n)
pii = matrix(0,nrow=FF)
lambda = matrix(rep(ini_marg_house,FF),ncol=FF)
pr_M_post = matrix(0,ncol=SS,nrow=N)
omega = matrix(0,nrow=FF,ncol=SS)
phi = matrix(0,nrow=length(ini_marg_indiv),ncol=FF*SS) #make phi matrix and not array for c++
for(gm in 1:(FF*SS)){
    phi[,gm] = ini_marg_indiv
}
U = matrix(rbeta(FF,1,alpha),nrow=FF)
V = matrix(rbeta((FF*SS),1,beta),nrow=FF,ncol=SS)
U[FF]=1
V[,SS]=1
one_min_U = 1L-U
one_min_U = c(1,cumprod(one_min_U[1:(FF-1)]))
one_min_V = 1L-V
one_min_V = cbind(1,t(apply(one_min_V[,-SS],1,cumprod)))
pii = U*one_min_U
omega = V*one_min_V
n_iter = 10000
burn_in = n_iter/2
dp_imput_house = dp_imput_indiv = dp_imput = NULL
MM = 10
M_to_use = seq((burn_in+1), n_iter, (length(c((burn_in+1):n_iter))/(5*MM)))
M_to_use_mc = sample(M_to_use,MM,replace=FALSE)
FFF_indiv = matrix(rep(cumsum(c(0,d_k_indiv[,-p])),each=N),ncol=p)
FFF_house = matrix(rep(cumsum(c(0,d_k_house[,-q])),each=n),ncol=q)
house_index = rep(c(1:n),n_i)
ALPHA = BETA = PII = G_CLUST = M_CLUST = N0_STAR_02 = NULL
LAMBDA = matrix(0,ncol=ncol(lambda),nrow=nrow(lambda))
OMEGA = matrix(0,ncol=ncol(omega),nrow=nrow(omega))
cn_i = c(n_i)
n_i_index = rep(n_i,n_i)
d_k_gm = 1+cumsum(c(0,d_k_indiv[,-p]))
d_k = 1+cumsum(c(0,d_k_house[,-q]))

#List all possible combinations for households of size 2
All.poss.combs_02 = as.matrix(expand.grid(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2)))

#Now check structural zeros within combinations -- households of size 2
FFF_indiv_02 = matrix(rep(cumsum(c(0,d_k_indiv[,-p])),each=2*nrow(All.poss.combs_02)),ncol=p)
struc_zero_chk_counter_02 = matrix(0,nrow=nrow(All.poss.combs_02),ncol=1)
Data_indiv_02 = NULL
cn_i_02 = c(rep(2,nrow(All.poss.combs_02)))
for(str_ch_02 in 1:nrow(All.poss.combs_02)){
  hh_to_check_02 = matrix(All.poss.combs_02[str_ch_02,],byrow=T,ncol=p)
  hh_to_check_02 = hh_to_check_02[order(hh_to_check_02[,3]),]
  if(hh_to_check_02[1,1] == 1 & hh_to_check_02[1,2] == 2 |
     hh_to_check_02[1,1] == 2 & hh_to_check_02[1,3] == 1 |
     hh_to_check_02[2,1] == 1 & hh_to_check_02[2,2] == 2 |
     hh_to_check_02[2,1] == 2 & hh_to_check_02[2,3] == 1 |

     hh_to_check_02[1,3] == 1 & hh_to_check_02[2,3] == 1 |
     hh_to_check_02[1,1] == 2 & hh_to_check_02[1,3] == 1 &
     hh_to_check_02[2,1] == 2 & hh_to_check_02[2,3] == 2 |
     hh_to_check_02[1,2] == 1 & hh_to_check_02[1,3] == 1 &
     hh_to_check_02[2,2] == 2 & hh_to_check_02[2,3] == 2 |
     hh_to_check_02[1,1] == 1 & hh_to_check_02[1,2] == 2 &
     hh_to_check_02[2,1] == 2 & hh_to_check_02[2,2] == 1 ){

    struc_zero_chk_counter_02[str_ch_02] = 0
  } else{struc_zero_chk_counter_02[str_ch_02] = 1 }
  Data_indiv_02 = rbind(Data_indiv_02,hh_to_check_02)
}
phi_index_02 = data.matrix(Data_indiv_02)+FFF_indiv_02
n_star_02 = length(unique(house_index[which(n_i_index==2)]))
combs_min_hhsize_02 = as.matrix(expand.grid(c(1,2),c(1,2))) #combs for other hh variables
lambda_index_02 = matrix(rep(cumsum(d_k_house[,-q]),each=nrow(combs_min_hhsize_02)),ncol=q-1)+
  combs_min_hhsize_02
lambda_index_02 = cbind(2,lambda_index_02) #where 2 is the household size


ptm = proc.time()
for(mc in 1:n_iter){
  #structural zeros
  pr_each_comb_02 = prEachComb02Func(phi_index_02,lambda_index_02,phi,lambda,omega,c(pii),FF,SS,cn_i_02)
  pr_each_comb_02 = pr_each_comb_02/sum(pr_each_comb_02) #renormalize before N_Mult sampling
  pr_struc_zero_02 = pr_each_comb_02[which(struc_zero_chk_counter_02==0)]
  struc_zero_to_use_02 = All.poss.combs_02[which(struc_zero_chk_counter_02==0),]
  n0_star_02 = rnbinom(1,n_star_02,prob=(1-sum(pr_struc_zero_02)))
  n1.to.nC_02 = rmultinom(1, n0_star_02, prob=pr_struc_zero_02)
  pr_G_post_02 = prGpost02Func(phi_index_02,lambda_index_02,phi,lambda,omega,c(pii),FF,SS,cn_i_02)
  pr_G_post_02 = pr_G_post_02[which(struc_zero_chk_counter_02==0),]
  rep_n1.to.nC_02 = rep(n1.to.nC_02,ncol(pr_G_post_02))
  rep_pr_G_post_02 = matrix(rep(c(pr_G_post_02),rep_n1.to.nC_02),ncol=ncol(pr_G_post_02))
  Ran_unif_G_02 = runif(nrow(rep_pr_G_post_02))
  cumul_G_02 = rep_pr_G_post_02%*%upper.tri(diag(ncol(rep_pr_G_post_02)),diag=TRUE)
  G_02 = rowSums(Ran_unif_G_02>cumul_G_02) + 1L
  rep_G_02 = rep(G_02,each=2)
  Data_indiv_struc_02 = matrix(rep(c(struc_zero_to_use_02),
             rep(n1.to.nC_02,ncol(struc_zero_to_use_02))),ncol=ncol(struc_zero_to_use_02))
  Data_indiv_struc_02 = t(matrix(c(t(Data_indiv_struc_02)),nrow=p))
  FFF_indiv_struc_02 = matrix(rep(cumsum(c(0,d_k_indiv[,-p])),each=nrow(Data_indiv_struc_02)),ncol=p)
  phi_index_struc_02 = data.matrix(Data_indiv_struc_02)+FFF_indiv_struc_02
  pr_M_post_02 = prMpostFunc(phi_index_struc_02,phi,omega,rep_G_02,FF,SS)
  Ran_unif_M_02 = runif(nrow(pr_M_post_02))
  cumul_M_02 = pr_M_post_02%*%upper.tri(diag(ncol(pr_M_post_02)),diag=TRUE)
  M_02 = rowSums(Ran_unif_M_02>cumul_M_02) + 1L
  Data_indiv_struc_02 = apply(Data_indiv_struc_02,2,as.factor)
  Data_indiv_struc_02 = data.frame(Data_indiv_struc_02)
  Data_house_struc_02 = matrix(0,nrow=length(G_02),ncol=q)
  Data_house_struc_02[,1] = 2
  lambda_g_02 = t(lambda[,G_02])
  for(kjk in 2:q){
    pr_X_q_02 = lambda_g_02[,d_k[kjk]:cumsum(d_k_house)[kjk]]
    Ran_unif_q_02 = runif(nrow(pr_X_q_02))
    cumul_q_02 = pr_X_q_02%*%upper.tri(diag(ncol(pr_X_q_02)),diag=TRUE)
    Data_house_struc_02[,kjk] = rowSums(Ran_unif_q_02>cumul_q_02) + 1L
  }
  Data_house_struc_02 = apply(Data_house_struc_02,2,as.factor)
  Data_house_struc_02 = data.frame(Data_house_struc_02)
  Data_house_struc_02[,1] = factor(Data_house_struc_02[,1],levels=c(1:d_k_house[1]))

  #sample G
  phi_index = data.matrix(Data_indiv)+FFF_indiv
  lambda_index = data.matrix(Data_house)+FFF_house
  pr_G_post = prGpostFunc(phi_index,lambda_index,phi,lambda,omega,c(pii),FF,SS,cn_i)
  Ran_unif_G = runif(nrow(pr_G_post))
  cumul_G = pr_G_post%*%upper.tri(diag(ncol(pr_G_post)),diag=TRUE)
  G = rowSums(Ran_unif_G>cumul_G) + 1L

  #sample M
  rep_G = rep(G,cn_i)
  pr_M_post = prMpostFunc(phi_index,phi,omega,rep_G,FF,SS)
  Ran_unif_M = runif(nrow(pr_M_post))
  cumul_M = pr_M_post%*%upper.tri(diag(ncol(pr_M_post)),diag=TRUE)
  M = rowSums(Ran_unif_M>cumul_M) + 1L
  #M = as.matrix(M_true)

  #sample phi
  for(gg in 1:SS){
    for(ggg in 1:p){
      phi[d_k_gm[ggg]:cumsum(d_k_indiv)[ggg],(c(1:FF)+((gg-1)*FF))] =
        t(rdirichlet(FF,matrix(a_kdk + table(factor(rep_G[which(M==gg)],levels=c(1:FF)),
           Data_indiv[which(M==gg),ggg])+
             table(factor(rep_G_02[which(M_02==gg)],levels=c(1:FF)),
              Data_indiv_struc_02[which(M_02==gg),ggg]),nrow=FF)))
    }
  }

  #sample lambda
  for(kk in 1:q){
    lambda[d_k[kk]:cumsum(d_k_house)[kk],] =
      t(rdirichlet(FF,matrix(a_kdk + table(factor(G,levels=c(1:FF)),Data_house[,kk])
             +table(factor(G_02,levels=c(1:FF)),Data_house_struc_02[,kk]),nrow=FF)))
  }

  #sample U and pii
  n_f = matrix(summary(factor(G,levels=c(1:FF))),ncol=1)+
    matrix(summary(factor(G_02,levels=c(1:FF))),ncol=1)
  U[FF]=1
  U[1:(FF-1),1] = rbeta((FF-1),(1L+n_f[1:(FF-1)]),(alpha+(sum(n_f)-cumsum(n_f[-FF]))))
  if(length(which(U[-FF]==1))>0){
    U[which(U[-FF]==1)] = 0.99999
  }
  one_min_U = 1L-U
  one_min_U_prod = c(1,cumprod(one_min_U[1:(FF-1)]))
  pii = U*one_min_U_prod
  #if(length(which(pii==0))>0){
  #  pii[which(pii==0)] = 0.000001
  #}

  #sample V and omega
  M_G = table(factor(rep_G,levels=c(1:FF)),factor(M,levels=c(1:SS))) +
    table(factor(rep_G_02,levels=c(1:FF)),factor(M_02,levels=c(1:SS)))
  n_gm = as.data.frame(M_G)$Freq
  V[,SS]=1
  no_V_to_sim = (FF*(SS-1))
  V[,1:(SS-1)] = rbeta(no_V_to_sim,(1L+n_gm[1:no_V_to_sim]),
                       (beta + c( matrix(rowSums(M_G),ncol=SS-1,nrow=FF)-
                                    t(apply(M_G,1,cumsum))[,1:SS-1])))
  if(length(which(V[,-SS]==1))>0){
    V[which(V[,-SS]==1)] = 0.99999
  }
  one_min_V = 1L-V
  one_min_V_prod = cbind(1,t(apply(one_min_V[,-SS],1,cumprod)))
  omega = V*one_min_V_prod
  #if(length(which(omega==0))>0){
  #  omega[which(omega==0)] = 0.000001
  #}

  #sample alpha
  alpha = rgamma(1,shape=(a_alpha+FF-1),rate=(b_alpha-log(pii[FF])))

  #sample beta
  beta = rgamma(1,shape=(a_beta+(FF*(SS-1))),rate=(b_beta-sum(log(omega[,SS]))))

  #check number of occupied clusters
  S.occup = NULL
  for(occ in sort(unique(G))){
    S.occup = rbind(S.occup,dim(table(rep_G[which(rep_G==occ)],M[which(rep_G==occ)]))[2])
  }

  #print
  cat("Iter=", mc,"\t"," ",
      "F=",formatC(length(unique(G)), width=2, flag=" "),"\t",
      "S=",formatC(max(S.occup), width=2, flag=" "),"\t",
      "alp=",alpha,"\t",
      "bet=",beta,"\t",
      "n_02=",n0_star_02,"\n",sep = " ")

  #save and sample synthetic data
  if(mc > burn_in){
    PII = rbind(PII,c(pii))
    ALPHA = rbind(ALPHA,alpha)
    G_CLUST = rbind(G_CLUST,length(unique(G)))
    M_CLUST = rbind(G_CLUST,max(S.occup))
    N0_STAR_02 = rbind(N0_STAR_02,n0_star_02)
    BETA = rbind(BETA,beta)
    LAMBDA = LAMBDA + lambda
    OMEGA = OMEGA + omega
    if(sum(mc==M_to_use_mc)==1){
      #sample missing data for household
      Data_house_new = Data_house
      lambda_g = t(lambda[,G])
      for(kkk in 2:q){
        pr_X_miss_q = lambda_g[,d_k[kkk]:cumsum(d_k_house)[kkk]]
        Ran_unif_miss_q = runif(nrow(pr_X_miss_q))
        cumul_miss_q = pr_X_miss_q%*%upper.tri(diag(ncol(pr_X_miss_q)),diag=TRUE)
        Data_house_new[,kkk] = rowSums(Ran_unif_miss_q>cumul_miss_q) + 1L
      }
      Data_house_new_fuse = matrix(rep(c(data.matrix(Data_house_new)), rep(cn_i,q)),ncol=q)

      #sample missing data for indiv without struc zeros
      Data_indiv_new = Data_indiv
      phi_m_g = matrix(0,nrow=N,ncol=dim(phi)[1])
      for(jjj in 1:N){
        phi_m_g[jjj,] = phi[,(rep_G[jjj]+((M[jjj]-1)*FF))]
      }
      for(kkkk in 1:p){
        pr_X_miss_p = phi_m_g[-which(n_i_index==2),d_k_gm[kkkk]:cumsum(d_k_indiv)[kkkk]]
        Ran_unif_miss_p = runif(nrow(pr_X_miss_p))
        cumul_miss_p = pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
        Data_indiv_new[-which(n_i_index==2),kkkk] = rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L
      }

      #sample missing data for indiv with struc zeros HH2
      Data_house_miss_02 = Data_house_new[unique(house_index[which(n_i_index==2)]),]
      FFF_house_miss_02 = matrix(rep(cumsum(c(0,d_k_house[,-q])),
                                          each=nrow(Data_house_miss_02)),ncol=q)
      lambda_index_house_miss_02 = data.matrix(Data_house_miss_02)+FFF_house_miss_02
      M_X_miss_02 = M[which(n_i_index==2)]
      G_X_miss_02 = G[unique(house_index[which(n_i_index==2)])]
      struc_zero_pass_02 = All.poss.combs_02[which(struc_zero_chk_counter_02==1),]
      Data_struc_zero_pass_02 = t(matrix(c(t(struc_zero_pass_02)),nrow=p))
      FFF_struc_zero_pass_02 = matrix(rep(cumsum(c(0,d_k_indiv[,-p])),
                             each=nrow(Data_struc_zero_pass_02)),ncol=p)
      phi_index_struc_zero_pass_02 = data.matrix(Data_struc_zero_pass_02)+
        FFF_struc_zero_pass_02
      prob_X_miss_02_p = sampleXMiss02Func(lambda_index_house_miss_02,
              phi_index_struc_zero_pass_02,lambda,pii,phi,omega,G_X_miss_02,M_X_miss_02,FF)
      Ran_unif_X_miss_02_p = runif(nrow(prob_X_miss_02_p))
      cumul_X_miss_02_p = prob_X_miss_02_p%*%upper.tri(diag(ncol(prob_X_miss_02_p)),diag=TRUE)
      X_miss_02_p_poss_index = rowSums(Ran_unif_X_miss_02_p>cumul_X_miss_02_p) + 1L
      Data_indiv_new[which(n_i_index==2),] =
        t(matrix(c(t(struc_zero_pass_02[X_miss_02_p_poss_index,])),nrow=p))

      Data_house_new_fuse = cbind(Data_house_new_fuse,Data_indiv_new)
      Data_house_new_fuse = apply(Data_house_new_fuse,2,as.factor)
      Data_house_new_fuse = data.frame(Data_house_new_fuse)

      dp_imput_indiv = rbind(dp_imput_indiv,Data_indiv_new)
      dp_imput_house = rbind(dp_imput_house,Data_house_new)
      dp_imput = rbind(dp_imput,Data_house_new_fuse)
    }
  }
}
(proc.time() - ptm)/60

mean(ALPHA)
mean(BETA)
round(sort(colMeans(PII),decreasing = T),4)
sort(pi_true[,1],decreasing = T)
LAMBDA/(n_iter-burn_in)
lambda_true
OMEGA/(n_iter-burn_in)
omega_true

ALPHA.mcmc = mcmc(ALPHA)
plot(ALPHA.mcmc)
effectiveSize(ALPHA.mcmc)
geweke.diag(ALPHA.mcmc)
autocorr.plot(ALPHA.mcmc)

BETA.mcmc = mcmc(BETA)
plot(BETA.mcmc)
effectiveSize(BETA.mcmc)
geweke.diag(BETA.mcmc)
autocorr.plot(BETA.mcmc)

N0_STAR_02.mcmc = mcmc(N0_STAR_02)
plot(N0_STAR_02.mcmc)
effectiveSize(N0_STAR_02.mcmc)
geweke.diag(N0_STAR_02.mcmc)
autocorr.plot(N0_STAR_02.mcmc)

write.csv(dp_imput_house, file = "dp_imput_house.csv",row.names = FALSE)
write.csv(dp_imput_indiv, file = "dp_imput_indiv.csv",row.names = FALSE)
write.csv(dp_imput, file = "dp_imput.csv",row.names = FALSE)

#Combined Data Probabilities
Data_house_fuse = matrix(rep(c(data.matrix(Data_house)), rep(cn_i,q)),ncol=q)
Data_house_fuse = cbind(Data_house_fuse,Data_indiv)
Data_house_fuse = apply(Data_house_fuse,2,as.factor)
Data_house_fuse = data.frame(Data_house_fuse)
Data.nomiss = Data_house_fuse
mm=MM
all.biprcomb = combn(p+q,2)
all.triprcomb = combn(p+q,3)
dp.imput = read.csv("dp_imput.csv",header=TRUE)
# Compute marginal probs
dp.qbarMmarg = dp.bMmarg = dp.ubarMmarg = NULL
margprnomiss = margvnomiss = NULL
for(j in 1:(p+q)){
  dp.margpr.j = dp.margv.j = matrix(0,nrow=mm,ncol=length(unique(Data.nomiss[,j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((n*(k-1))+1):(n*k),]
    dp.margpr.j[k,] = as.data.frame(table(k.imp3[,j]))$Freq/dim(k.imp3)[1]
    dp.margv.j[k,] = (dp.margpr.j[k,]*(1-dp.margpr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMmarg = rbind(dp.qbarMmarg,matrix(apply(dp.margpr.j,2,mean),ncol=1))
  dp.bMmarg = rbind(dp.bMmarg,matrix(apply(dp.margpr.j,2,var),ncol=1))
  dp.ubarMmarg = rbind(dp.ubarMmarg,matrix(apply(dp.margv.j,2,mean),ncol=1))

  margprnomiss.j = as.data.frame(table(Data.nomiss[,j]))$Freq/dim(Data.nomiss)[1]
  margprnomiss.j = matrix(margprnomiss.j,ncol=1)
  margvnomiss.j = (margprnomiss.j*(1-margprnomiss.j))/dim(Data.nomiss)[1]
  margvnomiss.j = matrix(margvnomiss.j,ncol=1)
  margprnomiss = rbind(margprnomiss,margprnomiss.j)
  margvnomiss = rbind(margvnomiss,margvnomiss.j)
}

# Compute Bivariate Probs
dp.qbarMbi = dp.bMbi = dp.ubarMbi = NULL
biprnomiss = bivnomiss = NULL
for(j in 1:dim(all.biprcomb)[2]){
  comb.j = all.biprcomb[,j]
  dp.bipr.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,comb.j])))
  dp.biv.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,comb.j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((n*(k-1))+1):(n*k),]
    dp.bipr.j[k,] = as.data.frame(table(k.imp3[,comb.j]))$Freq/dim(k.imp3)[1]
    dp.biv.j[k,] = (dp.bipr.j[k,]*(1-dp.bipr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMbi = rbind(dp.qbarMbi,matrix(apply(dp.bipr.j,2,mean),ncol=1))
  dp.bMbi = rbind(dp.bMbi,matrix(apply(dp.bipr.j,2,var),ncol=1))
  dp.ubarMbi = rbind(dp.ubarMbi,matrix(apply(dp.biv.j,2,mean),ncol=1))

  biprnomiss.j = as.data.frame(table(Data.nomiss[,comb.j]))$Freq/dim(Data.nomiss)[1]
  biprnomiss.j = matrix(biprnomiss.j,ncol=1)
  bivnomiss.j = (biprnomiss.j*(1-biprnomiss.j))/dim(Data.nomiss)[1]
  bivnomiss.j = matrix(bivnomiss.j,ncol=1)
  biprnomiss = rbind(biprnomiss,biprnomiss.j)
  bivnomiss = rbind(bivnomiss,bivnomiss.j)
}

# Compute Trivariate Probs
dp.qbarMtri = dp.bMtri = dp.ubarMtri = NULL
triprnomiss = trivnomiss = NULL
for(j in 1:dim(all.triprcomb)[2]){
  combtri.j = all.triprcomb[,j]
  dp.tripr.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,combtri.j])))
  dp.triv.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,combtri.j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((n*(k-1))+1):(n*k),]
    dp.tripr.j[k,] = as.data.frame(table(k.imp3[,combtri.j]))$Freq/dim(k.imp3)[1]
    dp.triv.j[k,] = (dp.tripr.j[k,]*(1-dp.tripr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMtri = rbind(dp.qbarMtri,matrix(apply(dp.tripr.j,2,mean),ncol=1))
  dp.bMtri = rbind(dp.bMtri,matrix(apply(dp.tripr.j,2,var),ncol=1))
  dp.ubarMtri = rbind(dp.ubarMtri,matrix(apply(dp.triv.j,2,mean),ncol=1))

  triprnomiss.j = as.data.frame(table(Data.nomiss[,combtri.j]))$Freq/dim(Data.nomiss)[1]
  triprnomiss.j = matrix(triprnomiss.j,ncol=1)
  trivnomiss.j = (triprnomiss.j*(1-triprnomiss.j))/dim(Data.nomiss)[1]
  trivnomiss.j = matrix(trivnomiss.j,ncol=1)
  triprnomiss = rbind(triprnomiss,triprnomiss.j)
  trivnomiss = rbind(trivnomiss,trivnomiss.j)
}

margprnomisssim = margprnomiss
margvnomisssim = margvnomiss
biprnomisssim = biprnomiss
bivnomisssim = bivnomiss
triprnomisssim = triprnomiss
trivnomisssim = trivnomiss

dp.qbarmargsim = dp.qbarMmarg
dp.bmargsim = dp.bMmarg
dp.ubarmargsim = dp.ubarMmarg
dp.qbarbisim = dp.qbarMbi
dp.bbisim = dp.bMbi
dp.ubarbisim = dp.ubarMbi
dp.qbartrisim = dp.qbarMtri
dp.btrisim = dp.bMtri
dp.ubartrisim = dp.ubarMtri

## Check Convergence to truth
comparemarg = round(cbind(dp.qbarmargsim,margprnomisssim),4)
colnames(comparemarg) = c("DP","Truth")
write.csv(comparemarg,"comparemarg.csv", row.names = FALSE)

comparebi = round(cbind(dp.qbarbisim,biprnomisssim),4)
colnames(comparebi) = c("DP","Truth")
write.csv(comparebi,"comparebi.csv", row.names = FALSE)

comparetri = round(cbind(dp.qbartrisim,triprnomisssim),4)
colnames(comparetri) = c("DP","Truth")
write.csv(comparetri,"comparetri.csv", row.names = FALSE)

comparemarg = read.csv("comparemarg.csv",header=TRUE)
comparebi = read.csv("comparebi.csv",header=TRUE)
comparetri = read.csv("comparetri.csv",header=TRUE)

par(mfrow=c(1,3),oma=c(0,0,2,0))
plot(comparemarg[,2],comparemarg[,1],col="red",lwd=1,
     xlab="Truth",ylab="DP",main="Marginal Probabilities"); grid(nx=20)
plot(comparebi[,2],comparebi[,1],col="dark blue",lwd=1,
     xlab="Truth",ylab="DP",main="Bivariate Probabilities"); grid(nx=20)
plot(comparetri[,2],comparetri[,1],col="black",bg = "gold",lwd=1,
     xlab="Truth",ylab="DP",main="Trivariate Probabilities"); grid(nx=20)
title("Combined Data Probabilities", outer=TRUE)


#Household Probabilities
Data.nomiss = Data_house
mm=MM
all.biprcomb = combn(q,2)
all.triprcomb = combn(q,3)
dp.imput = read.csv("dp_imput_house.csv",header=TRUE)
# Compute marginal probs
dp.qbarMmarg = dp.bMmarg = dp.ubarMmarg = NULL
margprnomiss = margvnomiss = NULL
for(j in 1:q){
  dp.margpr.j = dp.margv.j = matrix(0,nrow=mm,ncol=length(unique(Data.nomiss[,j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((n*(k-1))+1):(n*k),]
    dp.margpr.j[k,] = as.data.frame(table(k.imp3[,j]))$Freq/dim(k.imp3)[1]
    dp.margv.j[k,] = (dp.margpr.j[k,]*(1-dp.margpr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMmarg = rbind(dp.qbarMmarg,matrix(apply(dp.margpr.j,2,mean),ncol=1))
  dp.bMmarg = rbind(dp.bMmarg,matrix(apply(dp.margpr.j,2,var),ncol=1))
  dp.ubarMmarg = rbind(dp.ubarMmarg,matrix(apply(dp.margv.j,2,mean),ncol=1))

  margprnomiss.j = as.data.frame(table(Data.nomiss[,j]))$Freq/dim(Data.nomiss)[1]
  margprnomiss.j = matrix(margprnomiss.j,ncol=1)
  margvnomiss.j = (margprnomiss.j*(1-margprnomiss.j))/dim(Data.nomiss)[1]
  margvnomiss.j = matrix(margvnomiss.j,ncol=1)
  margprnomiss = rbind(margprnomiss,margprnomiss.j)
  margvnomiss = rbind(margvnomiss,margvnomiss.j)
}

# Compute Bivariate Probs
dp.qbarMbi = dp.bMbi = dp.ubarMbi = NULL
biprnomiss = bivnomiss = NULL
for(j in 1:dim(all.biprcomb)[2]){
  comb.j = all.biprcomb[,j]
  dp.bipr.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,comb.j])))
  dp.biv.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,comb.j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((n*(k-1))+1):(n*k),]
    dp.bipr.j[k,] = as.data.frame(table(k.imp3[,comb.j]))$Freq/dim(k.imp3)[1]
    dp.biv.j[k,] = (dp.bipr.j[k,]*(1-dp.bipr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMbi = rbind(dp.qbarMbi,matrix(apply(dp.bipr.j,2,mean),ncol=1))
  dp.bMbi = rbind(dp.bMbi,matrix(apply(dp.bipr.j,2,var),ncol=1))
  dp.ubarMbi = rbind(dp.ubarMbi,matrix(apply(dp.biv.j,2,mean),ncol=1))

  biprnomiss.j = as.data.frame(table(Data.nomiss[,comb.j]))$Freq/dim(Data.nomiss)[1]
  biprnomiss.j = matrix(biprnomiss.j,ncol=1)
  bivnomiss.j = (biprnomiss.j*(1-biprnomiss.j))/dim(Data.nomiss)[1]
  bivnomiss.j = matrix(bivnomiss.j,ncol=1)
  biprnomiss = rbind(biprnomiss,biprnomiss.j)
  bivnomiss = rbind(bivnomiss,bivnomiss.j)
}

# Compute Trivariate Probs
dp.qbarMtri = dp.bMtri = dp.ubarMtri = NULL
triprnomiss = trivnomiss = NULL
for(j in 1:dim(all.triprcomb)[2]){
  combtri.j = all.triprcomb[,j]
  dp.tripr.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,combtri.j])))
  dp.triv.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,combtri.j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((n*(k-1))+1):(n*k),]
    dp.tripr.j[k,] = as.data.frame(table(k.imp3[,combtri.j]))$Freq/dim(k.imp3)[1]
    dp.triv.j[k,] = (dp.tripr.j[k,]*(1-dp.tripr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMtri = rbind(dp.qbarMtri,matrix(apply(dp.tripr.j,2,mean),ncol=1))
  dp.bMtri = rbind(dp.bMtri,matrix(apply(dp.tripr.j,2,var),ncol=1))
  dp.ubarMtri = rbind(dp.ubarMtri,matrix(apply(dp.triv.j,2,mean),ncol=1))

  triprnomiss.j = as.data.frame(table(Data.nomiss[,combtri.j]))$Freq/dim(Data.nomiss)[1]
  triprnomiss.j = matrix(triprnomiss.j,ncol=1)
  trivnomiss.j = (triprnomiss.j*(1-triprnomiss.j))/dim(Data.nomiss)[1]
  trivnomiss.j = matrix(trivnomiss.j,ncol=1)
  triprnomiss = rbind(triprnomiss,triprnomiss.j)
  trivnomiss = rbind(trivnomiss,trivnomiss.j)
}

margprnomisssim = margprnomiss
margvnomisssim = margvnomiss
biprnomisssim = biprnomiss
bivnomisssim = bivnomiss
triprnomisssim = triprnomiss
trivnomisssim = trivnomiss

dp.qbarmargsim = dp.qbarMmarg
dp.bmargsim = dp.bMmarg
dp.ubarmargsim = dp.ubarMmarg
dp.qbarbisim = dp.qbarMbi
dp.bbisim = dp.bMbi
dp.ubarbisim = dp.ubarMbi
dp.qbartrisim = dp.qbarMtri
dp.btrisim = dp.bMtri
dp.ubartrisim = dp.ubarMtri

## Check Convergence to truth
comparemarghouse = round(cbind(dp.qbarmargsim,margprnomisssim),4)
colnames(comparemarghouse) = c("DP","Truth")
write.csv(comparemarghouse,"comparemarghouse.csv", row.names = FALSE)

comparebihouse = round(cbind(dp.qbarbisim,biprnomisssim),4)
colnames(comparebihouse) = c("DP","Truth")
write.csv(comparebihouse,"comparebihouse.csv", row.names = FALSE)

comparetrihouse = round(cbind(dp.qbartrisim,triprnomisssim),4)
colnames(comparetrihouse) = c("DP","Truth")
write.csv(comparetrihouse,"comparetrihouse.csv", row.names = FALSE)

comparemarghouse = read.csv("comparemarghouse.csv",header=TRUE)
comparebihouse = read.csv("comparebihouse.csv",header=TRUE)
comparetrihouse = read.csv("comparetrihouse.csv",header=TRUE)

par(mfrow=c(1,3),oma=c(0,0,2,0))
plot(comparemarghouse[,2],comparemarghouse[,1],col="red",lwd=1,
     xlab="Truth",ylab="DP",main="Marginal Probabilities"); grid(nx=20)
plot(comparebihouse[,2],comparebihouse[,1],col="dark blue",lwd=1,
     xlab="Truth",ylab="DP",main="Bivariate Probabilities"); grid(nx=20)
plot(comparetrihouse[,2],comparetrihouse[,1],col="black",bg = "gold",lwd=1,
     xlab="Truth",ylab="DP",main="Trivariate Probabilities"); grid(nx=20)
title("Household Probabilities", outer=TRUE)


#Individual Probabilities
Data.nomiss = Data_indiv
mm=MM
all.biprcomb = combn(p,2)
all.triprcomb = combn(p,3)
dp.imput = read.csv("dp_imput_indiv.csv",header=TRUE)
# Compute marginal probs
dp.qbarMmarg = dp.bMmarg = dp.ubarMmarg = NULL
margprnomiss = margvnomiss = NULL
for(j in 1:p){
  dp.margpr.j = dp.margv.j = matrix(0,nrow=mm,ncol=length(unique(Data.nomiss[,j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((N*(k-1))+1):(N*k),]
    dp.margpr.j[k,] = as.data.frame(table(k.imp3[,j]))$Freq/dim(k.imp3)[1]
    dp.margv.j[k,] = (dp.margpr.j[k,]*(1-dp.margpr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMmarg = rbind(dp.qbarMmarg,matrix(apply(dp.margpr.j,2,mean),ncol=1))
  dp.bMmarg = rbind(dp.bMmarg,matrix(apply(dp.margpr.j,2,var),ncol=1))
  dp.ubarMmarg = rbind(dp.ubarMmarg,matrix(apply(dp.margv.j,2,mean),ncol=1))

  margprnomiss.j = as.data.frame(table(Data.nomiss[,j]))$Freq/dim(Data.nomiss)[1]
  margprnomiss.j = matrix(margprnomiss.j,ncol=1)
  margvnomiss.j = (margprnomiss.j*(1-margprnomiss.j))/dim(Data.nomiss)[1]
  margvnomiss.j = matrix(margvnomiss.j,ncol=1)
  margprnomiss = rbind(margprnomiss,margprnomiss.j)
  margvnomiss = rbind(margvnomiss,margvnomiss.j)
}

# Compute Bivariate Probs
dp.qbarMbi = dp.bMbi = dp.ubarMbi = NULL
biprnomiss = bivnomiss = NULL
for(j in 1:dim(all.biprcomb)[2]){
  comb.j = all.biprcomb[,j]
  dp.bipr.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,comb.j])))
  dp.biv.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,comb.j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((N*(k-1))+1):(N*k),]
    dp.bipr.j[k,] = as.data.frame(table(k.imp3[,comb.j]))$Freq/dim(k.imp3)[1]
    dp.biv.j[k,] = (dp.bipr.j[k,]*(1-dp.bipr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMbi = rbind(dp.qbarMbi,matrix(apply(dp.bipr.j,2,mean),ncol=1))
  dp.bMbi = rbind(dp.bMbi,matrix(apply(dp.bipr.j,2,var),ncol=1))
  dp.ubarMbi = rbind(dp.ubarMbi,matrix(apply(dp.biv.j,2,mean),ncol=1))

  biprnomiss.j = as.data.frame(table(Data.nomiss[,comb.j]))$Freq/dim(Data.nomiss)[1]
  biprnomiss.j = matrix(biprnomiss.j,ncol=1)
  bivnomiss.j = (biprnomiss.j*(1-biprnomiss.j))/dim(Data.nomiss)[1]
  bivnomiss.j = matrix(bivnomiss.j,ncol=1)
  biprnomiss = rbind(biprnomiss,biprnomiss.j)
  bivnomiss = rbind(bivnomiss,bivnomiss.j)
}

# Compute Trivariate Probs
dp.qbarMtri = dp.bMtri = dp.ubarMtri = NULL
triprnomiss = trivnomiss = NULL
for(j in 1:dim(all.triprcomb)[2]){
  combtri.j = all.triprcomb[,j]
  dp.tripr.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,combtri.j])))
  dp.triv.j = matrix(0,nrow=mm,ncol=length(table(Data.nomiss[,combtri.j])))
  for(k in 1:mm){
    k.imp3 = dp.imput[((N*(k-1))+1):(N*k),]
    dp.tripr.j[k,] = as.data.frame(table(k.imp3[,combtri.j]))$Freq/dim(k.imp3)[1]
    dp.triv.j[k,] = (dp.tripr.j[k,]*(1-dp.tripr.j[k,]))/dim(k.imp3)[1]
  }
  dp.qbarMtri = rbind(dp.qbarMtri,matrix(apply(dp.tripr.j,2,mean),ncol=1))
  dp.bMtri = rbind(dp.bMtri,matrix(apply(dp.tripr.j,2,var),ncol=1))
  dp.ubarMtri = rbind(dp.ubarMtri,matrix(apply(dp.triv.j,2,mean),ncol=1))

  triprnomiss.j = as.data.frame(table(Data.nomiss[,combtri.j]))$Freq/dim(Data.nomiss)[1]
  triprnomiss.j = matrix(triprnomiss.j,ncol=1)
  trivnomiss.j = (triprnomiss.j*(1-triprnomiss.j))/dim(Data.nomiss)[1]
  trivnomiss.j = matrix(trivnomiss.j,ncol=1)
  triprnomiss = rbind(triprnomiss,triprnomiss.j)
  trivnomiss = rbind(trivnomiss,trivnomiss.j)
}

margprnomisssim = margprnomiss
margvnomisssim = margvnomiss
biprnomisssim = biprnomiss
bivnomisssim = bivnomiss
triprnomisssim = triprnomiss
trivnomisssim = trivnomiss

dp.qbarmargsim = dp.qbarMmarg
dp.bmargsim = dp.bMmarg
dp.ubarmargsim = dp.ubarMmarg
dp.qbarbisim = dp.qbarMbi
dp.bbisim = dp.bMbi
dp.ubarbisim = dp.ubarMbi
dp.qbartrisim = dp.qbarMtri
dp.btrisim = dp.bMtri
dp.ubartrisim = dp.ubarMtri

## Check Convergence to truth
comparemargindiv = round(cbind(dp.qbarmargsim,margprnomisssim),4)
colnames(comparemargindiv) = c("DP","NoMiss")
write.csv(comparemargindiv,"comparemargindiv.csv", row.names = FALSE)

comparebiindiv = round(cbind(dp.qbarbisim,biprnomisssim),4)
colnames(comparebiindiv) = c("DP","NoMiss")
write.csv(comparebiindiv,"comparebiindiv.csv", row.names = FALSE)

comparetriindiv = round(cbind(dp.qbartrisim,triprnomisssim),4)
colnames(comparetriindiv) = c("DP","NoMiss")
write.csv(comparetriindiv,"comparetriindiv.csv", row.names = FALSE)

comparemargindiv = read.csv("comparemargindiv.csv",header=TRUE)
comparebiindiv = read.csv("comparebiindiv.csv",header=TRUE)
comparetriindiv = read.csv("comparetriindiv.csv",header=TRUE)

par(mfrow=c(1,3),oma=c(0,0,2,0))
plot(comparemargindiv[,2],comparemargindiv[,1],col="red",lwd=1,
     xlab="Truth",ylab="DP",main="Marginal Probabilities"); grid(nx=20)
plot(comparebiindiv[,2],comparebiindiv[,1],col="dark blue",lwd=1,
     xlab="Truth",ylab="DP",main="Bivariate Probabilities"); grid(nx=20)
plot(comparetriindiv[,2],comparetriindiv[,1],col="black",bg = "gold",lwd=1,
     xlab="Truth",ylab="DP",main="Trivariate Probabilities"); grid(nx=20)
title("Individual {Within Household} Probabilities", outer=TRUE)






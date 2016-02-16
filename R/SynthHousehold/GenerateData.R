
#################
# I added more structural zeros for households of size two
## 1/26/2016 -- I also edit the probabilities for the structural zeros

##Sample Data
library(DirichletReg)
library(SynthHousehold)

n = 10000
F_true = 4
S_true = 3
pi_true = c(0.50,0.30,0.15,0.05)
omega_true = matrix(c(0.55,0.3,0.15,0.3,0.3,0.4,0.2,
                      0.5,0.3,0.25,0.6,0.15),ncol=3,byrow=T)
q = 3
p = 3

#First households
G = sample(c(1:F_true),n,replace=TRUE,prob=pi_true)
d_k_house = c(2,2,2) #First one is household size, 2 categories each
d_k = 1+cumsum(c(0,d_k_house[-q]))
lambda_true = matrix(0,nrow=sum(d_k_house),ncol=F_true)
for(g in 1:F_true){
  lambda_true[d_k[1]:cumsum(d_k_house)[1],g] = t(rdirichlet(1,c(1,50)))
  lambda_true[d_k[2]:cumsum(d_k_house)[2],g] = t(rdirichlet(1,sample(c(1,3),2,replace=TRUE)))
  lambda_true[d_k[3]:cumsum(d_k_house)[3],g] = t(rdirichlet(1,sample(c(70,5),2,replace=TRUE)))
}
X_house = matrix(0,ncol=q,nrow=n)
for(i in 1:n){
  for(k in 1:q){
    X_house[i,k] = sample(c(1:d_k_house[k]),1,replace=TRUE,
                          prob=lambda_true[d_k[k]:cumsum(d_k_house)[k],G[i]])
  }
}
#X_house[,1] = 2  #only households of size 2

#Now individuals
d_k_indiv = c(2,2,2) #2 categories each, let last one be the relationship variable in a sense
d_k_gm = 1+cumsum(c(0,d_k_indiv[-p]))
n_i = X_house[,1]
house_index = rep(c(1:n),n_i)
n_i_index = rep(n_i,n_i)
rep_G = rep(G,n_i)
N = length(house_index)
M = matrix(0,nrow=N,ncol=1)
for(j in 1:N){
  M[j] = sample(c(1:S_true),1,replace=TRUE,prob=omega_true[rep_G[j],])
}
phi_true = array(0,dim=c(sum(d_k_indiv),S_true,F_true))
for(m in 1:S_true){
  for(g in 1:F_true){
    phi_true[d_k_gm[1]:cumsum(d_k_indiv)[1],m,g] = t(rdirichlet(1,sample(c(10,100),2,replace=TRUE)))
    phi_true[d_k_gm[2]:cumsum(d_k_indiv)[2],m,g] = t(rdirichlet(1,sample(c(1,5),2,replace=TRUE)))
    phi_true[d_k_gm[3]:cumsum(d_k_indiv)[3],m,g] = t(rdirichlet(1,sample(c(200,5),2,replace=TRUE)))
  } }
X_indiv = matrix(0,ncol=p,nrow=N)
for(j in 1:N){
  for(k in 1:p){
    X_indiv[j,k] = sample(c(1:d_k_indiv[k]),1,replace=TRUE,
                          prob=phi_true[d_k_gm[k]:cumsum(d_k_indiv)[k],M[j],rep_G[j]])
  }
}

##Structural Zeros for household 1
#X_indiv_check_counter = X_house_check_counter = NULL
#for(str_ch in unique(house_index[which(n_i_index==1)])){
#  hh_to_check = X_indiv[which(house_index==str_ch),]
#  if(hh_to_check[1,1] == 1 & hh_to_check[1,2] == 2 |
#     hh_to_check[1,2] == 1 & hh_to_check[1,3] == 2 |
#     hh_to_check[1,1] == 2 & hh_to_check[1,3] == 1){
#    X_indiv_check_counter = c(X_indiv_check_counter,which(house_index==str_ch))
#    X_house_check_counter = c(X_house_check_counter,str_ch)
#  }
#}
#X_indiv = X_indiv[-X_indiv_check_counter,]
#X_house = X_house[-X_house_check_counter,]
#G = G[-X_house_check_counter]
#M = M[-X_indiv_check_counter]

#Structural Zeros for household 2
X_indiv_check_counter = X_house_check_counter = NULL
for(str_ch in unique(house_index[which(n_i_index==2)])){
  hh_to_check = X_indiv[which(house_index==str_ch),]
  hh_to_check = hh_to_check[order(hh_to_check[,3]),]
  if(hh_to_check[1,1] == 1 & hh_to_check[1,2] == 2 |
     hh_to_check[1,1] == 2 & hh_to_check[1,3] == 1 |
     hh_to_check[2,1] == 1 & hh_to_check[2,2] == 2 |
     hh_to_check[2,1] == 2 & hh_to_check[2,3] == 1 |

     hh_to_check[1,3] == 1 & hh_to_check[2,3] == 1 |
     hh_to_check[1,1] == 2 & hh_to_check[1,3] == 1 & hh_to_check[2,1] == 2 & hh_to_check[2,3] == 2 |
     hh_to_check[1,1] == 1 & hh_to_check[1,3] == 2 & hh_to_check[2,1] == 1 & hh_to_check[2,3] == 2 |
     hh_to_check[1,2] == 1 & hh_to_check[1,3] == 1 & hh_to_check[2,2] == 2 & hh_to_check[2,3] == 2 |
     hh_to_check[1,1] == 1 & hh_to_check[1,2] == 2 & hh_to_check[2,1] == 2 & hh_to_check[2,2] == 1 ){

    X_indiv_check_counter = c(X_indiv_check_counter,which(house_index==str_ch))
    X_house_check_counter = c(X_house_check_counter,str_ch)
  }
}
X_indiv = X_indiv[-X_indiv_check_counter,]
X_house = X_house[-X_house_check_counter,]
G = G[-X_house_check_counter]
M = M[-X_indiv_check_counter]

write.csv(X_house, file = "X_house.csv",row.names = FALSE)
write.csv(X_indiv, file = "X_indiv.csv",row.names = FALSE)
write.csv(pi_true, file = "pi_true.csv",row.names = FALSE)
write.csv(omega_true, file = "omega_true.csv",row.names = FALSE)
write.csv(lambda_true, file = "lambda_true.csv",row.names = FALSE)
write.csv(phi_true, file = "phi_true.csv",row.names = FALSE)
write.csv(G, file = "G.csv",row.names = FALSE)
write.csv(M, file = "M.csv",row.names = FALSE)

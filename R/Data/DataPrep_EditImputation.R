

################### Phase One: One Time Data Preparation ##################
rm(list = ls())
set.seed(1234)
Rcpp::sourceCpp('checkSZ.cpp')
###### 1: Import Data
House <- read.csv("House.csv",header=T)
Indiv <- read.csv("Indiv.csv",header=T)


###### 2: Remove Households with size < 2 and > 6
House <- House[which(House$NP >= 2 & House$NP <= 6),]


###### 3: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
House$TEN[which(House$TEN == 2)] <- 1
House$TEN[which(House$TEN == 3)] <- 2


###### 4: Take a sample of size 3000 Households
sample_size <- 3000
samp_index <- sort(sample(1:nrow(House),sample_size,replace=F))
House <- House[samp_index,]


###### 5: Pick the same households in the indiv data
pick_index <- is.element(Indiv$SERIALNO,House$SERIALNO)
Indiv <- Indiv[pick_index,]


###### 6: Recode within-household variables
###### 6a: First, the relationship variable
Indiv$RELP[which(Indiv$RELP == 11 | Indiv$RELP == 12 | Indiv$RELP == 13)] <- 12 #Boarder, roommate or partner
Indiv$RELP[which(Indiv$RELP == 14 | Indiv$RELP == 15)] <- 13 #Other non-relative or foster child
Indiv$RELP[which(Indiv$RELP == 10)] <- 11 #Other relative
Indiv$RELP[which(Indiv$RELP == 9)] <- 10 #Child-in-law
Indiv$RELP[which(Indiv$RELP == 8)] <- 9 #Parent-in-law
Indiv$RELP[which(Indiv$RELP == 7)] <- 8 #Grandchild
Indiv$RELP[which(Indiv$RELP == 6)] <- 7 #Parent
Indiv$RELP[which(Indiv$RELP == 5)] <- 6 #Sibling
Indiv$RELP[which(Indiv$RELP == 4)] <- 5 #Stepchild
Indiv$RELP[which(Indiv$RELP == 3)] <- 4 #Adopted child
Indiv$RELP[which(Indiv$RELP == 2)] <- 3 #Biological child
Indiv$RELP[which(Indiv$RELP == 1)] <- 2 #Spouse
Indiv$RELP[which(Indiv$RELP == 0)] <- 1 #Household head
###### 6b: Next, the race variable
Indiv$RAC3P[which(Indiv$RAC3P == 4 | Indiv$RAC3P == 8| Indiv$RAC3P == 9 | Indiv$RAC3P == 10)] <- 6
Indiv$RAC3P[which(Indiv$RAC3P == 5)] <- 4
Indiv$RAC3P[which(Indiv$RAC3P == 7)] <- 5
Indiv$RAC3P[which(Indiv$RAC3P >= 11 & Indiv$RAC3P <= 15)] <- 7
Indiv$RAC3P[which(Indiv$RAC3P >= 16 & Indiv$RAC3P <= 59)] <- 8
Indiv$RAC3P[which(Indiv$RAC3P >= 60 & Indiv$RAC3P <= 100)] <- 9
###### 6c: Next, the hisp variable
Indiv$HISP[which(Indiv$HISP >= 5 & Indiv$HISP <= 24)] <- 5
###### 6d: Lastly, age
Indiv$AGEP <- Indiv$AGEP + 1L


###### 7: Make household head into household level data
HHhead_data <- Indiv[which(Indiv$SPORDER==1),]
Indiv_minHH <- Indiv[-which(Indiv$SPORDER==1),]


###### 8: Combine Household and within-household data using the following ordering:
origdata <- data.frame(HHIndex = rep(c(1:sample_size),(House$NP-1L)),
                       WithinHHIndex = Indiv_minHH$SPORDER,
                       Gender = Indiv_minHH$SEX,Race = Indiv_minHH$RAC3P,Hisp = Indiv_minHH$HISP,
                       Age = Indiv_minHH$AGEP,Relate = Indiv_minHH$RELP,
                       Owner = rep(House$TEN,(House$NP-1L)),
                       HHGender = rep(HHhead_data$SEX,(House$NP-1L)),
                       HHRace = rep(HHhead_data$RAC3P,(House$NP-1L)),
                       HHHisp = rep(HHhead_data$HISP,(House$NP-1L)),
                       HHAge = rep(HHhead_data$AGEP,(House$NP-1L)),
                       HHRelate = rep(HHhead_data$RELP,(House$NP-1L)))


###### 9: Save and load back
#write.table(origdata,"Data/origdata.txt",row.names = FALSE)
#origdata <- read.table("Data/origdata.txt",header=T)


###### 10: Separate household and individual data
n_all <- length(unique(origdata$HHIndex))
X_indiv <- NULL
X_house <- NULL
for(i in 1:n_all){
  which_indiv <- which(origdata$HHIndex==i)
  X_indiv <- rbind(X_indiv,origdata[which_indiv,c("Gender","Race","Hisp","Age","Relate")])
  X_house <- rbind(X_house,cbind(length(which_indiv),origdata[
    which_indiv[length(which_indiv)],c("Owner","HHGender","HHRace","HHHisp","HHAge","HHRelate")]))
}
colnames(X_house) <- c("HHSize","Owner","HHGender","HHRace","HHHisp","HHAge","HHRelate")
X_house <- as.data.frame(X_house)


###### 11: Now create erroneous data using the measurement error model
N <- nrow(X_indiv)
n <- nrow(X_house)
n_i <- as.numeric(as.character(X_house[,1]))
house_index <- rep(c(1:n),n_i)
p <- ncol(X_indiv)
q <- ncol(X_house)
level_indiv <- c(1:p)
level_indiv <- lapply(level_indiv, function(x) c(min(X_indiv[,x]):max(X_indiv[,x])) )
level_house <- c(1:q)
level_house <- lapply(level_house, function(x) c(min(X_house[,x]):max(X_house[,x])) )
Y_house <- X_house; Y_indiv <- X_indiv
struc_zero_variables_house <- which(is.element(colnames(X_house),c("HHGender","HHAge"))) ##gender is still included because I am still using 2012 data
struc_zero_variables_indiv <- which(is.element(colnames(X_indiv),c("Gender","Age","Relate"))) ##gender is still included because I am still using 2012 data
epsilon_indiv <- c(0.70,0.85,0.90)
epsilon_house <- c(0.65,0.80)
#epsilon_indiv <- rbeta(length(struc_zero_variables_indiv),30,5)
#epsilon_house <- rbeta(length(struc_zero_variables_house),30,5)
#subst_prob_house <- lapply(struc_zero_variables_house,function(x) rgamma(length(level_house[[x]]),50,1))
#subst_prob_house <- lapply(subst_prob_house,function(x) x/sum(x))
#subst_prob_indiv <- lapply(struc_zero_variables_indiv,function(x) rgamma(length(level_indiv[[x]]),50,1))
#subst_prob_indiv <- lapply(subst_prob_indiv,function(x) x/sum(x))
gamma <- 0.20
z_i <- rbinom(n,1,gamma)
Error_index_house <- which(z_i == 1)
E_house <- matrix(0,ncol=q,nrow=n)
E_indiv <- matrix(0,ncol=p,nrow=N)
for(i in Error_index_house){
  Error_index_indiv_i <- which(is.element(house_index,i)==TRUE)
  check_counter <- 1
  while(check_counter == 1){
    for(kq in struc_zero_variables_house){
      E_house[i,kq] <- rbinom(1,1,epsilon_house[which(struc_zero_variables_house==kq)])
      if(E_house[i,kq]==1){
        if(length(level_house[[kq]])<=2){
          Y_house[i,kq] <- (level_house[[kq]])[-(1+X_house[i,kq]-min(level_house[[kq]]))]
        } else {
          Y_house[i,kq] <- 
            #sample((level_house[[kq]])[-(1+X_house[i,kq]-min(level_house[[kq]]))],1,
            #prob=rgamma((length(level_house[[kq]])-1),50,1))
            sample((level_house[[kq]])[-(1+X_house[i,kq]-min(level_house[[kq]]))],1)
        }
      }
    }
    for(kp in struc_zero_variables_indiv){
      E_indiv[Error_index_indiv_i,kp] <- rbinom(length(Error_index_indiv_i),1,
                                                epsilon_indiv[which(struc_zero_variables_indiv==kp)])
      for(ii in 1:length(Error_index_indiv_i)){
        if(E_indiv[Error_index_indiv_i[ii],kp]==1){
          if(length(level_indiv[[kp]])<=2){
            Y_indiv[Error_index_indiv_i[ii],kp] <- 
              (level_indiv[[kp]])[-(1+X_indiv[Error_index_indiv_i[ii],kp]-min(level_indiv[[kp]]))]
          } else {
            Y_indiv[Error_index_indiv_i[ii],kp] <- 
              #sample((level_indiv[[kp]])[-(1+X_indiv[Error_index_indiv_i[ii],kp]-min(level_indiv[[kp]]))],1,
              #prob=rgamma((length(level_indiv[[kp]])-1),50,1))
              sample((level_indiv[[kp]])[-(1+X_indiv[Error_index_indiv_i[ii],kp]-min(level_indiv[[kp]]))],1)
          }
        }
      }
    }
    comb_to_check <- matrix(t(Y_indiv[Error_index_indiv_i,]),byrow=T,nrow=1)
    comb_to_check <- cbind(as.matrix(Y_house[i,(q-p+1):q]),comb_to_check) #add the household head before check
    check_counter <- checkSZ(comb_to_check,(length(Error_index_indiv_i)+1))
  }
}
E_house <- data.matrix(X_house)- data.matrix(Y_house)
E_house[E_house!=0] <- 1
E_indiv <- data.matrix(X_indiv)- data.matrix(Y_indiv)
E_indiv[E_indiv!=0] <- 1

colSums(E_house)/length(Error_index_house)
#0.0000000 0.0000000 0.9206612 0.0000000 0.0000000 0.9371901 0.0000000 
colSums(E_indiv)/length(which(is.element(house_index,Error_index_house)==TRUE))
#0.9517358 0.0000000 0.0000000 0.9364945 0.9237934 

colSums(E_house)/nrow(E_house)
colSums(E_indiv)/nrow(E_indiv)


###### 12: Add missing data
O_house <- matrix(1,ncol=q,nrow=n)
colnames(O_house) <- colnames(Y_house)
nonstruc_zero_variables_house <- c(1:ncol(Y_house))[-c(struc_zero_variables_house,which(colnames(X_house)=="HHSize"))]
O_house[,nonstruc_zero_variables_house] <- rbinom((n*length(nonstruc_zero_variables_house)),1,0.80)
O_indiv <- matrix(1,ncol=p,nrow=N)
colnames(O_indiv) <- colnames(Y_indiv)
nonstruc_zero_variables_indiv <- c(1:ncol(Y_indiv))[-struc_zero_variables_indiv]
O_indiv[,nonstruc_zero_variables_indiv] <- rbinom((N*length(nonstruc_zero_variables_indiv)),1,0.80)
Y_house[O_house==0] <- NA; Y_house$HHRelate <- 1;
Y_indiv[O_indiv==0] <- NA



origdata <- origdata[colnames(origdata)!="HHRelate"]
colnames(origdata) <- c("Hhindex","pernum","sex","race","hisp","age","relate",
                        "ownership","headsex","headrace","headhisp","headage")
origdata_truth <- origdata
origdata[,c("sex","race","hisp","age","relate")] <- Y_indiv
origdata$ownership <- rep(Y_house$Owner,n_i)
origdata$headsex <- rep(Y_house$HHGender,n_i)
origdata$headrace <- rep(Y_house$HHRace,n_i)
origdata$headhisp <- rep(Y_house$HHHisp,n_i)
origdata$headage <- rep(Y_house$HHAge,n_i)



###### 13: Save!!!
write.table(origdata,"origdata_EI.txt",row.names = FALSE)
write.table(origdata_truth,"origdata_EI_truth.txt",row.names = FALSE)
#write.table(E_matrix,"E_matrix_truth.txt",row.names = FALSE)
#write.table(E_house, file = "E_house_truth.txt",row.names = FALSE)
#write.table(E_indiv, file = "E_indiv_truth.txt",row.names = FALSE)
############################ End of Phase One #############################










###########################################################################
###########################################################################
####### MI using NDPMPM: Use rejection sampler to sample missing data #####
####### Add weighting option for sampling impossibles #####################
####### Add the Rejection sampler/sampling proposals hybrid option ########
####### Also make household head into a household level variable ##########
###########################################################################
###########################################################################

###################################### START ######################################
########################## Step 1: One Time Data Preparation ##########################
rm(list = ls())
set.seed(1234)
Rcpp::sourceCpp('checkSZ.cpp')
###### 1: Import Data
House <- read.csv("House.csv",header=T)
Indiv <- read.csv("Indiv.csv",header=T)


###### 2: Remove Households with size < 2 and > 6
House <- House[which(House$NP >= 2 & House$NP <= 6),]


###### 3: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
House$TEN[which(House$TEN == 2)] <- 1
House$TEN[which(House$TEN == 3)] <- 2


###### 4: Take a sample of size 3000 Households
sample_size <- 3000
samp_index <- sort(sample(1:nrow(House),sample_size,replace=F))
House <- House[samp_index,]


###### 5: Pick the same households in the indiv data
pick_index <- is.element(Indiv$SERIALNO,House$SERIALNO)
Indiv <- Indiv[pick_index,]


###### 6: Recode within-household variables
###### 6a: First, the relationship variable
Indiv$RELP[which(Indiv$RELP == 11 | Indiv$RELP == 12 | Indiv$RELP == 13)] <- 12 #Boarder, roommate or partner
Indiv$RELP[which(Indiv$RELP == 14 | Indiv$RELP == 15)] <- 13 #Other non-relative or foster child
Indiv$RELP[which(Indiv$RELP == 10)] <- 11 #Other relative
Indiv$RELP[which(Indiv$RELP == 9)] <- 10 #Child-in-law
Indiv$RELP[which(Indiv$RELP == 8)] <- 9 #Parent-in-law
Indiv$RELP[which(Indiv$RELP == 7)] <- 8 #Grandchild
Indiv$RELP[which(Indiv$RELP == 6)] <- 7 #Parent
Indiv$RELP[which(Indiv$RELP == 5)] <- 6 #Sibling
Indiv$RELP[which(Indiv$RELP == 4)] <- 5 #Stepchild
Indiv$RELP[which(Indiv$RELP == 3)] <- 4 #Adopted child
Indiv$RELP[which(Indiv$RELP == 2)] <- 3 #Biological child
Indiv$RELP[which(Indiv$RELP == 1)] <- 2 #Spouse
Indiv$RELP[which(Indiv$RELP == 0)] <- 1 #Household head
###### 6b: Next, the race variable
Indiv$RAC3P[which(Indiv$RAC3P == 4 | Indiv$RAC3P == 8| Indiv$RAC3P == 9 | Indiv$RAC3P == 10)] <- 6
Indiv$RAC3P[which(Indiv$RAC3P == 5)] <- 4
Indiv$RAC3P[which(Indiv$RAC3P == 7)] <- 5
Indiv$RAC3P[which(Indiv$RAC3P >= 11 & Indiv$RAC3P <= 15)] <- 7
Indiv$RAC3P[which(Indiv$RAC3P >= 16 & Indiv$RAC3P <= 59)] <- 8
Indiv$RAC3P[which(Indiv$RAC3P >= 60 & Indiv$RAC3P <= 100)] <- 9
###### 6c: Next, the hisp variable
Indiv$HISP[which(Indiv$HISP >= 5 & Indiv$HISP <= 24)] <- 5
###### 6d: Lastly, age
Indiv$AGEP <- Indiv$AGEP + 1L


###### 7: Make household head into household level data
HHhead_data <- Indiv[which(Indiv$SPORDER==1),]
Indiv_minHH <- Indiv[-which(Indiv$SPORDER==1),]


###### 8: Combine Household and within-household data using the following ordering:
origdata <- data.frame(Hhindex = rep(c(1:sample_size),(House$NP-1L)),
                       pernum = Indiv_minHH$SPORDER,
                       sex = Indiv_minHH$SEX, race = Indiv_minHH$RAC3P, hisp = Indiv_minHH$HISP,
                       age = Indiv_minHH$AGEP, relate = Indiv_minHH$RELP,
                       ownership = rep(House$TEN,(House$NP-1L)),
                       headsex = rep(HHhead_data$SEX,(House$NP-1L)),
                       headrace = rep(HHhead_data$RAC3P,(House$NP-1L)),
                       headhisp = rep(HHhead_data$HISP,(House$NP-1L)),
                       headage = rep(HHhead_data$AGEP,(House$NP-1L)))
origdata_truth <- origdata


###### 9: Now create erroneous data
N <- nrow(origdata)
n <- length(unique(origdata$Hhindex))
house_variables <- c("ownership","headsex","headrace","headhisp","headage")
indiv_variables <- c("sex","race","hisp","age","relate")
p <- length(indiv_variables)
q <- length(house_variables)
variable_levels <- c(1:ncol(origdata))
variable_levels <- lapply(variable_levels, function(x) c(min(origdata[,x]):max(origdata[,x])) )
struc_zero_variables_indiv <- which(is.element(colnames(origdata),c("sex","age","relate"))) ##sex is still included because I am still using 2012 data
struc_zero_variables_house <- which(is.element(colnames(origdata),c("headsex","headage"))) ##sex is still included because I am still using 2012 data
epsilon_indiv <- c(0.70,0.85,0.90)
epsilon_house <- c(0.65,0.80)
#epsilon <- rbeta(length(struc_zero_variables),30,5)
gamma <- 0.20
z_i <- rbinom(n,1,gamma)
hh_error_index <- which(z_i == 1)
E_matrix <- matrix(0,ncol=ncol(origdata),nrow=nrow(origdata))
for(i in hh_error_index){
  hh_error_index_i <- which(is.element(origdata$Hhindex,i)==TRUE)
  check_counter <- 1
  while(check_counter == 1){
    for(kq in struc_zero_variables_house){
      E_matrix[hh_error_index_i,kq] <- rbinom(1,1,epsilon_house[which(struc_zero_variables_house==kq)])
      if(E_matrix[hh_error_index_i[1],kq]==1){
        if(length(variable_levels[[kq]])==2){
          origdata[hh_error_index_i,kq] <- (variable_levels[[kq]])[-(1+origdata_truth[hh_error_index_i[1],kq]-min(variable_levels[[kq]]))]
        } else {
          origdata[hh_error_index_i,kq] <-
            sample((variable_levels[[kq]])[-(1+origdata_truth[hh_error_index_i[1],kq]-min(variable_levels[[kq]]))],1)
            #sample((variable_levels[[kq]])[-(1+origdata_truth[hh_error_index_i[1],kq]-min(variable_levels[[kq]]))],1)
            #prob=rgamma((length(variable_levels[[kq]])-1),50,1))
        }
      }
    }
    for(kp in struc_zero_variables_indiv){
      E_matrix[hh_error_index_i,kp] <- rbinom(length(hh_error_index_i),1,epsilon_indiv[which(struc_zero_variables_indiv==kp)])
      for(ii in 1:length(hh_error_index_i)){
        if(E_matrix[hh_error_index_i[ii],kp]==1){
          if(length(variable_levels[[kp]])==2){
            origdata[hh_error_index_i[ii],kp] <- (variable_levels[[kp]])[-(1+origdata_truth[hh_error_index_i[ii],kp]-min(variable_levels[[kp]]))]
          } else {
            origdata[hh_error_index_i[ii],kp] <-
              sample((variable_levels[[kp]])[-(1+origdata_truth[hh_error_index_i[ii],kp]-min(variable_levels[[kp]]))],1)
              #sample((variable_levels[[kp]])[-(1+origdata_truth[hh_error_index_i[ii],kp]-min(variable_levels[[kp]]))],1,
              #prob=rgamma((length(variable_levels[[kp]])-1),50,1))
          }
        }
      }
    }
    comb_to_check <- matrix(t(origdata[hh_error_index_i,indiv_variables]),byrow=T,nrow=1)
    comb_to_check <- cbind(as.matrix(origdata[i,house_variables[-1]]),1,comb_to_check) #add the household head before check
    colnames(comb_to_check) <- NULL
    check_counter <- checkSZ(comb_to_check,(length(hh_error_index_i)+1))
  }
}
E_matrix <- data.matrix(origdata_truth)- data.matrix(origdata)
E_matrix[E_matrix!=0] <- 1


colSums(E_matrix[,indiv_variables])/length(which(is.element(origdata$Hhindex,hh_error_index)==TRUE))
#0.9517358 0.0000000 0.0000000 0.9364945 0.9237934
colSums(E_matrix[!duplicated(origdata$Hhindex),house_variables])/length(hh_error_index)
#0.0000000 0.9206612 0.0000000 0.0000000 0.9371901
colSums(E_matrix)/nrow(E_matrix)


###### 10: Add missing data
O_matrix <- matrix(1,ncol=ncol(origdata),nrow=nrow(origdata))
colnames(O_matrix) <- colnames(origdata)
nonstruc_zero_variables <- which(is.element(colnames(origdata),c("race","hisp")))
O_matrix[,nonstruc_zero_variables] <- rbinom((N*length(nonstruc_zero_variables)),1,0.80)
O_matrix[,"ownership"] <- rep(rbinom(n,1,0.80),as.data.frame(table(origdata$Hhindex))$Freq)
O_matrix[,"headrace"] <- rep(rbinom(n,1,0.80),as.data.frame(table(origdata$Hhindex))$Freq)
O_matrix[,"headhisp"] <- rep(rbinom(n,1,0.80),as.data.frame(table(origdata$Hhindex))$Freq)
origdata[O_matrix==0] <- NA
origdata[1:20,]
origdata_truth[1:20,]

###### 11: Save!!!
write.table(origdata,"origdata_EI.txt",row.names = FALSE)
write.table(origdata_truth,"origdata_EI_truth.txt",row.names = FALSE)
write.table(E_matrix,"E_matrix_truth.txt",row.names = FALSE)
########################## End of Step 1 ##########################


###########################################################################
###########################################################################
#############################***Blank Space***#############################
###########################################################################
###########################################################################
###########################################################################




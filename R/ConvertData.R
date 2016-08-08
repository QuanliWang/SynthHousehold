
####### Create two functions for data_to_check

ConvertDataForward <- function(data_to_check,hh_size,p,lambda){
  n_lambdas <- length(lambda)
  ### Create the new data_to_check matrix
  new_data_to_check <- matrix(0,nrow=nrow(data_to_check), ncol=((2+n_lambdas-1)*(hh_size+1) + (hh_size+1+1)))
  ### First assign the individual level indicators, dummy for household head
  new_data_to_check[,(ncol(new_data_to_check)-hh_size+1):ncol(new_data_to_check)] <-
    data_to_check[,(ncol(data_to_check)-hh_size+1):ncol(data_to_check)]
  ### Then the household level indicator
  new_data_to_check[,(ncol(new_data_to_check)-(hh_size+1))] <- data_to_check[,(ncol(data_to_check)-hh_size)]
  ### Now assign household level indexes
  new_data_to_check[,seq.int(from=1,to=(ncol(new_data_to_check)-(hh_size+1+1)),by=(2+n_lambdas-1))] <-
    data_to_check[,1]
  ### Assign individual level indexes
  new_data_to_check[,seq.int(from=2,to=(ncol(new_data_to_check)-(hh_size+1+1)),by=(2+n_lambdas-1))] <-
    cbind(1,(data_to_check[,seq.int(from=2,to=(ncol(data_to_check)-(hh_size+1)),by=(2+p+n_lambdas-1))]+1))
  ### Next, individual level variables
  for(i in 1:p){
    ### First the household head
    new_data_to_check[,(2+i)] <- data_to_check[,(2+p+1+i)]
    ### Now everyone else
    new_data_to_check[,seq.int(from=(2+n_lambdas+1+i),to=(ncol(new_data_to_check)-(hh_size+1+1)),
                               by=(2+n_lambdas-1))] <-
      data_to_check[,seq.int(from=(2+i),to=(ncol(data_to_check)-(hh_size+1)),by=(2+p+n_lambdas-1))]
  }
  ### Lastly, household level variables
  for(j in 1:(n_lambdas-p-1)){
    ### First the household head
    new_data_to_check[,(2+p+j)] <- data_to_check[,(2+p+j)]
    ### Now everyone else
    new_data_to_check[,seq.int(from=(2+n_lambdas+p+1+j),to=(ncol(new_data_to_check)-(hh_size+1+1)),
                               by=(2+n_lambdas-1))] <- data_to_check[,(2+p+j)]
  }
  return(new_data_to_check)
}


ConvertDataBack <- function(new_data_to_check,hh_size,p,lambda){
  n_lambdas <- length(lambda)
  ### Create the new data_to_check matrix
  old_data_to_check <- matrix(0,nrow=nrow(new_data_to_check),ncol=((2+p+n_lambdas-1)*(hh_size) + (hh_size+1)))
  ### First assign the individual level indicators
  old_data_to_check[,(ncol(old_data_to_check)-hh_size+1):ncol(old_data_to_check)] <-
    new_data_to_check[,(ncol(new_data_to_check)-hh_size+1):ncol(new_data_to_check)]
  ### Then the household level indicator
  old_data_to_check[,(ncol(old_data_to_check)-hh_size)] <-
    new_data_to_check[,(ncol(new_data_to_check)-(hh_size+1))]
  ### Now assign household level indexes
  old_data_to_check[,seq.int(from=1,to=(ncol(old_data_to_check)-(hh_size+1)),by=(2+p+n_lambdas-1))] <-
    new_data_to_check[,1]
  ### Now assign individual level indexes
  old_data_to_check[,seq.int(from=2,to=(ncol(old_data_to_check)-(hh_size+1)),by=(2+p+n_lambdas-1))] <-
    new_data_to_check[,seq.int(from=(2+2+n_lambdas-1),to=(ncol(new_data_to_check)-(hh_size+1+1)),
                               by=(2+n_lambdas-1))]-1
  ### Next, individual level variables
  for(i in 1:p){
    old_data_to_check[,seq.int(from=(2+i),to=(ncol(old_data_to_check)-(hh_size+1)),by=(2+p+n_lambdas-1))] <-
      new_data_to_check[,seq.int(from=(2+n_lambdas+1+i),to=(ncol(new_data_to_check)-(hh_size+1+1)),
                                 by=(2+n_lambdas-1))]
  }
  ### Lastly, household level variables
  for(j in 1:(n_lambdas-1)){
    if(j <= (n_lambdas-p-1)){
      old_data_to_check[,seq.int(from=(2+p+j),to=(ncol(old_data_to_check)-(hh_size+1)),
                                 by=(2+p+n_lambdas-1))] <- cbind(new_data_to_check[,(2+p+j)])
    } else{
      old_data_to_check[,seq.int(from=(2+p+j),to=(ncol(old_data_to_check)-(hh_size+1)),
                                 by=(2+p+n_lambdas-1))] <- cbind(new_data_to_check[,(2+j-1)])
    }
  }
  return(old_data_to_check)
}



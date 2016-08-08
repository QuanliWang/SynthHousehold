.GenerateData <- function(hh_size,lambda, w, phi,pi, d, total.batch,possiblehhcount,howmany,synindex) {
  #save(hh_size,lambda, w, phi,pi, d, total.batch,possiblehhcount,howmany,synindex, file = "last2.RData")
  #return(NULL)
  Individuals_extra <- list()
  z_HH_extra <- list()
  HHData_extra <- list()
  batch.index <- 0
  n_possible_household <- 0
  p <- length(d)
  synIndividuals <- list()

  while (n_possible_household< possiblehhcount) {
    batch.index <- batch.index + 1
    #print(batch.index)
    #generate a batch of 10K household
    data_to_check <- samplehouseholds(phi,w, pi, d, lambda,batch.index+total.batch, howmany,hh_size)

    #impossible household
    hh_size_real <- hh_size + 1
    data_to_check_old_format <- ConvertDataForward(data_to_check,hh_size,p,lambda)

    checked.households <- checkconstraints(data_to_check,possiblehhcount-n_possible_household, hh_size+1)
    n_possible_household <- n_possible_household + checked.households$possible

    Individuals_extra[[batch.index]] <- households2individuals(checked.households$Households)
    z_HH_extra[[batch.index]] <- checked.households$Households[hh_size * 8 +1,]
    HHData_extra[[batch.index]] <- checked.households$Households[8,]
    if (synindex > 0) {
      synIndividuals[[batch.index]]  <- households2individuals(checked.households$synHouseholds)
    }
  }

  Individuals_extra <- do.call(cbind, Individuals_extra)
  z_HH_extra <- unlist(z_HH_extra)
  HHData_extra <- unlist(HHData_extra)
  if (synindex > 0) {
    synIndividuals <- do.call(cbind, synIndividuals)
  }
  batch.index <- batch.index + total.batch
  return(list(Individuals_extra = Individuals_extra,
              z_HH_extra = z_HH_extra,
              HHData_extra = HHData_extra,
              synIndividuals = synIndividuals,
              batch.index = batch.index))
}

GetImpossibleHouseholds <- function(d,ACS_count,lambda,w,phi,pi,howmany,n,synindex) {
  #save(d,ACS_count,lambda,w,phi,pi,howmany,n,synindex, file = "last.RData")
  #return(NULL)
  cumsize <- 0
  hh_size_new <-  matrix(0,nrow = length(ACS_count), ncol = 1)
  hh_index <- list()
  ImpossibleIndividuals <- list()
  z_HH_extra <- list()
  HHdata_extra <- list()
  synIndividuals_all <- list()

  ##
  total.batch <- 0
  for (hh_size in  1:(length(ACS_count))) {
    batch <- .GenerateData(hh_size,lambda, w, phi,pi, d, total.batch,ACS_count[hh_size],howmany,synindex)

    hh_size_new[hh_size] <- length(batch$z_HH_extra)
    hh_index[[hh_size]] <- cumsize + rep(1:hh_size_new[hh_size], each = hh_size)
    cumsize <- cumsize + hh_size_new[hh_size]
    ImpossibleIndividuals[[hh_size]] <- batch$Individuals_extra
    z_HH_extra[[hh_size]] <-  batch$z_HH_extra
    HHdata_extra[[hh_size]] <- rbind(batch$HHData_extra,rep(hh_size, times = hh_size_new[hh_size])) # 2 by ...
    if (synindex > 0) {
      synIndividuals_all[[hh_size]] <- batch$synIndividuals
    }
    total.batch <- batch$batch.index
  }

  ##
  hh_index <- unlist(hh_index)
  ImpossibleIndividuals <- do.call(cbind, ImpossibleIndividuals)
  z_HH_extra <- unlist(z_HH_extra)
  HHdata_extra <- do.call(cbind, HHdata_extra)
  if (synindex > 0) {
    synIndividuals_all <- do.call(cbind, synIndividuals_all)
  }

  ImpossibleIndividuals[1,] <- n + hh_index
  IndividualData_extra <- ImpossibleIndividuals[1:8,]
  z_HHdata_individual_extra <- ImpossibleIndividuals[9:10,]
  return(list(z_HHdata_individual_extra = z_HHdata_individual_extra,
              z_HH_extra = z_HH_extra,
              IndividualData_extra = IndividualData_extra,
              HHdata_extra = HHdata_extra,
              hh_size_new = hh_size_new,
              synIndividuals_all = synIndividuals_all))
}

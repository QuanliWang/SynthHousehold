GetImpossibleHouseholds <- function(d,n_star_h,lambda,omega,phi,pi,blocksize,n,synindex,HHhead_at_group_level) {
  #save(d,n_star_h,lambda,omega,phi,pi,blocksize,n,synindex, file = "last.RData")
  cumsize <- 0
  hh_size_new <-  matrix(0,nrow = length(n_star_h), ncol = 1)
  hh_index <- list()
  ImpossibleIndividuals <- list()
  G_extra <- list()
  HHdata_extra <- list()
  synIndividuals_all <- list()

  ##
  batches.done <- 0
  for (hh_size in  1:(length(n_star_h))) {
    if (HHhead_at_group_level) {
      hh_size_real <- hh_size
    } else {
      hh_size_real <- hh_size + 1
    }
    batch <- GenerateData(hh_size_real,lambda, omega, phi,pi, d, batches.done,
                          n_star_h[hh_size],blocksize,synindex,HHhead_at_group_level)

    hh_size_new[hh_size] <- length(batch$G_extra)
    hh_index[[hh_size]] <- cumsize + rep(1:hh_size_new[hh_size], each = hh_size_real)
    cumsize <- cumsize + hh_size_new[hh_size]
    ImpossibleIndividuals[[hh_size]] <- batch$Individuals_extra
    G_extra[[hh_size]] <-  batch$G_extra
    HHdata_extra[[hh_size]] <- rbind(batch$HHData_extra,rep(hh_size, times = hh_size_new[hh_size])) # length(lambda) by ...
    if (synindex > 0) {
      synIndividuals_all[[hh_size]] <- batch$synIndividuals
    }
    batches.done <- batch$batch.index
  }

  ##
  hh_index <- unlist(hh_index)
  ImpossibleIndividuals <- do.call(cbind, ImpossibleIndividuals)
  G_extra <- unlist(G_extra)
  HHdata_extra <- do.call(cbind, HHdata_extra)
  if (synindex > 0) {
    synIndividuals_all <- do.call(cbind, synIndividuals_all)
  }

  ImpossibleIndividuals[1,] <- n + hh_index
  DIM <- dim(ImpossibleIndividuals)[1] - 2
  IndividualData_extra <- ImpossibleIndividuals[1:DIM,]
  G_Individuals_and_M_extra <- ImpossibleIndividuals[(DIM+1):(DIM+2),]
  return(list(G_Individuals_and_M_extra = G_Individuals_and_M_extra,
              G_extra = G_extra,
              IndividualData_extra = IndividualData_extra,
              HHdata_extra = HHdata_extra,
              hh_size_new = hh_size_new,
              synIndividuals_all = synIndividuals_all))
}

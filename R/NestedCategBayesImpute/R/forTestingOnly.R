GenerateDataR <- function(hh_size,lambda, omega, phi,pi, d, batches_done,
                         valid_hh_needed,blocksize,synindex,HHhead_at_group_level, Parallel) {
  Individuals_extra <- list()
  G_extra <- list()
  HHData_extra <- list()
  batch.index <- 0
  valid_hh_found <- 0
  p <- length(d)
  synIndividuals <- list()

  while (valid_hh_found< valid_hh_needed) {

    parallel = 0
    if (Parallel) {parallel = 1;}

    batch.index <- batch.index + 1
    #generate a batch of 10K household
    if (HHhead_at_group_level) {
      data_to_check <- samplehouseholds(phi,omega, pi, d, lambda,
                                                             batch.index+batches_done, blocksize,hh_size, 1, parallel)
      data_to_check1 <- sampleHH(phi,omega, pi, d, lambda,
                                                       batch.index+batches_done, blocksize,hh_size,1)

      checked_households <- checkconstraints_HHhead_at_group_level(data_to_check,valid_hh_needed-valid_hh_found, hh_size)
    } else {
      data_to_check <- samplehouseholds(phi,omega, pi, d, lambda,batch.index+batches_done, blocksize,hh_size, 0, parallel)
      checked_households <- checkconstraints(data_to_check,valid_hh_needed-valid_hh_found, hh_size)
    }
    if (length(checked_households$Households) > 0) {
      DIM <- p + length(lambda) + 1
      G_extra[[batch.index]] <- checked_households$Households[hh_size * DIM +1,]
      HHData_extra[[batch.index]] <- checked_households$Households[(p+3): DIM,]
      Individuals_extra[[batch.index]] <- households2individuals(checked_households$Households, hh_size)
    }

    valid_hh_found <- valid_hh_found + checked_households$possible
    if (synindex > 0) {
      if (length(checked_households$synHouseholds) > 0) {
        synIndividuals[[batch.index]]  <- households2individuals(checked_households$synHouseholds,hh_size)
      }
    }
  }

  Individuals_extra <- do.call(cbind, Individuals_extra)
  G_extra <- unlist(G_extra)
  if (HHhead_at_group_level) {
    HHData_extra <- do.call(cbind, HHData_extra)
  } else {
    HHData_extra <- unlist(HHData_extra)
  }

  if (synindex > 0) {
    synIndividuals <- do.call(cbind, synIndividuals)
  }
  batch.index <- batch.index + batches_done
  return(list(Individuals_extra = Individuals_extra,
              G_extra = G_extra,
              HHData_extra = HHData_extra,
              synIndividuals = synIndividuals,
              batch.index = batch.index))
}

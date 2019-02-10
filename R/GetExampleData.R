PrepareData <- function(options) {
  ### Use data included in package; prepare data and specify variable indexes
  if (options$HHhead_at_group_level) {
    orig.file <- system.file("extdata","origdata.txt",package="NestedCategBayesImpute")
    orig.data <- read.table(orig.file,header = TRUE, sep = " ")
    orig.data$relate <- orig.data$relate - 1L #recode relate variable to 11 levels
    household.size <- as.data.frame(table(orig.data$Hhindex))
    household.size[,1] <- as.numeric(household.size[,1]) #conver factor to numeric
    names(household.size) <- c("Hhindex", 'householdsize')
    household <- orig.data %>% inner_join(household.size)

    individual_variable_index = c(3:7)
    household_variable_index = c(8:13) #make sure the last column represents household size

  } else {
    orig.file <- system.file("extdata","origdata_oldFormat.txt",package="NestedCategBayesImpute")
    orig.data <- read.table(orig.file,header = TRUE, sep = " ")
    orig.data$Hhindex
    household.size <- as.data.frame(table(orig.data$Hhindex))
    household.size[,1] <- as.numeric(household.size[,1])
    names(household.size) <- c("Hhindex", 'householdsize')
    household <- orig.data %>% inner_join(household.size)

    individual_variable_index = c(3:7)
    household_variable_index = c(8,9) #make sure the last column represents household size
  }
  return(list(household = household,
              individual_variable_index = individual_variable_index,
              household_variable_index = household_variable_index))
}

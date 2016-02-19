#this code prepare the example data used by both Matlab and R package.
#assume the working directory is package root.
household = read.csv("../../data/household.csv",header=TRUE)
save(household,file = "data/household.RData")

# Submit script for dissertation simulations
args <- commandArgs(trailingOnly=TRUE)

print(args)

rep <- args[1]

# use expand.grid to generate a table of simulation conditions.
# and give a seed to each replication based on a range 
# Run 100 replications in each cell
# List SCAD and boot cells 
Ns <- c(6, 18, 54 ,100) # Will need to change this
ps <- c(5, 10, 20) # This less so, to handle the ranks
rs <- c(0.2, 0.5)
wts <- c('Eq', 'Skw') # Delete this 
model <- c('Fix', 'Rand', 'Hetero') # And this

cells_table <- expand.grid(Ns, ps, rs, wts, model)
source("../runsims.R")
for (i in 1:nrow(cells_table)){

    bootvar <- (cells_table[i,1] >= 18)*( cells_table[i,2] <= 10)
    scadvar <- (cells_table[i,1] >= 18)*( cells_table[i,2] >= 10)
    print(paste("N = ", cells_table[i,1], "p = ", cells_table[i,2], "r = ", cells_table[i,3], "wt = ", cells_table[i,4], "model = ", cells_table[i,5], "seed = ", round(1e6*runif(1)), "rep = ", rep, "boots = ", bootvar, "scad = ", scadvar), sep = "")
    diss_run(N = cells_table[i,1], p = cells_table[i,2], r = cells_table[i,3], wt = cells_table[i,4], model = cells_table[i,5], seed = round(1e6*runif(1)), rep = rep, boots = bootvar, scad = scadvar)
    print("*****")
    
}

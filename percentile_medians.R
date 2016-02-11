# percentile_medians.R

#   summary: Given a matrix of results, finds the column-based rank associated with a given row, returns median and MAD of ranks aggregated across columns


### input:
#     res.M, a results matrix with results from each Monte Carlo replication stored in each column

### output:
#     medMAD.M, a two-column matrix with the same number of rows as res.M, with the first column being the median rank, and the second column the MAD


### it does:
#    1. Rank the entries in each column of res.M from largest (1) to smallest (nrows of res.M), store in ranks.M
#    2. Calculate the median of each ROW of ranks.M
#    3. Calculate the MAD of each row of ranks.M
#    4. Bind the vector of medians to the vector of MADs and return as results.M

################################################################################

percentileMedians <- function(res.M){
#    1. Rank the entries in each column of res.M from largest (1) to smallest (nrows of res.M), store in ranks.M
  ranks.M <- 1 + nrow(res.M) - apply(res.M, 2, rank)

#    2. Calculate the median of each ROW of ranks.M
  medianVec <- apply(ranks.M, 1, median)
  
#    3. Calculate the MAD of each row of ranks.M
  MADVec <- apply(ranks.M, 1, mad)
#    4. Bind the vector of medians to the vector of MADs and return as medMAD.M
  medMAD.M <- cbind(medianVec, MADVec)
  return(medMAD.M)
}

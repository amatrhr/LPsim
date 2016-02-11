# cor_resolve.R
#  summary: a function that takes a matrix of covariates and a correlation coefficient and returns a pruned version of the matrix having a correlation matrix with off-diagonal entries less than the input correlation value


### input:
##   a matrix of covariates
##   a correlation coefficient R

### output:
##   the matrix of covariates, pruned to have a correlatin matrix with off-diagonal values less than R


### it does:
##   0. Check that R >= 0.05, generate a warning that the function is likely to eliminate all variables
##   1. Calculate the absolute correlation matrix, correl.M, and set the diagonal to 0
##   2. Count the number of coefficients greater than R in each column
##   3. Remove the variable corresponding to the column with the largest count
##   4. Repeat the process recursively until no column in correl.M exceeds R
##   5. Return the pruned matrix 

################################################################################

corResolve <- function(inputMatrix, R){
##   0. Check that R >= 0.05, generate a warning that the function is likely to eliminate all variables
  if ( R <= 0.05 ){
    warning( paste("Supplied correlation coefficient,",R, ", is less than 0.05. This function is likely to eliminate all variables but one from consideration.")  )
  }

##   1. Calculate the absolute correlation matrix, correl.M, and set the diagonal to 0
  correl.M <- abs( cor(inputMatrix) - diag( 1, nrow = ncol(inputMatrix), ncol = ncol(inputMatrix) ) )
##   2. Count the number of coefficients greater than R in each column

  correlCounts <- colSums( correl.M >= R  )

  if ( max(correlCounts) == 0 ) { # termination condition
    return(inputMatrix)
  }
  else {
    ##   3. Remove the variable corresponding to the column with the largest count
    toRemove <- which.max(correlCounts)
    secondMatrix <- inputMatrix[, -toRemove]
    ##   4. Repeat the process recursively until no column in correl.M exceeds R
    corResolve(secondMatrix, R)
  }

}


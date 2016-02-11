## Script to aggregate results, which are presented in two forms:
# 1. Selection statistics--can be aggregated by cat because they are per-replication stats arranged in rows, with columns:
# 2. Coefficient statistics--stats are arranged in 13 rows per replication, with each column representing a coefficient, last representing the lambda metaparameter; this means that aggregating them using cat requires aggregation over rows with indices that are equivalent modulo 13

# Prior to running this script, comparable sets of results should be aggregated using cat
# This script can then run through each set of results and get the standard mean/sd/median/mad/IQR/fivenum stuff, and then a median polish or smth applied to get row/column effects
# 

library(plyr)
args <- commandArgs(TRUE)
# convert to the estimate and the selection file names

#args <- "polyG.simple.0.01.txt"

### Take as input the name that lists the settings of each replication, then aggregate all ver=sions of that replication; because method and model are the main determinants, do various investigations based on the results  

# Get paths right for 'production version' 

selname <- paste0("../combsel.", args)
selframe <- read.table(selname, header = FALSE, as.is = TRUE)
# VSP TPR FPR TNR
# qVSP qTPR qFPR qTNR ntVSP ntTPR ntFPR ntTNR
estname <- paste0("../combest.", args)
estframe <- read.table(estname, header = FALSE, as.is = TRUE)




sel_file_aggregate <- function(inframe, inframename){
    # input: df of results, and name of same
    # output: 4x7 or 8x7 dataframe of results, then to be written to file and fed to sel_comparison function in analyze.R 
    
output <- data.frame( mean = colMeans(inframe), sd = apply(inframe, 2, sd), median = apply(inframe, 2, median) , IQR = apply(inframe, 2, IQR), mad = apply(inframe, 2, mad), min = apply(inframe, 2, min), max = apply(inframe, 2, max))
    if(length(grep("sim", inframename)) > 0){
        rownames(output) <- c("VSP", "TPR", "FPR", "TNR")
    } else {
           rownames(output) <- c("qVSP", "qTPR", "qFPR", "qTNR", "ntVSP", "ntTPR", "ntFPR", "ntTNR")
        }
    return(output)
}



est_file_aggregate <- function(inframe, inframename){
    
    ## Input: df of results, and name of same 
    ## Output: list of data frames, corresponding to each of the stats, by coefficients.
    ## This needs to be fed into an est_compare() function, which looks at each, by condition

    
    ## 
    modframe <- inframe[, -(ncol(inframe))]
    rowidxs <- inframe[, ncol(inframe)]
    ## Loop over the different 'outcomes of the 13'
    ## STORE THE results lists as lists for each stat and return a kind of combined thing
    ## test whether it's simple or boot colnames(blasmat) <- c('mean', 'sd','min','1Q', 'Med','3Q','max' ,'qciL', 'qciU', 'ntciL', 'ntciU', "qciLen", "ntciLen")
    ## # -> this is for boot/nested
    
    if(length(grep("sim", inframename)) > 0){
        # columns are just the different coefficients, so you get the 
        output <- t(data.frame( mean = colMeans(modframe), sd = apply(modframe, 2, sd), median = apply(modframe, 2, median) , IQR = apply(modframe, 2, IQR), mad = apply(modframe, 2, mad), min = apply(modframe, 2, min), max = apply(modframe, 2, max))) # rownames would be the columns 
        
    } else {
        output <- list()
        for ( stat in c('mean', 'sd','min','X1Q', 'Med','X3Q','max' ,'qciL', 'qciU', 'ntciL', 'ntciU', "qciLen", "ntciLen")){
            output[[stat]] <- t(data.frame(mean = colMeans(modframe[rowidxs == stat,]), sd = apply(modframe[rowidxs == stat,], 2, sd), median = apply(modframe[rowidxs == stat,], 2, median) , IQR = apply(modframe[rowidxs == stat,], 2, IQR), mad = apply(modframe[rowidxs == stat,], 2, mad), min = apply(modframe[rowidxs == stat,], 2, min), max = apply(modframe[rowidxs == stat,], 2, max)))
                                        #rownames would be the column
        }
    }
return(output)
}
# How to organize results in a sensible way?

sel_results <- sel_file_aggregate(inframe = selframe, inframename = selname)

write.table(sel_results, paste0("../summ.sel.", args), quote = FALSE)


est_results <- est_file_aggregate(inframe = estframe, inframename = estname)

if(length(grep("sim", estname)) > 0){#
    est_return <- est_results
} else {
    est_return <- rbind.fill(data.frame(est_results$mean),data.frame(est_results$sd),data.frame(est_results$min),data.frame(est_results$X1Q),data.frame(est_results$Med),data.frame(est_results$X3Q),data.frame(est_results$max),data.frame(est_results$qciL),data.frame(est_results$qciU),data.frame(est_results$ntciL),data.frame(est_results$ntciU),data.frame(est_results$qciLen),data.frame(est_results$ntciLen))

    est_return$estimator <- rep(c("mean", "sd", "median", "IQR", "mad", "min", "max"), 13)
est_return$result <- unlist(lapply( c('mean', 'sd','min','X1Q', 'Med','X3Q','max' ,'qciL', 'qciU', 'ntciL', 'ntciU', "qciLen", "ntciLen"), rep,7))
}
write.table(est_return, paste0("../summ.est.", args), quote = FALSE)
                                        # figure out how to produce estimates results 



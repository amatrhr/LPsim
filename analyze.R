## R script containing helper functions to compare lasso selection stats under different conditions
## 
### Needs to do a couple of things:
#### 1. Report variances of the coefficients (or SD's, monotone trans, so w/e):
####       in the regular vs boot vs nested scenarios, can do for each sig one and the null ones
####
#### 2. Report five-number summaries of lambda values in each condition (avg 5NS for the REs; this is for every cell)
####
#### 3. Report average(SD) CI lengths for the coefficients in each scenario ^^^ This goes together with task #1


get_stat_fmests <- function(dataFrame, stat, variable, method, model, getstat){
    # inputs: the name of a dataframe, the statistic to get, and the variable (significant or null) to compute the statistic for the variable in the appropriate model (as the simple analyses have data laid out in a different way than the combined data are)
    ### Namely. the simple analysis data are 7 rows plus a header, where the first col gives the stat
    ### and the combined data are in 7 * 13 rows, with the last two cols giving the stats

    # Make a table to line up the models with which variables are important, and which unimportant
    models <- c("Fix", "ranks", "polyG")
    vtypes <- c("signal", "noise", "lambda")
    mv <-  expand.grid(models, vtypes)
    colnames(mv) <- c("model", "vtype")
    mv$loc <- list(1,1:5, 1:10,2:100,6:104,10:3000,101,105,3001)
    userow <- (mv$model == model & mv$vtype == variable)
    useloc <- unlist(mv[userow, 'loc'])
        useloc <- ifelse( (length(useloc) > 1),sample(useloc,1),useloc)
    if (method == "simple") {
    
        outstat <- dataFrame[match(stat,rownames(dataFrame)),useloc ]  # only one occurrence of the stat, so first match will be acceptable
    } else {
        statloc <- intersect(which(dataFrame$result == stat), which(dataFrame$estimator == getstat))
        outstat <- dataFrame[statloc, useloc] 
    }
    return(outstat)
}
args <- commandArgs(TRUE)
args <- "../combest.Fix.simple.0.0033.txt"
frameName <- read.table(args, header =FALSE, as.is = TRUE)

signif( sd(frameName[,1]), 3)
signif(2*qnorm(p=0.975)*sd(frameName[,1]),3)
signif(quantile(x = frameName[,1], probs = 0.975) - quantile(x = frameName[,1], probs = 0.025) , 3)

args <- "../summ.est.Fix.boot.0.0033.txt"
frameName <- read.table(args, header = TRUE, as.is = TRUE)

btestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="signal", method = "boot", model = "Fix", getstat = "mean")
signif(btestvars,3)
nbtestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="noise", method = "boot", model = "Fix", getstat = "mean")
signif(nbtestvars,3)


args <- "../summ.est.Fix.nested.0.0033.txt"
frameName <- read.table(args, header = TRUE, as.is = TRUE)

#lambdasnwo <- get_stat_fmests(dataFrame = frameName, stat = "min", variable ="lambda", method = "nested", model = "Fix", getstat = "mean")
testfivenum <- sapply(X = c('min','X1Q','Med','X3Q','max'), get_stat_fmests, dataFrame = frameName,  variable ="lambda", method = "nested", model = "Fix", getstat = "mean")

testvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="signal", method = "nested", model = "Fix", getstat = "mean")
signif(testvars,3)
ntestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="noise", method = "nested", model = "Fix", getstat = "mean")
signif(ntestvars,3)

## RANKS
args <- "../combest.ranks.simple.0.01.txt"
frameName <- read.table(args, header =FALSE, as.is = TRUE)
signif(fivenum(frameName[,105]), 3)

signif( sd(frameName[,1]), 3)
signif(2*qnorm(p=0.975)*sd(frameName[,1]),3)
signif(quantile(x = frameName[,1], probs = 0.975) - quantile(x = frameName[,1], probs = 0.025) , 3)

args <- "../summ.est.ranks.boot.0.01.txt"
frameName <- read.table(args, header = TRUE, as.is = TRUE)

btestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="signal", method = "boot", model = "ranks", getstat = "mean")
signif(btestvars,3)
nbtestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="noise", method = "boot", model = "ranks", getstat = "mean")
signif(nbtestvars,3) 


args <- "../summ.est.ranks.nested.0.01.txt"
frameName <- read.table(args, header = TRUE, as.is = TRUE)


rankfivenum <- sapply(X = c('min','X1Q','Med','X3Q','max'), get_stat_fmests, dataFrame = frameName,  variable ="lambda", method = "nested", model = "ranks", getstat = "mean")
signif(rankfivenum,3)
bnestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="signal", method = "nested", model = "ranks", getstat = "mean")
signif(bnestvars,3)
nbnestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="noise", method = "nested", model = "ranks", getstat = "mean")
signif(nbnestvars,3) 


## polyG
args <- "../combest.polyG.simple.0.01.txt"
frameName <- read.table(args, header =FALSE, as.is = TRUE)
signif(fivenum(frameName[,3001]), 3)

signif(fivenum(frameName[,105]), 3)

signif( sd(frameName[,1]), 3)
signif(2*qnorm(p=0.975)*sd(frameName[,1]),3)
signif(quantile(x = frameName[,1], probs = 0.975) - quantile(x = frameName[,1], probs = 0.025) , 3)




args <- "../summ.est.polyG.boot.0.01.txt"
frameName <- read.table(args, header = TRUE, as.is = TRUE)

btestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="signal", method = "boot", model = "polyG", getstat = "mean")
signif(btestvars,3)
nbtestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="noise", method = "boot", model = "polyG", getstat = "mean")
signif(nbtestvars,3) 

args <- "../summ.est.polyG.nested.0.01.txt"
frameName <- read.table(args, header = TRUE, as.is = TRUE)


pgfivenum <- sapply(X = c('min','X1Q','Med','X3Q','max'), get_stat_fmests, dataFrame = frameName,  variable ="lambda", method = "nested", model = "polyG", getstat = "mean")
signif(pgfivenum,3)

bnestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="signal", method = "nested", model = "polyG", getstat = "mean")
signif(bnestvars,3)
nbnestvars <- sapply(X = c('sd','qciLen','ntciLen'), get_stat_fmests, dataFrame = frameName,  variable ="noise", method = "nested", model = "polyG", getstat = "mean")
signif(nbnestvars,3) 

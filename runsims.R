## run script for simulations

## called using Rscipt
### Rscript runsims.R N p r model nsnps hsq method batchsize outfile1 outfile2
### N: number of observations
### p: number of simulated snps, under a non-polygenic model
### r: r-squared of association between important predictors and outcome under a non-polygenic model
### model:  "Fix"=single important predictor, "ranks"=five important predictors with fixed R2s (see paper), "polyG"=polygenic model
### nsnps: number of simulated snps under a polygenic model
### hsq: heritability of simulated polygenic trait
### method: method to select predictors; "simple"=select lambda via 5-fold cross-validation once on entire data, "boot"=bootstrap after cross-validated selection of lambda, "nested"=cross-validation nested within bootstrapping

## functions: 
## data_gen -- generates data according to the true model, gives achieved R^2, weights, true predictors, noise predictors, finalizes data 
## simple_reg -- simple regression for each column -- produces set of variable selections
## lasso_reg -- lasso regression to produce a set of variable selections
## boot_rep_las -- a single bootstrap replicate using the Lasso
## select_res -- convert coefficients from a single replicate to selection stats (TPR, FPR, TNR, VSP), using CIs if the replicate used bootstrapping
## est_res -- reorder coefficients from a single replicate so that the important predictors are first
## sim_init -- function to initialize a simulation object for one or several replicates under the same conditions
## boot_run -- wrapper function to perform bootstrap resampling for a single simulation replicate
## one_replication -- perform a single replicate of an initialized simulation



library(methods)
library(plyr) 
library(glmnet)
library(boot)
library(cvTools)
library(MASS)

data_gen <- function(N, psig, pnoise, model, rs , nsnps = 10000, hsq = 0.10){ 
    if(model == "ranks"){ psig <- 5}   
    ## model = 'solitary,' 'ranks,' 'polyG'
    if(model != "polyG"){
        
        totalP <- pnoise + psig
        fixdata <- matrix(rbinom(N*psig, size = 2, prob = 0.5), nrow = N, ncol = psig)
        noiseX <- matrix( rbinom(N * pnoise, size = 2, prob = 0.5), nrow = N, ncol = pnoise)
        fullX <- matrix(nrow = N, ncol = totalP)
        
        cols <- sample(totalP) 
        
        siglist <- cols[1:psig]
        fullX[ , siglist] <- fixdata[, 1:psig]
        
        nlist <- cols[(psig + 1):totalP ]
        fullX[, nlist] <- noiseX 
        
        rm(noiseX)    

        if (model == "ranks"){ 
            
            b <-  matrix(sqrt(2*c(0.01, 0.0067, 0.0033, 0.0022, 0.001)) , nrow = psig, ncol = 1) 
            errors <- sqrt(1 - sum(c(0.01, 0.0067, 0.0033, 0.0022, 0.001)))*rnorm(N) ## so that v(y) = 1 and coefs are readily interpretable

        } else {
            
            b <- matrix(sqrt(2*rs), ncol = 1) 
            errors <- sqrt(1 - sum(rs))*rnorm(N) ## so that v(y) = 1 and coefs are readily interpretable
            
        }  	

        Y <- fullX[,siglist]%*%b + errors
        achR2 <- summary(lm(Y~fixdata))$r.squared ## achieved R^2 for diagnostics
        
    } else {
        ## polyG models
        fullX <- matrix( rbinom( n = N * nsnps, size = 2, prob = 0.5), nrow = N, ncol = nsnps )
        b <- rnorm(nsnps, mean = 0, sd = sqrt(hsq/nsnps)) ## hsq prop of variance accounted for 
        b <- b[order(abs(b), decreasing = TRUE)]
        errors <- rnorm(N, mean = 0, sd = sqrt(1 - hsq) )
        lc <- scale(fullX)%*%b ## linear combination=>polygenic effect
        Y<- lc + errors

        psig <- floor(nsnps/300) ## Arbitrary definition of 'important' or 'significant' SNPs
        pnoise <- nsnps - psig
        totalP <- nsnps
        siglist <- 1:psig
        nlist <- (psig+1) : nsnps
        achR2 <- var(lc)/var(Y) ## achieved R^2 in polygenic model
    } 

    output <- list(Y = Y, fullX = fullX,  achR2 = achR2, siglist = siglist, nlist = nlist, p = totalP, psig = psig)
    
    return(output) 

}

lasso_reg <- function(X, Y,  nestL = FALSE, givenL = 0.01){
    if(nestL == TRUE){ ## Nested selection of lambda: resample each time
        lassomodel <- cv.glmnet(x = X, y = Y  )
        coefs <- coefficients(lassomodel, s = lassomodel$lambda.min)[-1]
        
        output <- matrix(c(coefs = coefs,lassomodel$lambda.min),nrow = 1)
        colnames(output) <- c(paste("coef", 1:length(coefs), sep = ""), "lambda")
        
    } else { ## non-nested: use givenL (given lambda value)
        lassomodel <-  glmnet(x = X, y = Y)
        output <- vector(mode = "list")
        output$coefs <- predict.glmnet(lassomodel, newx = X, s = givenL, type = "coef")[-1]
    }
    
    return(output)
    
}

boot_rep_las <- function(data, idx, nested = FALSE, givenL = 0.01 ){

    replicate <- lasso_reg(X = data[idx,-1], Y = data[idx,1], nestL = nested, givenL = givenL)
    if(nested == TRUE){ return(replicate) } else { return(replicate$coefs)}

}

select_res <- function(sim_example, model_results){
#### select_res() -> turns coefficients stats into selection stats
### coefs is the vector of coefficients in the simple case but the matrix of bootstrap results in the fixboot and nested cases
### lambda should be listed as one of the elements in the matrix of nested bootstrap results 
    
    basevec <- rep(0, (length(sim_example$modelData$siglist) + length(sim_example$modelData$nlist)))
    
### convert lists of indices to logical vectors (probably a clumsy way to do this)
    sigvec <- basevec
    sigvec[sim_example$modelData$siglist] <- 1 
    
    nvec <- basevec
    nvec[sim_example$modelData$nlist] <- 1 

    if(sim_example$method == "simple"){
### If given a single result, produce: FPR, TPR, FNR, VSP
        coefs <- model_results[-nrow(model_results), ] ### lambda is always to be the bottom row!
        selects <- as.numeric(coefs != 0) ### or abs() > 1e-9    
        TPs <- sum(selects*sigvec)
        TPR <- TPs/sum(sigvec + 1e-16)
        FPs <- sum(selects*nvec)
        FPR <- FPs/sum(nvec + 1e-16)
        TNR <- 1 - FPR
        
        VSP <- TPs /(TPs + FPs + 1e-16)
        
        output <- data.frame(list(VSP = VSP,  TPR = TPR, FPR = FPR, TNR = TNR))
        
    } else if (sim_example$method == "boot"){

        qselects <- as.numeric(model_results$qciL*model_results$qciU > 0) ### or abs() > 1e-9
        
        qTPs <- sum(qselects*sigvec)
        qTPR <- qTPs/sum(sigvec + 1e-16)
        
        qFPs <- sum(qselects*nvec)
        qFPR <- qFPs/sum(nvec + 1e-16)
        qTNR <- 1 - qFPR
        
        qVSP <- qTPs /(qTPs + qFPs + 1e-16)
###        qciLen <- coefs$qciU - coefs$qciL

        ntselects <- as.numeric(model_results$ntciL*model_results$ntciU > 0) ### or abs() > 1e-9

        ntTPs <- sum(ntselects*sigvec)
        ntTPR <- ntTPs/sum(sigvec + 1e-16)
        
        ntFPs <- sum(ntselects*nvec)
        ntFPR <- ntFPs/sum(nvec + 1e-16)
        ntTNR <- 1 - ntFPR
        
        ntVSP <- ntTPs /(ntTPs + ntFPs + 1e-16)
###       ntciLen <- coefs$ntciU - coefs$ntciL

        
        output <- data.frame(list(qVSP = qVSP,  qTPR = qTPR, qFPR = qFPR, qTNR = qTNR, ntVSP = ntVSP,  ntTPR = ntTPR, ntFPR = ntFPR, ntTNR = ntTNR))#, ntciLen = ntciLen , qciLen = qciLen)
        
### If given fixed-lambda bootstrap, produce: R/P & lengths for each kind of CI
    } else if (sim_example$method == "nested"){
        coefs <- model_results[-nrow(model_results), ] ### lambda is always to be the bottom row!
        qselects <- as.numeric(coefs$qciL*coefs$qciU > 0) ### or abs() > 1e-9

        qTPs <- sum(qselects*sigvec)
        qTPR <- qTPs/sum(sigvec + 1e-16)
        
        qFPs <- sum(qselects*nvec)
        qFPR <- qFPs/sum(nvec + 1e-16)
        qTNR <- 1 - qFPR
        
        qVSP <- qTPs /(qTPs + qFPs + 1e-16)
###        qciLen <- coefs$qciU - coefs$qciL
        
        ntselects <- as.numeric(coefs$ntciL*coefs$ntciU > 0) ### or abs() > 1e-9

        ntTPs <- sum(ntselects*sigvec)
        ntTPR <- ntTPs/sum(sigvec + 1e-16)
        
        ntFPs <- sum(ntselects*nvec)
        ntFPR <- ntFPs/sum(nvec + 1e-16)
        ntTNR <- 1 - ntFPR
        
        ntVSP <- ntTPs /(ntTPs + ntFPs + 1e-16)
        
        output <- data.frame(list(qVSP = qVSP,  qTPR = qTPR, qFPR = qFPR, qTNR = qTNR, ntVSP = ntVSP,  ntTPR = ntTPR, ntFPR = ntFPR, ntTNR = ntTNR))

    }

    return(output)
}

est_res <- function(sim_example, model_results){
##### est_res() -> unscramble coefficient stats
##### and save to a data structure of some kind
    if(sim_example$method == "simple"){

        results <- matrix(c(model_results[sim_example$modelData$siglist,], model_results[-sim_example$modelData$siglist,]), ncol = 1)
        colnames(results) <- "simpleEst"

    } else {
        results <- data.matrix(rbind(model_results[sim_example$modelData$siglist,], model_results[-sim_example$modelData$siglist,]))
    }
    return(data.frame(t(results))) ### transpose so that replications can be rbind() together
}



sim_init <- function(N = 2500, p = 100, r=0.5, model = 'Fix',nsnps = 10000, hsq = 0.10 ,seed = 4237, method = "simple"){

    pnoise <- p - 1 ## ugly hardcoding for scenarios in the paper

    if(model == "ranks"){ ## ugly hardcoding for scenarios in the paper
        p <- p + 5
        r <- 0.0232
    }
    
    
    modelData <- data_gen(N = N, psig = 1, pnoise = pnoise, model = model, rs = r, nsnps = nsnps, hsq = hsq)

    lasso <- cv.glmnet(x = modelData$fullX, y = modelData$Y)

    return(list(lasso = lasso, modelData = modelData, N = N, model = model, rs = r, nsnps = nsnps, hsq = hsq, nsnps = nsnps, p = p,  method = method))

}

boot_run <- function(sim_example){
    if(sim_example$method == "simple"){

        output <- data.frame(matrix(c(coef(sim_example$lasso, s = sim_example$lasso$lambda.min)[-1], sim_example$lasso$lambda.min), ncol = 1))
        
    } else if (sim_example$method == "boot"){

        nested <- FALSE
        bdata <- matrix(cbind(sim_example$modelData$Y, sim_example$modelData$fullX), nrow = length(sim_example$modelData$Y), ncol = (ncol(sim_example$modelData$fullX) + 1) )
        
        bootlas <- boot(bdata, boot_rep_las, R = 1000, nested = nested, givenL = sim_example$lasso$lambda.min)
        blasmat <- matrix( nrow = ncol(sim_example$modelData$fullX), ncol = 13)
        colnames(blasmat) <- c('mean', 'sd','min','1Q', 'Med','3Q','max' ,'qciL', 'qciU', 'ntciL', 'ntciU', "qciLen", "ntciLen")

        for (j in 1:nrow(blasmat)){
            cires <- boot.ci(bootlas, index = j, type = c("norm","perc"))
            blasmat[j, 1] <- mean(bootlas$t[,j])
            blasmat[j, 2] <- sd(bootlas$t[,j])
            blasmat[j, 3:7] <- fivenum(bootlas$t[,j])
            blasmat[j,8:9] <- cires$percent[,4:5]
            blasmat[j,10:11] <- cires$normal[,2:3]
            blasmat[j, 12] <- cires$percent[,5] - cires$percent[,4]
            blasmat[j, 13] <- cires$normal[,3] - cires$normal[,2]
        }

        output <- data.frame(blasmat)

    } else if (sim_example$method == "nested"){

        nested <- TRUE ##so ugly
        bdata <- matrix(cbind(sim_example$modelData$Y, sim_example$modelData$fullX), nrow = length(sim_example$modelData$Y), ncol = (ncol(sim_example$modelData$fullX) +1) )
        
        bootlas <- boot(bdata, boot_rep_las, R = 1000, nested = nested)
        blasmat <- matrix( nrow = (ncol(sim_example$modelData$fullX) + 1), ncol = 13)
        colnames(blasmat) <- c('mean', 'sd','min','1Q', 'Med','3Q','max' ,'qciL', 'qciU', 'ntciL', 'ntciU', "qciLen", "ntciLen")
        
        for (j in 1:nrow(blasmat)){
            cires <- boot.ci(bootlas, index = j, type = c("norm","perc"))
            blasmat[j, 1] <- mean(bootlas$t[,j])
            blasmat[j, 2] <- sd(bootlas$t[,j])
            blasmat[j, 3:7] <- fivenum(bootlas$t[,j])
            blasmat[j, 8:9] <- cires$percent[,4:5]
            blasmat[j, 10:11] <- cires$normal[,2:3]
            blasmat[j, 12] <- cires$percent[,5] - cires$percent[,4]
            blasmat[j, 13] <- cires$normal[,3] - cires$normal[,2]
        }
        
        output <- data.frame(blasmat)
        
    }

    return(output) 

}


one_replication <- function(conditions = list(N = 1000, p = 150, r = 0.16, model = "Fix", nsnps = 1000, hsq = 0.4, seed = 31, method = "simple")) {
    
    rep_init <- sim_init(N = conditions$N, p = conditions$p, r = conditions$r, model = conditions$model, nsnps = conditions$nsnps, hsq = conditions$hsq, seed = floor(1e7* runif(1)), method = conditions$method)
    rep_boot <- boot_run(sim_example = rep_init)
    rep_select <- select_res(sim_example = rep_init, model_results = rep_boot)
    rep_est <- est_res(sim_example = rep_init, model_results = rep_boot)
    
    return(list(est = rep_est, sel = rep_select))

}




args <- commandArgs(TRUE)

## debugging args <- list(500, 12, 0.5, "polyG", 1000, 0.5, "boot", 5, "../selcta.txt", "../ests.txt")

## Convert args to a list of numeric values
conditi = list(N = as.numeric(args[[1]]), p = as.numeric(args[[2]]), r = as.numeric(args[[3]]), model = args[[4]], nsnps = as.numeric(args[[5]]), hsq = as.numeric(args[[6]]), seed = floor(1e9*runif(1)), method = args[[7]])
batch <- args[[8]]
print(conditi)
set.seed(conditi$seed) ## problems if many jobs start on same procs at same times 

sel_obj <- paste("sel", 1:batch, sep = "")
est_obj <- paste("est", 1:batch, sep = "")


for (i in 1:batch){

    rep <- one_replication(conditions = conditi)
    assign(sel_obj[i], value = rep$sel)
    assign(est_obj[i], value = rep$est)
    
}



selects <- rbind.fill(lapply(sel_obj, get))

estimates <- rbind.fill(lapply(est_obj, get))

estimates[,ncol(estimates) + 1] <- rep(row.names(est1), batch)

ransims <- list(selects = selects, estimates = estimates) ## this object is here because it might be sensible to use this list rather than the names 


write.table(x = ransims$selects, file = args[[9]], row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(x = ransims$estimates, file = args[[10]], row.names = FALSE, col.names = FALSE, quote = FALSE)

# run script for simulations


# functions: 
# data_gen -- generates data according to the true model, gives achieved R^2, weights, true predictors, noise predictors, finalizes data 
# simple_reg -- simple regression for each column -- produces set of variable selections
# lasso_reg -- lasso regression to produce a set of variable selections
# boot_rep -- bootstrap replicates 
# res_fun
# Some other things to handle ranked predictors


# Other things that are required (loop over lambda and collect statistics on it)

require(glmnet)
require(boot)
require(cvTools)
require(MASS)

data_gen <- function(N, psig, pnoise, model, rs , nsnps = 10000, hsq = 0.10){ 
    if(model == "ranks"){ psig <- 5}   
    # model = 'solitary,' 'ranks,' 'polyG'
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
            errors <- sqrt(1 - sum(c(0.01, 0.0067, 0.0033, 0.0022, 0.001)))*rnorm(N) # so that v(y) = 1 and coefs are readily interpretable

        } else {
            
            b <- matrix(sqrt(2*rs), ncol = 1) 
            errors <- sqrt(1 - sum(rs))*rnorm(N) # so that v(y) = 1 and coefs are readily interpretable

        }  	

	Y <- fullX[,siglist]%*%b + errors
        achR2 <- summary(lm(Y~fixdata))$r.squared
   
    } else {
  # polyG models
        fullX <- matrix( rbinom( n = N * nsnps, size = 2, prob = 0.5), nrow = N, ncol = nsnps )
        b <- rnorm(nsnps, mean = 0, sd = sqrt(hsq/nsnps)) # hsq prop of variance accounted for by the 10000 loci with EAF 0.5; 
        b <- b[order(abs(b), decreasing = TRUE)]
        errors <- rnorm(N, mean = 0, sd = sqrt(1 - hsq) )
        lc <- scale(fullX)%*%b
        Y<- lc + errors
        siglist <- b[1:floor(nsnps/1000)] 
        nlist <- b[(floor(nsnps/1000)+1) : nsnps]
        achR2 <- var(lc)/var(Y) 
    } 

    output <- list(Y =Y, fullX = fullX,  achR2 = achR2, siglist = siglist, nlist = nlist)

    return(output) 

}

lasso_reg <- function(X, Y,  nestL = FALSE, givenL = 0.01){
    if(nestL == TRUE){ 
        lassomodel <- cv.glmnet(x = X, y = Y  )
        coefs <- coefficients(lassomodel, s = lassomodel$lambda.min)[-1]
     
        output <- matrix(c(coefs = coefs,lassomodel$lambda.min),nrow = 1)
        colnames(output) <- c(paste("coef", 1:length(coefs), sep = ""), "lambda")
    
    } else{
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

select_res <- function(siglist, nlist, coefs, method = "simple"){
### select_res() -> turns coefficients stats into selection stats
                                        # coefs is the vector of coefficients in the simple case but the matrix of bootstrap results in the fixboot and nested cases
    # lambda should be listed as one of the elements in the matrix of nested bootstrap results 
    
    basevec <- rep(0, (length(siglist) + length(nlist)))
    
       # convert lists of indices to logical vectors (probably a clumsy way to do this)
    sigvec <- basevec
    sigvec[siglist] <- 1 
        
    nvec <- basevec
    nvec[nlist] <- 1 

    if(method == "simple"){
    # If given a single result, produce: FPR, TPR, FNR, VSP
        selects <- as.numeric(coefs != 0) # or abs() > 1e-9    
        TPs <- sum(selects*sigvec)
        TPR <- TPs/sum(sigvec + 1e-16)
        FPs <- sum(selects*nvec)
        FPR <- FPs/sum(nvec + 1e-16)
        TNR <- 1 - FPR
    
        VSP <- TPs /(TPs + FPs + 1e-16)

        output <- list(VSP = VSP,  TPR = TPR, FPR = FPR, TNR = TNR)
   
    } else if (method == "boot"){

        qselects <- as.numeric(coefs$qciL*coefs$qciU > 0) # or abs() > 1e-9

        qTPs <- sum(qselects*sigvec)
        qTPR <- qTPs/sum(sigvec + 1e-16)
    
        qFPs <- sum(qselects*nvec)
        qFPR <- qFPs/sum(nvec + 1e-16)
        qTNR <- 1 - qFPR
    
        qVSP <- qTPs /(qTPs + qFPs + 1e-16)
#        qciLen <- coefs$qciU - coefs$qciL

        ntselects <- as.numeric(coefs$ntciL*coefs$ntciU > 0) # or abs() > 1e-9

        ntTPs <- sum(ntselects*sigvec)
        ntTPR <- ntTPs/sum(sigvec + 1e-16)
    
        ntFPs <- sum(ntselects*nvec)
        ntFPR <- ntFPs/sum(nvec + 1e-16)
        ntTNR <- 1 - ntFPR
    
        ntVSP <- ntTPs /(ntTPs + ntFPs + 1e-16)
 #       ntciLen <- coefs$ntciU - coefs$ntciL

        
        output <- list(qVSP = qVSP,  qTPR = qTPR, qFPR = qFPR, qTNR = qTNR, ntVSP = ntVSP,  ntTPR = ntTPR, ntFPR = ntFPR, ntTNR = ntTNR)#, ntciLen = ntciLen , qciLen = qciLen)
        
    # If given fixed-lambda bootstrap, produce: R/P & lengths for each kind of CI
    } else if (method == "nested"){
        qselects <- as.numeric(coefs$qciL*coefs$qciU > 0) # or abs() > 1e-9

        qTPs <- sum(qselects*sigvec)
        qTPR <- qTPs/sum(sigvec + 1e-16)
    
        qFPs <- sum(qselects*nvec)
        qFPR <- qFPs/sum(nvec + 1e-16)
        qTNR <- 1 - qFPR
    
        qVSP <- qTPs /(qTPs + qFPs + 1e-16)
#        qciLen <- coefs$qciU - coefs$qciL

        ntselects <- as.numeric(coefs$ntciL*coefs$ntciU > 0) # or abs() > 1e-9

        ntTPs <- sum(ntselects*sigvec)
        ntTPR <- ntTPs/sum(sigvec + 1e-16)
    
        ntFPs <- sum(ntselects*nvec)
        ntFPR <- ntFPs/sum(nvec + 1e-16)
        ntTNR <- 1 - ntFPR
    
        ntVSP <- ntTPs /(ntTPs + ntFPs + 1e-16)
 #       ntciLen <- coefs$ntciU - coefs$ntciL
  #      fivenumL <- fivenum(lambda)
        
        output <- list(qVSP = qVSP,  qTPR = qTPR, qFPR = qFPR, qTNR = qTNR, ntVSP = ntVSP,  ntTPR = ntTPR, ntFPR = ntFPR, ntTNR = ntTNR)#,, qciLen = qciLen ntciLen = ntciLen, mL = mL, sdL = sdL, fivenumL = fivenumL)
        
     # If given random-lambda bootstrap, produce: "" plus m/sd/5num lambda stats 
    }
                                        # Needs to only produce ROC-type selection stats and CI lengths.
    #                     1          2          3          15, 4, 5       6,7        8     9          10,11,12,13,14
    return(output)
}

est_res <- function(){
### est_res() -> unscramble coefficient stats
}

replications <- function(){}
=
sim_init <- function(N = 2500, p = 100, r=0.5, model = 'Fix',nsnps = 10000, hsq = 0.10 ,seed = 4237, method = "simple",  rep = 0){
### SPLIT THIS INTO A 'INITIAL MODEL' FUNCTION and a BOOTSTRAPPING FUNCTION
#                                        # problem is that 'nested' keywoed doesn't need to be reused
    pnoise <- p - 1

    if(model == "ranks"){
      p <- p + 5
      r <- 0.0232
    }
        
    set.seed(seed)
                                        
    modelData <- data_gen(N = N, psig = 1, pnoise = pnoise, model = model, rs = r, nsnps = nsnps, hsq = hsq)

    lasso <- cv.glmnet(x = modelData$fullX, y = modelData$Y)
    return(list(lasso = lasso, modelData = modelData, N = N, model = model, rs = r, nsnps = nsnps, hsq = hsq, nsnps = nsnps, p = p, rep = rep, method = method))
}

boot_run <- function(sim_example){
    if(sim_example$method == "simple"){

        output <- matrix(c(coef(sim_example$lasso, s = sim_example$lasso$lambda.min)[-1], sim_example$lasso$lambda.min))
        
    } else if (sim_example$method == "boot"){
     # any kind of bootstrapping, not just nested, so will need to add logic to handle such a case
     # make a new data matrix for the sake of bootstrapping?
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

        output <- blasmat

    } else if (sim_example$method == "nested"){
        nested <- TRUE
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
                
        output <- blasmat
        
    }

    return(output) # or write.table(output)

}

multi_rep_output <- function(output_list, method = "simple"){
    # Aggregate summary stats over all replications? How to do this?
    # It's wasteful to save to text file but the replicate function isn't going to handle a list with a big matrix element
    # Selection stats

    # Estimate summary stats

    # Yeah, the above code is not right! Need to return a list, p + 1 coefficeients, or similarly a p+1 matrix of bootstrap results; THEN aggregate all of those and later run res_fun() to get VSPs and such 
    # Should then, using some kind of cleverness, be able to call res.fun on the replicate() data
}

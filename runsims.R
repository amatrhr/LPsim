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

data_gen <- function(N, psig, pnoise, model, rs , nsnp = 10000, hsq = 0.10){ 
   
    # model = 'solitary,' 'ranks,' 'polyG'
    if(model != "polyG"){
        
        totalP <- pnoise + psig
        fixdata <- matrix(rbinom(N*psig, size = 2, prob = 0.5), nrow = N, ncol = psig)
        noiseX <- matrix( rbinom(N * pnoise, size = 2, prob = 0.5), nrow = N, ncol = pnoise)
    
        fullX <- matrix(nrow = N, ncol = totalP)
    
        cols <- sample(totalP) siglist <- cols[1:psig]

        fullX[ , siglist] <- X[, 1:psig]
    
        nlist <- cols[(psig + 1):(pnoise + psig) ]

        fullX[, nlist] <- noiseX 
        
        rm(noiseX)    

        if (model == "ranks"){ 
   
            b <-  matrix(c(0.01, 0.0067, 0.0033, 0.0022, 0.001) , nrow = psig, ncol = 1) 
            errors <- sqrt(1 - 0.25*sum(c(0.01, 0.0067, 0.0033, 0.0022, 0.001)))*rnorm(N) # so that v(y) = 1 and coefs are readily interpretable

        } else {
            b <- matrix(rs, ncol = 1) 
            errors <- sqrt(1 - 0.25*sum(rs))*rnorm(N) # so that v(y) = 1 and coefs are readily interpretable

        }  	

	Y <- fixdata%*%b + errors
    achR2 <- summary(lm(Y~fixdata))$r.squared
   
   } else {
  # polyG models
        polyGdata <- replicate( nsnps, rbinom( n = N, size = 2, prob = 0.5))
        b <- sqrt(hsq/(nsnps*0.25))*rnorm(nsnps) # hsq prop of variance accounted for by the 10000 loci with EAF 0.5; 
        b <- b[order(abs(b), decreasing = TRUE)]
        errors <-  (1 - sqrt(hsq/(nsnps*0.25)))*rnorm(nsnps)
        Y <- polyGdata%*%b + errors

  } 

    output <- list(fullX = fullX,  achR2 = achR2, siglist = siglist, nlist = nlist)

    return(output) 

}

lasso_reg <- function(X, Y,  nestL = FALSE){
    
    lassomodel <- cv.glmnet(x = X, y = Y  )
    coefs <- coefficients(lassomodel, s = lassomodel$lambda.min)[-1]
    lasselect <- as.numeric(coefs != 0)
   
    output <- list(select = lasselect, coefs = coefs, lambda = lassomodel$lambda.min, lambda1SE = lassomodel$lambda.1se)
    
    if(nestL == TRUE){ output <- list(lambda = lassomodel$lambda.min, lambda1SE = lassomodel$lambda.1se)} 
    
    return(output)

}

boot_rep_las <- function(data, idx, nested = FALSE ){

    replicate <- lasso_reg(X = data[idx,-1], Y = data[idx,1], nestL = nested)
       if(nested == TRUE){ return(replicate) } else { return(replicate$coefs)}

}



res_fun <- function(selects, siglist, nlist){
    basevec <- selects*0 
    
    # convert lists of indices to logical vectors (probably a clumsy way to do this)
    sigvec <- basevec
    sigvec[siglist] <- 1 
    
    nvec <- basevec
    nvec[nlist] <- 1 

    TPs <- sum(selects*sigvec)
    TPR <- TPs/sum(sigvec + 1e-16)
    
    FPs <- sum(selects*nvec)
    FPR <- FPs/sum(nvec + 1e-16)
    
    fVSP <- TPs /(TPs + FPs + 1e-16)

	# also need CI handling here
	# quantile
	# normal theory
    # lambda-distribution stuff
    output <- list(VSP = fVSP,  TPR = TPR, FPR = FPR)
    return(output)
}


diss_run <- function(N = 7, p = 10, r=0.5, model = 'Fix', seed = 4231, boots = FALSE,  rep = 0){
    pn <- p - 1
    

    set.seed(seed)
    modelData <- data_gen(N = N, psig = 1, wts = wt, model = model, phet = ph, rs = r)
    noisyData <- noise_gen(pnoise = pn, pnhet = pnh, psig = 1, X = modelData$X, phet = ph, Z = modelData$Z)

    simple <- simple_reg(X = noisyData$fullX, Y = modelData$Y, wts = (1/modelData$fweights))
    lasso <- lasso_reg(X = noisyData$fullX, Y = modelData$Y, wts = (1/modelData$fweights))
    outrows <- 2
    outcols <- 14

    output_array <- matrix(nrow = outrows, ncol = outcols)

    output_array[,1] <- N
    output_array[,2] <- p
    output_array[,3] <- r
    output_array[,4] <- wt
    output_array[,5] <- model
    output_array[,6] <- seed
    output_array[,7] <- modelData$achR
    output_array[1, 8:14] <-  matrix(c("simp", unlist(res_fun(simple$select, rep(0,length(c(noisyData$hetlist, noisyData$nhetlist))), siglist = noisyData$siglist, nlist = noisyData$nlist, hetlist = noisyData$hetlist, nhetlist = noisyData$nhetlist))))

    output_array[2, 8:14] <-  c("las", unlist(res_fun(lasso$select, rep(0,length(c(noisyData$hetlist, noisyData$nhetlist))), siglist = noisyData$siglist, nlist = noisyData$nlist, hetlist = noisyData$hetlist, nhetlist = noisyData$nhetlist)))

                                
colnames(output_array) <- c("N", "p", "R2", "wts", "model", "seed","achR","meth","vsp", "tpr", "fpr", "hvsp", "htpr", "hfpr")

    if (model == "Hetero"){
        simH <- simple_het(X = noisyData$fullX, Y = modelData$Y, Z = noisyData$fullZ, wts = (1/modelData$fweights), hetlist = c(noisyData$hetlist, noisyData$nhetlist))
        simH_res<-  unlist(res_fun(simple$select,simH$select, siglist = noisyData$siglist, nlist = noisyData$nlist, hetlist = noisyData$hetlist, nhetlist = noisyData$nhetlist))
        output_array[1, 12:14] <- simH_res[4:6]
        lasH <- lasso_het(X = noisyData$fullX, Y = modelData$Y, Z = noisyData$fullZ, wts = (1/modelData$fweights), hetlist = c(noisyData$hetlist, noisyData$nhetlist))
        lasH_res<-  unlist(res_fun(lasso$select, lasH$select, siglist = noisyData$siglist, nlist = noisyData$nlist, hetlist = noisyData$hetlist, nhetlist = noisyData$nhetlist))
        output_array[2, 12:14] <- lasH_res[4:6]
    }


    
    if (scad == TRUE){
        scad <- scad_reg(X = noisyData$fullX, Y = modelData$Y, wts = (1/modelData$fweights))
        scadout <-  c(output_array[1,1:7], unlist(c("scad", res_fun(scad$select, rep(0,length(c(noisyData$hetlist, noisyData$nhetlist))),siglist = noisyData$siglist, nlist = noisyData$nlist, hetlist = noisyData$hetlist, nhetlist = noisyData$nhetlist))))
        output_array <- rbind(output_array, scadout)
    }

    if (boots == TRUE){
        bdata <- matrix(cbind(modelData$Y, noisyData$fullX), nrow = length(modelData$Y), ncol = (ncol(noisyData$fullX) +1) )

        T <- max(ncol(noisyData$fullZ),0)

        bytpe <- "ordinary"
        if (nrow(bdata < 10)) {
            btype <- "balanced"
        }
        bootsimp <- boot(bdata, boot_rep_simp, sim = btype, R = 1000, wts = (1/modelData$fweights), hetlist=c(noisyData$hetlist, noisyData$nhetlist), Z = noisyData$fullZ)
        bsimpmat <- matrix(nrow = (ncol(noisyData$fullX) + is.matrix(noisyData$fullZ)*T), ncol = 3)
        for (j in 1:nrow(bsimpmat)){
            bsimpmat[j,1:2] <- boot.ci(bootsimp,index = j, type = "perc")$percent[4:5]
            bsimpmat[j,3] <- (sign(bsimpmat[j,1])*sign(bsimpmat[j,2]) > 0)
        }
        QQ <- 0
        if(model =="Hetero"){

            QQ <- bsimpmat[(ncol(noisyData$fullX) + 1):nrow(bsimpmat),3] 
            
        }
        
        bsout <-  c(output_array[1,1:7], c("bsimp", unlist(res_fun(bsimpmat[1:ncol(noisyData$fullX),3], QQ, siglist = noisyData$siglist, nlist = noisyData$nlist, hetlist = noisyData$hetlist, nhetlist = noisyData$nhetlist))))

        bootlas <- boot(bdata, boot_rep_las, R = 1000, wts = (1/modelData$fweights), hetlist=c(noisyData$hetlist, noisyData$nhetlist), Z = noisyData$fullZ)
        blasmat <- matrix( nrow = (ncol(noisyData$fullX) + is.matrix(noisyData$fullZ)*T), ncol = 3)
        for (j in 1:nrow(blasmat)){
            blasmat[j,1:2] <- boot.ci(bootlas, index = j, type = "perc")$percent[4:5]
            blasmat[j,3] <- (sign(blasmat[j,1])*sign(blasmat[j,2]) > 0)
        }
        blasout <-  c(output_array[1,1:7], c("blas", unlist(res_fun(blasmat[1:ncol(noisyData$fullX),3], QQ,siglist = noisyData$siglist, nlist = noisyData$nlist, hetlist = noisyData$hetlist, nhetlist = noisyData$nhetlist))))

        output_array <- rbind(output_array, bsout, blasout)
    }

    write.table(output_array, paste("ds", N, ".", p,".", r,".", wt, ".",model,".", "rep", rep, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}


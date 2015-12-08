# run script for dissertation simulations


# functions: 
# data_gen -- generates data according to the true model, gives achieved R^2, weights, true predictors, 
# noise_gen -- generates noise predictors, finalizes data 
# simple_reg -- simple regression for each column -- produces set of variable selections
# lasso_reg -- lasso regression to produce a set of variable selections
# scad_reg -- scad regression for a set of variable selections
# het_lasso -- lasso regression to select linear contributors to heteroskedasticity
# het_scad -- scad regression to select linear contributors to heteroskedasticity
# boot_rep -- bootstrap replicates 
# res_fun
# Some other things to handle ranked predictors

require(glmnet)
require(boot)
require(cvTools)
require(MASS)
data_gen <- function(N, psig, wts, rs, model, phet){
    fixdata <- matrix(rnorm(N*psig), nrow = N, ncol = psig)
    hetpred <- matrix(0, nrow = N, ncol = 1)
    weights <- rep(1,N)
    X <- fixdata
    Z <- 0
    if (phet > 0){
        Z <- matrix(runif(N*phet), nrow = N, ncol = phet)%*%matrix(2, nrow = phet, ncol = 1)
        
    }

    if (wts == "Skw"){
        weights <- matrix(rexp(N, 1), nrow = N)
    }

    fweights <- weights
    errors <- rnorm(N, 0, sqrt(weights)) ### sqrt wts gives the errors 

    
    if (model != "Fix"){
        q <- 0.80
        random_effects <- (1-q)*rnorm( N, 0, sqrt(mean(weights)) ) + (q)*rnorm(N, 0, sqrt(Z)) 
        errors <- errors + random_effects ####  = 1.04 sigma in RE case, 1.04 sigma + 1.28u
        weights <- weights + matrix( mean(weights), nrow = N) + Z

    }

    Rmult <- rs/(1-rs)
    W <- diag(N)
    diag(W) <- 1/weights
   
    Q <- matrix(1,ncol = psig)%*%t(fixdata) %*% W %*%fixdata%*%matrix(1,nrow = psig)
    b <-  matrix(sqrt(Rmult * N / sum(diag(Q))), nrow = psig, ncol = 1)
    Y <- fixdata%*%b + errors

    output <- list(X=X, Y=Y, tweights=weights, fweights = fweights, achR = summary(lm(Y~fixdata, weights = 1/weights))$r.squared, Z=Z)
    return(output)
}

noise_gen <- function(pnoise, pnhet, psig, phet, X, Z){
    N <- nrow(X)

    noiseX <- matrix( rnorm(N * pnoise), nrow = N, ncol = pnoise)
    
    totalP <- pnoise + psig
    totalH <- phet + pnhet
    
    fullX <- matrix(nrow = N, ncol = totalP)
    fullZ <- 0
    hetlist <- NA
    nhetlist <- NA
    
    cols <- sample(totalP)
    siglist <- cols[1:psig]

    fullX[ , siglist] <- X[, 1:psig]
    
    nlist <- cols[(psig + 1):(pnoise + psig) ]

    fullX[, nlist] <- noiseX
    
    if (phet > 0){
        fullZ <- matrix(nrow = N, ncol = totalH)
        zcols <- sample(totalH)
        hetlist <- zcols[1:phet]
        fullZ[,hetlist] <- Z

        if  (pnhet > 0){
            noiseP <- matrix(runif(N*pnhet), nrow = N, ncol = pnhet)
            nhetlist <- zcols[(phet + 1):totalH]
            fullZ[, nhetlist] <- noiseP
        }
    }
    output <- list(fullX = fullX, fullZ = fullZ, siglist = siglist,  hetlist = hetlist, nlist = nlist,  nhetlist = nhetlist )

    return(output)
}

simple_reg <- function(X, Y, wts){
    P <- ncol(X)
    fixselect <- vector(mode = "numeric", length = P)
    coefs <- vector(mode = "numeric", length = P)
    
    for (j in 1: P){
        templm <- lm(Y~X[,j], weights = wts)
        tempp <- summary(templm)$coefficients[2,"Pr(>|t|)"]
        fixselect[j] <- tempp <= 0.05
        coefs[j] <- summary(templm)$coefficients[2,"Estimate"]
    }

    output <- list(select = fixselect, coefs = coefs)
    return(output)
}

lasso_reg <- function(X, Y, wts, nestL = FALSE){

    

    if(nestL == TRUE){

    } else{
      #
      lassomodel <- cv.glmnet(x = X, y = Y, weights = wts)
      coefs <- coefficients(lassomodel, s = lambda.min)[-1]
    }
    
    

    lasselect <- as.numeric(coefs != 0)
    output <- list(select = lasselect, coefs = coefs)
    return(output)

}



scad_reg <- function(X, Y,  wts){

    folds <- cvFolds(n = nrow(X), K = 5, type = "random")
    lambdarange <- seq(1e-6, max(abs(cov(X,Y))), length.out = 100)
    cvMatrix <- matrix(nrow = length(lambdarange), ncol =8)

    for (j in 1:5){
        train <- folds$subsets[folds$which == j]
        test <-  folds$subsets[folds$which != j]
        xtrain <- X[train, ]
        ytrain <- Y[train ]
        xtest <- X[test, ]
        ytest <- Y[test ]
        for (l in 1:100){
        trainmodel <- fullscadglm(x = xtrain, y = ytrain, weight = wts[train], family = gaussian(), lambda = lambdarange[l])
        testresid <- cbind(1,xtest)%*%trainmodel - ytest
        mspe <- mean(testresid^2)
        cvMatrix[l,(j + 1)] <- mspe
        cvMatrix[l,1] <- lambdarange[l]
    }
    }

    cvMatrix[,7] <- rowMeans(cvMatrix[,2:6])
    cvMatrix[,8] <- apply(cvMatrix[,2:6], 1, sd)

    lambda.min <- cvMatrix[which.min(cvMatrix[,7]), 1]
    scadmodel <- fullscadglm(x = X, y = Y, weight = wts, family = gaussian(), lambda = lambda.min)
    coefs <- scadmodel[-1]
    scadselect <- as.numeric(coefs != 0)
    # add cross-validation # cvtools
    
    output <- list(select = scadselect, coefs = coefs)
    return(output)
    
}




simple_het <- function(X, Y, Z, wts, hetlist){
    # read amemiya paper again
    P <- length(hetlist)
    fhetselect <- vector(mode = "numeric", length = P)
    coefs <- vector(mode = "numeric", length = P)


    fixR <-    Y - X%*%ginv(t(X)%*%X)%*%t(X)%*%Y
    fixR2 <- fixR^2
    for (j in 1:P){
        templm <- lm(fixR2~Z[,j])
        D <- (templm$fitted.values)^2
        newWt <- 1/D
        tempwlm <- lm(fixR2~Z[,j], weights = newWt)
        tempp <- summary(tempwlm)$coefficients[2,"Pr(>|t|)"]
        fhetselect[j] <- tempp <= 0.05
        coefs[j] <- summary(tempwlm)$coefficients[2,"Estimate"]
    }

    output <- list(select = fhetselect, coefs = coefs)
    return(output)
}

lasso_het <- function(X, Y, Z, wts, hetlist){
    # read amemiya paper again
    P <- length(hetlist)
    fhetselect <- vector(mode = "numeric", length = P)
    coefs <- vector(mode = "numeric", length = P)

    fixR <-    Y - X%*%ginv(t(X)%*%X)%*%t(X)%*%Y
    fixR2 <- fixR^2
    tempglm <- cv.glmnet(y = fixR2, x = Z)
    D <- (predict(tempglm, newx = Z, s = tempglm$lambda.min))^2
    newWt <- 1/D
    tempwlm <- cv.glmnet(y = fixR2, x = Z, weights = newWt)
    coefs <- coefficients(tempwlm, s = tempwlm$lambda.min)[-1]
    fhetselect <- as.numeric(coefs !=0)

    output <- list(select = fhetselect, coefs = coefs)
    return(output)
}


scad_het <- function(X, Y, Z, wts, hetlist){
    P <- length(hetlist)
    fhetselect <- vector(mode = "numeric", length = P)
    coefs <- vector(mode = "numeric", length = P)

    fixR <- lm(Y~X, weights = wts)$residuals
    fixR2 <- fixR^2
    tempscad <- scad_reg(Y = fixR2, X = Z)
    D <- (hetX%*%tempscad$coefs)^2
    newWt <- 1/D
    tempwlm <- scad_reg(Y = fixR2, X = Z, wts = newWt)
    coefs <- tempwlm$coefs
    fhetselect <- as.numeric(coefs !=0)

    output <- list(select = fhetselect, coefs = coefs)
    return(output)

}


boot_rep_simp <- function(data, idx, wts, hetlist = FALSE, Z = NA){
    if (!is.na(hetlist)){
        replicate <- simple_reg(X = data[idx, -1], Y = data[idx,1], wts = wts[idx])
        hreplicate <- simple_het(X = data[idx,-1], Y = data[idx,1], Z=Z, wts = wts[idx], hetlist = hetlist)
    } else {replicate <- simple_reg(X = data[idx,-1], Y = data[idx,1], wts = wts[idx])
            hreplicate <- list(NA, coefs=NA)
        }

    return(c(replicate$coefs, hreplicate$coefs))
}


boot_rep_las <- function(data, idx, wts, hetlist = FALSE, Z = NA){

      if (!is.na(hetlist)){
        replicate <- lasso_reg(X = data[idx, -1], Y = data[idx,1], wts = wts[idx])
        hreplicate <- lasso_het(X = data[idx,-1], Y = data[idx,1], Z=Z, wts = wts[idx], hetlist = hetlist)
    } else {replicate <- lasso_reg(X = data[idx,-1], Y = data[idx,1], wts = wts[idx])
            hreplicate <- list(NA, coefs=NA)
        }

    return(c(replicate$coefs, hreplicate$coefs))

}


boot_rep_scad <- function(data, idx, wts, hetlist = FALSE, Z = NA){

    if (!is.na(hetlist)){
        replicate <- scad_reg(X = data[idx, -1], Y = data[idx,1], wts = wts[idx])
        hreplicate <- scad_het(X = data[idx,-1], Y = data[idx,1], Z=Z, wts = wts[idx], hetlist = hetlist)
    } else {replicate <- scad_reg(X = data[idx,-1], Y = data[idx,1], wts = wts[idx])
            hreplicate <- list(NA, coefs=NA)
        }

    return(c(replicate$coefs, hreplicate$coefs))

}



res_fun <- function(selects, hselects, siglist, nlist, hetlist, nhetlist){
    basevec <- selects*0 
    hbasevec <- hselects*0 
    # convert lists of indices to logical vectors 
    sigvec <- basevec
    sigvec[siglist] <- 1 
    
    nvec <- basevec
    nvec[nlist] <- 1 

    hetvec <- hbasevec
    hetvec[hetlist] <- 1 

    nhetvec <- hbasevec
    nhetvec[nhetlist] <- 1 


    TPs <- sum(selects*sigvec)
    TPR <- TPs/sum(sigvec + 1e-16)
    
    FPs <- sum(selects*nvec)
    FPR <- FPs/sum(nvec + 1e-16)
    
    fVSP <- TPs /(TPs + FPs + 1e-16)

    hTP <- sum(hselects*hetvec)
    hTPR <- hTP/sum(hetvec + 1e-16)
    
    hFP <- sum(hselects*nhetvec)
    hFPR <- hFP/sum(nhetvec + 1e-16)

    hVSP <- hTP/(hTP + hFP + 1e-16)
    output <- list(VSP = fVSP,  TPR = TPR, FPR = FPR, hVSP = hVSP, hTPR = hTPR, hFPR = hFPR)
    return(output)
}


diss_run <- function(N = 7, p = 10, r=0.5, wt='Eq', model = 'Fix', seed = 4231, boots = FALSE, scad = FALSE, rep = 0){
    pn <- p - 1
    ph <- 0
    pnh <- 0
    
    if (model == "Hetero"){
        ph <- 1
        pnh <- 1
        pn <- p - 3
    }

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


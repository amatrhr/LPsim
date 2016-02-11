# Interaction Analysis of SNP effects on PaiBOR sumscore 

### input files:
#     top 125 SNPs, and covariates in the NTR sample (random selection of MZ twins, N = 5802, no IDs)--ntr_noIDs_data_matrix_Apr28.txt
#     top 125 SNPs, and covariates in the NESDA sample (random selection of MZ twins, N = 1323, no IDs)--nes_noIDs_data_matrix_Apr28.txt

### output files:
#     text files---
#              NTR_includedSTATS.txt: included predictors, ranked by coefficient of variation, NTR sample (with full-sample coefficient estimate, and fixed/random lambda values bootstrap mean estimate, bootstrap SE estimate, and percentile bootstrap CI)
#              NESDA_includedSTATS.txt: included predictors, ranked by coefficient of variation, NESDA sample (with full-sample coefficient estimate, and fixed/random lambda values bootstrap mean estimate, bootstrap SE estimate, and percentile bootstrap CI)
#              COMBINED_includedSTATS.txt: included predictors, ranked by coefficient of variation, combined sample (with full-sample coefficient estimate, and fixed/random lambda values bootstrap mean estimate, bootstrap SE estimate, and percentile bootstrap CI)
#
#
#               NTR_excludedSTATS.txt: excluded predictors, ranked by bootstrap SE estimate, NTR sample (with full-sample coefficient estimate, and fixed/random lambda values bootstrap mean estimate, bootstrap SE estimate, and percentile bootstrap CI)
#               NESDA_excludedSTATS.txt: excluded predictors, ranked by bootstrap SE estimate, NESDA sample (with full-sample coefficient estimate, and fixed/random lambda values bootstrap mean estimate, bootstrap SE estimate, and percentile bootstrap CI)
#               COMBINED_excludedSTATS.txt: excluded predictors, ranked by bootstrap SE estimate, combined sample (with full-sample coefficient estimate, and fixed/random lambda values bootstrap mean estimate, bootstrap SE estimate, and percentile bootstrap CI)

#     png graphics files---
#               bootstrap_moments.png: plot of fixed/random bootstrap mean estimates against bootstrap SE estimates, NTR, NESDA, and combined on same scale (2 x 3 array of scatter plots)

### it does:
#     1. Read in the 132 predictors and 7125 individuals from the two sample files to ntr.data.M and nesda.data.M
#     2. Calculate interaction terms for the 125 SNPs, within each sample, storing results in ntr.int.M, and nesda.int.M, cbind to the SNPs, getting ntr.pred.M, nesda.pred.M, combined.pred.M
#     3. Calculate the outcome as the residual of PAIBOR average score on age, sex, age*sex interaction, and PC1_NL, represening the geographic origin of each individual. Store this residual in ntr.outcome, nesda.outcome, and overall.outcome
#     4. Fit a lasso model to each of the three sets of predictors and outcomes
#     5. Bootstrap model fit using the regularization parameter selected in step 4
#       5.1. Using 1000 bootstrap replications, estimate the mean, SE, CV, inverted-t CI, and percentile CI of the coefficient for each predictor

#     6. Bootstrap model fit, re-selecting the regularization parameter in each resample
#       6.1. Using 1000 bootstrap replications, estimate the mean, SE, CV, inverted-t CI, and percentile CI of the coefficient for each predictor
#     7. Generate output files
###################################################################################################
library(glmnet)
library(doMC)
#registerDoMC(8)
registerDoMC()
library(ggplot2)
source("percentile_medians.R")
source("cor_resolve.R")


B <- 1000# Number of Bootstrap replications


#     1. Read in the 132 predictors and 7125 individuals from the two sample files to ntr.data.M and nesda.data.M
ntr.data.M <- read.table("ntr_noIDs_data_matrix_Apr28.txt", header = TRUE, as.is = TRUE)
ntr.data.M <- na.omit(ntr.data.M)
ntr.data.M <- corResolve(ntr.data.M, 0.60) # clear out SNPs that are highly correlated


## List the columns containing SNPs
snp.columns <- grep("X", colnames(ntr.data.M)) 

ntr.SNPs.M <- ntr.data.M[, snp.columns]
snp.names <- colnames(ntr.SNPs.M)

nesda.data.M <- read.table("nes_noIDs_data_matrix_Apr28.txt", header = TRUE, as.is = TRUE)
nesda.removSNPs <- setdiff( colnames(nesda.data.M), colnames(ntr.data.M) )
nesda.removCols <- match(nesda.removSNPs, colnames(nesda.data.M))

nesda.data.M <- nesda.data.M[, -nesda.removCols]

combined.data.M <- rbind(ntr.data.M, nesda.data.M)

nesda.SNPs.M <- nesda.data.M[, snp.columns]


# list all pairs of SNP columns, each in a row of a matrix
combination.M <- expand.grid(1:length(snp.columns), 1:length(snp.columns))

# remove rows where a column is paired with itself or duplicate pairings
# self-pairings
selfPairs <- combination.M[,1] == combination.M[,2]

#duplicate pairings
# 3 steps:
## 1. sort the rows of the combination matrix so that rows having the same pairs are in the same order using apply(X, 1, sort)
###    e.g. [1,5] and [5,1] are now both [1,5]
## 2. transpose this because apply returns each sorted row as a column
## 3. The function duplicated() produces a logical vector that is TRUE for rows that are duplicates of others


dupePairs <- duplicated( t( apply(combination.M, 1, sort ) ) )

excludeRows <- selfPairs + dupePairs # logical OR 

combination.M <- combination.M[!excludeRows, ]

print(dim(combination.M)) # should be 2926 by 2

# column names for the interaction vectors
names.vector <- vector(mode = "character", length = nrow(combination.M))

ntr.int.M <- matrix(NA, nrow = nrow(ntr.SNPs.M), ncol = nrow(combination.M))
nesda.int.M <- matrix(NA, nrow = nrow(nesda.SNPs.M), ncol = nrow(combination.M))

for (pair in 1:nrow(combination.M)){

  ntr.int.M[, pair] <- scale( ntr.SNPs.M[, combination.M[pair, 1]], scale = FALSE, center = TRUE ) * scale( ntr.SNPs.M[, combination.M[pair, 2]], scale = FALSE, center = TRUE )
  nesda.int.M[, pair] <- scale( nesda.SNPs.M[, combination.M[pair, 1]], scale = FALSE, center = TRUE ) * scale( nesda.SNPs.M[, combination.M[pair, 2]], scale = FALSE, center = TRUE )
  names.vector[pair] <- paste( snp.names[ combination.M[pair, 1] ], "X", snp.names[ combination.M[pair, 2] ], sep ="" )
  
}

colnames(ntr.int.M) <- names.vector
colnames(nesda.int.M) <- names.vector

#     2. Calculate interaction terms for the 125 SNPs, within each sample, storing results in ntr.int.M, and nesda.int.M, cbind to the SNPs, getting ntr.pred.M, nesda.pred.M, combined.pred.M

ntr.pred.M <- data.matrix(cbind(ntr.SNPs.M, ntr.int.M))
nesda.pred.M <- data.matrix(cbind(nesda.SNPs.M, nesda.int.M))
combined.pred.M <- rbind(ntr.pred.M, nesda.pred.M)
print(date() )
# clear up some space 
rm(list = c("ntr.SNPs.M", "ntr.int.M", "nesda.SNPs.M", "nesda.int.M"))

#     3. Calculate the outcome as the residual of PAIBOR average score on age, sex, age*sex interaction, and PC1_NL, represening the geographic origin of each individual. Store this residual in ntr.outcome, nesda.outcome, and overall.outcome
ntr.outcome <- lm(tot.m ~ sex + age7 + sex*age7 + PC1_NL, data = ntr.data.M)$residuals
nesda.outcome <- lm(tot.m ~ sex + age7 + sex*age7 + PC1_NL, data = nesda.data.M)$residuals
combined.outcome <- lm(tot.m ~ sex + age7 + sex*age7 + PC1_NL, data = combined.data.M)$residuals
print(date())
#     4. Fit a lasso model to each of the three sets of predictors and outcomes
ntr.fullsamp.lasso <- cv.glmnet( x = ntr.pred.M, y = ntr.outcome) # cross-validation
ntr.fullsamp.coefs <- coef ( ntr.fullsamp.lasso, s = ntr.fullsamp.lasso$lambda.min ) [ -1 ] # get coefficients associated with minimum lambda in entire sample, except intercept
print(date())
nesda.fullsamp.lasso <- cv.glmnet( x = nesda.pred.M, y = nesda.outcome) # cross-validation
nesda.fullsamp.coefs <- coef ( nesda.fullsamp.lasso, s = nesda.fullsamp.lasso$lambda.min ) [ -1 ] # get coefficients associated with minimum lambda in entire sample, except intercept
print(date())
combined.fullsamp.lasso <- cv.glmnet( x = combined.pred.M, y = combined.outcome ) # cross-validation 
combined.fullsamp.coefs <- coef ( combined.fullsamp.lasso, s = combined.fullsamp.lasso$lambda.min ) [ -1 ] # get coefficients associated with minimum lambda in entire sample, except intercept 

#     5. Bootstrap model fit using the regularization parameter selected in step 4
#       5.1. Using 1000 bootstrap replications, estimate the mean, SE, CV, inverted-t CI, and percentile CI of the coefficient for each predictor
print(date())
#### NTR BOOTSTRAPPING ####
ntr.boot.fix.out <- foreach( i = 1:B, .combine = 'cbind', .packages = 'glmnet', .noexport = c( "boot.fix.samples", "boot.fix.X", "boot.fix.Y" ) ) %dopar% {
  
  boot.fix.samples <- sample( nrow( ntr.pred.M ), replace = TRUE )
  boot.fix.X <- ntr.pred.M[ boot.fix.samples, ]
  boot.fix.Y <- ntr.outcome[ boot.fix.samples ]
  
  boot.fix.coefs <- coef( glmnet(x = boot.fix.X, y = boot.fix.Y, lambda = ntr.fullsamp.lasso$lambda.min, alpha = 1  )) [ -1 ] 

}
print(date())
# output matrix of bootstrap coefficient estimates and gzip -f it 
write.table( ntr.boot.fix.out, "ntr_out_boot_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep = "\t" )
system( "gzip -f ntr_out_boot_fix.txt" )


# Percentile bootstrap bounds
ntr.fix.lbs <- apply( ntr.boot.fix.out, 1, quantile, 0.025 ) # lower bound
ntr.fix.ubs <- apply( ntr.boot.fix.out, 1, quantile, 0.975 ) # upper bound

# coverage rate 
ntr.fix.coverage <- as.numeric ( sign ( ntr.fix.lbs * ntr.fix.ubs  ) <=  0 ) # does the percentile interval contain 0

# bootstrap estimate of the mean
ntr.fix.mean <- apply( ntr.boot.fix.out, 1, mean ) 

# bootstrap estimate of the sd
ntr.fix.sd <- apply( ntr.boot.fix.out, 1, sd )

ntr.fix.cv <- ntr.fix.sd/abs(ntr.fix.mean)

# output the data and gzip -f it
write.table( ntr.fix.lbs, "ntr_out_lb_boot_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_lb_boot_fix.txt")

write.table( ntr.fix.ubs, "ntr_out_ub_boot_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_ub_boot_fix.txt")

write.table( ntr.fix.coverage, "ntr_out_C_boot_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t") 
system("gzip -f ntr_out_C_boot_fix.txt")

write.table( ntr.fix.mean, "ntr_out_mn_boot_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_mn_boot_fix.txt")

write.table( ntr.fix.sd, "ntr_out_sd_boot_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_sd_boot_fix.txt")

write.table( ntr.fix.cv, "ntr_out_cv_boot_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_cv_boot_fix.txt")



# Naive T interval bootstrap bounds 
ntr.t.boot.fix.lb <- ntr.fullsamp.coefs + qnorm( .025 ) * ntr.fix.sd # lower bound
ntr.t.boot.fix.ub <- ntr.fullsamp.coefs + qnorm( .975 ) * ntr.fix.sd # upper bound

# coverage rate 
ntr.t.boot.fix.cover <- as.numeric ( sign ( ntr.t.boot.fix.lb * ntr.t.boot.fix.ub  ) <=  0 ) 

# output
write.table( ntr.t.boot.fix.lb, "ntr_out_lb_NTB_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_lb_NTB_fix.txt")

write.table( ntr.t.boot.fix.ub, "ntr_out_ub_NTB_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_ub_NTB_fix.txt")

write.table( ntr.t.boot.fix.cover, "ntr_out_C_NTB_fix.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_C_NTB_fix.txt")

#### NESDA BOOTSTRAPPING ####
nesda.boot.fix.out <- foreach( i = 1:B, .combine = 'cbind', .packages = 'glmnet', .noexport = c( "boot.fix.samples", "boot.fix.X", "boot.fix.Y" ) ) %dopar% {
  
  boot.fix.samples <- sample( nrow( nesda.pred.M ), replace = TRUE )
  boot.fix.X <- nesda.pred.M[ boot.fix.samples, ]
  boot.fix.Y <- nesda.outcome[ boot.fix.samples ]
  
  boot.fix.coefs <- coef( glmnet(x = boot.fix.X, y = boot.fix.Y, lambda = nesda.fullsamp.lasso$lambda.min, alpha = 1  )) [ -1 ] 

}
print(date())
# output matrix of bootstrap coefficient estimates and gzip -f it 
write.table( nesda.boot.fix.out, "nesda_out_boot_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep = "\t" )
system( "gzip -f nesda_out_boot_fix.txt" )


# Percentile bootstrap bounds
nesda.fix.lbs <- apply( nesda.boot.fix.out, 1, quantile, 0.025 ) # lower bound
nesda.fix.ubs <- apply( nesda.boot.fix.out, 1, quantile, 0.975 ) # upper bound

# coverage rate 
nesda.fix.coverage <- as.numeric ( sign ( nesda.fix.lbs * nesda.fix.ubs  ) <=  0 ) # does the percentile interval contain 0

# bootstrap estimate of the mean
nesda.fix.mean <- apply( nesda.boot.fix.out, 1, mean ) 

# bootstrap estimate of the sd
nesda.fix.sd <- apply( nesda.boot.fix.out, 1, sd )

nesda.fix.cv <- nesda.fix.sd/abs(nesda.fix.mean)

# output the data and gzip -f it
write.table( nesda.fix.lbs, "nesda_out_lb_boot_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_lb_boot_fix.txt")

write.table( nesda.fix.ubs, "nesda_out_ub_boot_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_ub_boot_fix.txt")

write.table( nesda.fix.coverage, "nesda_out_C_boot_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t") 
system("gzip -f nesda_out_C_boot_fix.txt")

write.table( nesda.fix.mean, "nesda_out_mn_boot_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_mn_boot_fix.txt")

write.table( nesda.fix.sd, "nesda_out_sd_boot_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_sd_boot_fix.txt")

write.table( nesda.fix.cv, "nesda_out_cv_boot_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_cv_boot_fix.txt")


# Naive T interval bootstrap bounds 
nesda.t.boot.fix.lb <- nesda.fullsamp.coefs + qnorm( .025 ) * nesda.fix.sd # lower bound
nesda.t.boot.fix.ub <- nesda.fullsamp.coefs + qnorm( .975 ) * nesda.fix.sd # upper bound

# coverage rate 
nesda.t.boot.fix.cover <- as.numeric ( sign ( nesda.t.boot.fix.lb * nesda.t.boot.fix.ub  ) <=  0 ) 

# output
write.table( nesda.t.boot.fix.lb, "nesda_out_lb_NTB_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_lb_NTB_fix.txt")

write.table( nesda.t.boot.fix.ub, "nesda_out_ub_NTB_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_ub_NTB_fix.txt")

write.table( nesda.t.boot.fix.cover, "nesda_out_C_NTB_fix.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_C_NTB_fix.txt")

#### COMBINED BOOTSTRAPPING ####
combined.boot.fix.out <- foreach( i = 1:B, .combine = 'cbind', .packages = 'glmnet', .noexport = c( "boot.fix.samples", "boot.fix.X", "boot.fix.Y" ) ) %dopar% {
  
  boot.fix.samples <- sample( nrow( combined.pred.M ), replace = TRUE )
  boot.fix.X <- combined.pred.M[ boot.fix.samples, ]
  boot.fix.Y <- combined.outcome[ boot.fix.samples ]
  
  boot.fix.coefs <- coef( glmnet(x = boot.fix.X, y = boot.fix.Y, lambda = combined.fullsamp.lasso$lambda.min, alpha = 1  )) [ -1 ] 

}
print(date())
# output matrix of bootstrap coefficient estimates and gzip -f it 
write.table( combined.boot.fix.out, "combined_out_boot_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep = "\t" )
system( "gzip -f combined_out_boot_fix.txt" )


# Percentile bootstrap bounds
combined.fix.lbs <- apply( combined.boot.fix.out, 1, quantile, 0.025 ) # lower bound
combined.fix.ubs <- apply( combined.boot.fix.out, 1, quantile, 0.975 ) # upper bound

# coverage rate 
combined.fix.coverage <- as.numeric ( sign ( combined.fix.lbs * combined.fix.ubs  ) <=  0 ) # does the percentile interval contain 0

# bootstrap estimate of the mean
combined.fix.mean <- apply( combined.boot.fix.out, 1, mean ) 

# bootstrap estimate of the sd
combined.fix.sd <- apply( combined.boot.fix.out, 1, sd )

combined.fix.cv <- combined.fix.sd/abs(combined.fix.mean)

# output the data and gzip -f it
write.table( combined.fix.lbs, "combined_out_lb_boot_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_lb_boot_fix.txt")

write.table( combined.fix.ubs, "combined_out_ub_boot_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_ub_boot_fix.txt")

write.table( combined.fix.coverage, "combined_out_C_boot_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t") 
system("gzip -f combined_out_C_boot_fix.txt")

write.table( combined.fix.mean, "combined_out_mn_boot_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_mn_boot_fix.txt")

write.table( combined.fix.sd, "combined_out_sd_boot_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_sd_boot_fix.txt")

write.table( combined.fix.cv, "combined_out_cv_boot_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_cv_boot_fix.txt")


# Naive T interval bootstrap bounds 
combined.t.boot.fix.lb <- combined.fullsamp.coefs + qnorm( .025 ) * combined.fix.sd # lower bound
combined.t.boot.fix.ub <- combined.fullsamp.coefs + qnorm( .975 ) * combined.fix.sd # upper bound

# coverage rate 
combined.t.boot.fix.cover <- as.numeric ( sign ( combined.t.boot.fix.lb * combined.t.boot.fix.ub  ) <=  0 ) 

# output
write.table( combined.t.boot.fix.lb, "combined_out_lb_NTB_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_lb_NTB_fix.txt")

write.table( combined.t.boot.fix.ub, "combined_out_ub_NTB_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_ub_NTB_fix.txt")

write.table( combined.t.boot.fix.cover, "combined_out_C_NTB_fix.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_C_NTB_fix.txt")


#     6. Bootstrap model fit, re-selecting the regularization parameter in each resample
#       6.1. Using 1000 bootstrap replications, estimate the mean, SE, CV, inverted-t CI, and percentile CI of the coefficient for each predictor

#### NTR BOOTSTRAPPING ####
ntr.boot.rand.out <- foreach( i = 1:B, .combine = 'cbind', .packages = 'glmnet', .noexport = c( "boot.rand.samples", "boot.rand.X", "boot.rand.Y", "boot.rand.crossval" ) ) %dopar% {
  
  boot.rand.samples <- sample( nrow( ntr.pred.M ), replace = TRUE )
  boot.rand.X <- ntr.pred.M[ boot.rand.samples, ]
  boot.rand.Y <- ntr.outcome[ boot.rand.samples ]
  
  boot.rand.crossval <- cv.glmnet( x = boot.rand.X, y = boot.rand.Y, alpha = 1 )

  boot.rand.lambda <- boot.rand.crossval$lambda.min

  boot.rand.coefs <- coef( boot.rand.crossval, s = boot.rand.lambda )[ -1 ]
  
}
print(date())
# output matrix of bootstrap coefficient estimates and gzip -f it 
write.table( ntr.boot.rand.out, "ntr_out_boot_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep = "\t" )
system( "gzip -f ntr_out_boot_rand.txt" )


# Percentile bootstrap bounds
ntr.rand.lbs <- apply( ntr.boot.rand.out, 1, quantile, 0.025 ) # lower bound
ntr.rand.ubs <- apply( ntr.boot.rand.out, 1, quantile, 0.975 ) # upper bound

# coverage rate 
ntr.rand.coverage <- as.numeric ( sign ( ntr.rand.lbs * ntr.rand.ubs  ) <=  0 ) # does the percentile interval contain 0

# bootstrap estimate of the mean
ntr.rand.mean <- apply( ntr.boot.rand.out, 1, mean ) 

# bootstrap estimate of the sd
ntr.rand.sd <- apply( ntr.boot.rand.out, 1, sd )

ntr.rand.cv <- ntr.rand.sd/abs(ntr.rand.mean)

# output the data and gzip -f it
write.table( ntr.rand.lbs, "ntr_out_lb_boot_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_lb_boot_rand.txt")

write.table( ntr.rand.ubs, "ntr_out_ub_boot_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_ub_boot_rand.txt")

write.table( ntr.rand.coverage, "ntr_out_C_boot_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t") 
system("gzip -f ntr_out_C_boot_rand.txt")

write.table( ntr.rand.mean, "ntr_out_mn_boot_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_mn_boot_rand.txt")

write.table( ntr.rand.sd, "ntr_out_sd_boot_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_sd_boot_rand.txt")

write.table( ntr.rand.cv, "ntr_out_cv_boot_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_cv_boot_rand.txt")



# Naive T interval bootstrap bounds 
ntr.t.boot.rand.lb <- ntr.fullsamp.coefs + qnorm( .025 ) * ntr.rand.sd # lower bound
ntr.t.boot.rand.ub <- ntr.fullsamp.coefs + qnorm( .975 ) * ntr.rand.sd # upper bound

# coverage rate 
ntr.t.boot.rand.cover <- as.numeric ( sign ( ntr.t.boot.rand.lb * ntr.t.boot.rand.ub  ) <=  0 ) 

# output
write.table( ntr.t.boot.rand.lb, "ntr_out_lb_NTB_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_lb_NTB_rand.txt")

write.table( ntr.t.boot.rand.ub, "ntr_out_ub_NTB_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_ub_NTB_rand.txt")

write.table( ntr.t.boot.rand.cover, "ntr_out_C_NTB_rand.txt", row.names = colnames(ntr.pred.M), quote = FALSE, sep ="\t")
system("gzip -f ntr_out_C_NTB_rand.txt")

#### NESDA BOOTSTRAPPING ####
nesda.boot.rand.out <- foreach( i = 1:B, .combine = 'cbind', .packages = 'glmnet', .noexport = c( "boot.rand.samples", "boot.rand.X", "boot.rand.Y", "boot.rand.crossval" ) ) %dopar% {
  
  boot.rand.samples <- sample( nrow( nesda.pred.M ), replace = TRUE )
  boot.rand.X <- nesda.pred.M[ boot.rand.samples, ]
  boot.rand.Y <- nesda.outcome[ boot.rand.samples ]
  
  boot.rand.crossval <- cv.glmnet( x = boot.rand.X, y = boot.rand.Y, alpha = 1 )

  boot.rand.lambda <- boot.rand.crossval$lambda.min

  boot.rand.coefs <- coef( boot.rand.crossval, s = boot.rand.lambda )[ -1 ]
  
}

print(date())
# output matrix of bootstrap coefficient estimates and gzip -f it 
write.table( nesda.boot.rand.out, "nesda_out_boot_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep = "\t" )
system( "gzip -f nesda_out_boot_rand.txt" )


# Percentile bootstrap bounds
nesda.rand.lbs <- apply( nesda.boot.rand.out, 1, quantile, 0.025 ) # lower bound
nesda.rand.ubs <- apply( nesda.boot.rand.out, 1, quantile, 0.975 ) # upper bound

# coverage rate 
nesda.rand.coverage <- as.numeric ( sign ( nesda.rand.lbs * nesda.rand.ubs  ) <=  0 ) # does the percentile interval contain 0

# bootstrap estimate of the mean
nesda.rand.mean <- apply( nesda.boot.rand.out, 1, mean ) 

# bootstrap estimate of the sd
nesda.rand.sd <- apply( nesda.boot.rand.out, 1, sd )

nesda.rand.cv <- nesda.rand.sd/abs(nesda.rand.mean)

# output the data and gzip -f it
write.table( nesda.rand.lbs, "nesda_out_lb_boot_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_lb_boot_rand.txt")

write.table( nesda.rand.ubs, "nesda_out_ub_boot_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_ub_boot_rand.txt")

write.table( nesda.rand.coverage, "nesda_out_C_boot_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t") 
system("gzip -f nesda_out_C_boot_rand.txt")

write.table( nesda.rand.mean, "nesda_out_mn_boot_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_mn_boot_rand.txt")

write.table( nesda.rand.sd, "nesda_out_sd_boot_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_sd_boot_rand.txt")

write.table( nesda.rand.cv, "nesda_out_cv_boot_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_cv_boot_rand.txt")


# Naive T interval bootstrap bounds 
nesda.t.boot.rand.lb <- nesda.fullsamp.coefs + qnorm( .025 ) * nesda.rand.sd # lower bound
nesda.t.boot.rand.ub <- nesda.fullsamp.coefs + qnorm( .975 ) * nesda.rand.sd # upper bound

# coverage rate 
nesda.t.boot.rand.cover <- as.numeric ( sign ( nesda.t.boot.rand.lb * nesda.t.boot.rand.ub  ) <=  0 ) 

# output
write.table( nesda.t.boot.rand.lb, "nesda_out_lb_NTB_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_lb_NTB_rand.txt")

write.table( nesda.t.boot.rand.ub, "nesda_out_ub_NTB_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_ub_NTB_rand.txt")

write.table( nesda.t.boot.rand.cover, "nesda_out_C_NTB_rand.txt", row.names = colnames(nesda.pred.M), quote = FALSE, sep ="\t")
system("gzip -f nesda_out_C_NTB_rand.txt")

#### COMBINED BOOTSTRAPPING ####

combined.boot.rand.out <- foreach( i = 1:B, .combine = 'cbind', .packages = 'glmnet', .noexport = c( "boot.rand.samples", "boot.rand.X", "boot.rand.Y" ) ) %dopar% {
  
  boot.rand.samples <- sample( nrow( combined.pred.M ), replace = TRUE )
  boot.rand.X <- combined.pred.M[ boot.rand.samples, ]
  boot.rand.Y <- combined.outcome[ boot.rand.samples ]
  
  boot.rand.coefs <- coef( glmnet(x = boot.rand.X, y = boot.rand.Y, lambda = combined.fullsamp.lasso$lambda.min, alpha = 1  )) [ -1 ] 

}
print(date())

# output matrix of bootstrap coefficient estimates and gzip -f it 
write.table( combined.boot.rand.out, "combined_out_boot_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep = "\t" )
system( "gzip -f combined_out_boot_rand.txt" )


# Percentile bootstrap bounds
combined.rand.lbs <- apply( combined.boot.rand.out, 1, quantile, 0.025 ) # lower bound
combined.rand.ubs <- apply( combined.boot.rand.out, 1, quantile, 0.975 ) # upper bound

# coverage rate 
combined.rand.coverage <- as.numeric ( sign ( combined.rand.lbs * combined.rand.ubs  ) <=  0 ) # does the percentile interval contain 0

# bootstrap estimate of the mean
combined.rand.mean <- apply( combined.boot.rand.out, 1, mean ) 

# bootstrap estimate of the sd
combined.rand.sd <- apply( combined.boot.rand.out, 1, sd )

combined.rand.cv <- combined.rand.sd/abs(combined.rand.mean)

# output the data and gzip -f it
write.table( combined.rand.lbs, "combined_out_lb_boot_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_lb_boot_rand.txt")

write.table( combined.rand.ubs, "combined_out_ub_boot_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_ub_boot_rand.txt")

write.table( combined.rand.coverage, "combined_out_C_boot_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t") 
system("gzip -f combined_out_C_boot_rand.txt")

write.table( combined.rand.mean, "combined_out_mn_boot_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_mn_boot_rand.txt")

write.table( combined.rand.sd, "combined_out_sd_boot_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_sd_boot_rand.txt")

write.table( combined.rand.cv, "combined_out_cv_boot_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_cv_boot_rand.txt")


# Naive T interval bootstrap bounds 
combined.t.boot.rand.lb <- combined.fullsamp.coefs + qnorm( .025 ) * combined.rand.sd # lower bound
combined.t.boot.rand.ub <- combined.fullsamp.coefs + qnorm( .975 ) * combined.rand.sd # upper bound

# coverage rate 
combined.t.boot.rand.cover <- as.numeric ( sign ( combined.t.boot.rand.lb * combined.t.boot.rand.ub  ) <=  0 ) 

# output
write.table( combined.t.boot.rand.lb, "combined_out_lb_NTB_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_lb_NTB_rand.txt")

write.table( combined.t.boot.rand.ub, "combined_out_ub_NTB_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_ub_NTB_rand.txt")

write.table( combined.t.boot.rand.cover, "combined_out_C_NTB_rand.txt", row.names = colnames(combined.pred.M), quote = FALSE, sep ="\t")
system("gzip -f combined_out_C_NTB_rand.txt")


#     7. Generate output files
# included predictors, ranked by CV (), coeff.estimate, mean estimates, se estimates, percentile bootstrap CIs 
output.names <- c("FullEst", "RandCV", "FixMean", "FixSD", "RandMean", "RandSD", "FixLB", "FixUB", "RandLB", "RandUB")

# NTR
ntr.included <- which(ntr.fullsamp.coefs != 0)
NTR_includedSTATS <- cbind(ntr.fullsamp.coefs[ntr.included ],  ntr.rand.cv[ntr.included ],  ntr.fix.mean[ntr.included ],  ntr.fix.sd[ntr.included ],  ntr.rand.mean[ntr.included ],  ntr.rand.sd[ntr.included ],  ntr.fix.lbs[ntr.included ],  ntr.fix.ubs[ntr.included ],  ntr.rand.lbs[ntr.included ],  ntr.rand.ubs[ntr.included ])

rownames(NTR_includedSTATS) <- colnames(ntr.pred.M)[ntr.included]
# sort by randm bootstrap CV
NTR_includedSTATS <- NTR_includedSTATS[ order(NTR_includedSTATS[, 2]), ]
colnames(NTR_includedSTATS) <- output.names
write.table(NTR_includedSTATS, "NTR_includedSTATS.txt", quote = FALSE, sep = "\t")

# NESDA
nesda.included <- which(nesda.fullsamp.coefs != 0)
NESDA_includedSTATS <- cbind(nesda.fullsamp.coefs[nesda.included ],  nesda.rand.cv[nesda.included ],  nesda.fix.mean[nesda.included ],  nesda.fix.sd[nesda.included ],  nesda.rand.mean[nesda.included ],  nesda.rand.sd[nesda.included ],  nesda.fix.lbs[nesda.included ],  nesda.fix.ubs[nesda.included ],  nesda.rand.lbs[nesda.included ],  nesda.rand.ubs[nesda.included ])
rownames(NESDA_includedSTATS) <- colnames(ntr.pred.M)[nesda.included]
# sort by randm bootstrap CV
NESDA_includedSTATS <- NESDA_includedSTATS[ order(NESDA_includedSTATS[, 2]), ]
colnames(NESDA_includedSTATS) <- output.names
write.table(NESDA_includedSTATS, "NESDA_includedSTATS.txt", quote = FALSE, sep = "\t")

# COMBINED
combined.included <- which(combined.fullsamp.coefs != 0)
COMBINED_includedSTATS <- cbind(combined.fullsamp.coefs[combined.included ],  combined.rand.cv[combined.included ],  combined.fix.mean[combined.included ],  combined.fix.sd[combined.included ],  combined.rand.mean[combined.included ],  combined.rand.sd[combined.included ],  combined.fix.lbs[combined.included ],  combined.fix.ubs[combined.included ],  combined.rand.lbs[combined.included ],  combined.rand.ubs[combined.included ])
rownames(COMBINED_includedSTATS) <- colnames(ntr.pred.M)[combined.included]
# sort by randm bootstrap CV
COMBINED_includedSTATS <- COMBINED_includedSTATS[ order(COMBINED_includedSTATS[, 2]), ]
colnames(COMBINED_includedSTATS) <- output.names
write.table(COMBINED_includedSTATS, "COMBINED_includedSTATS.txt", quote = FALSE, sep = "\t")



### STATS for excluded predictors

## NTR
ntr.excluded <- which(ntr.fullsamp.coefs == 0)
NTR_excludedSTATS<- cbind(ntr.fullsamp.coefs[ntr.excluded ],  ntr.rand.cv[ntr.excluded ],  ntr.fix.mean[ntr.excluded ],  ntr.fix.sd[ntr.excluded ],  ntr.rand.mean[ntr.excluded ],  ntr.rand.sd[ntr.excluded ],  ntr.fix.lbs[ntr.excluded ],  ntr.fix.ubs[ntr.excluded ],  ntr.rand.lbs[ntr.excluded ],  ntr.rand.ubs[ntr.excluded ])
rownames(NTR_excludedSTATS) <- colnames(ntr.pred.M)[ntr.excluded]
# reorder excluded variables by random SD
NTR_excludedSTATS <- NTR_excludedSTATS[ order(NTR_excludedSTATS[, 6]), ]
colnames(NTR_excludedSTATS) <- output.names
write.table(NTR_excludedSTATS, "NTR_excludedSTATS.txt", quote = FALSE, sep = "\t")

## NESDA
nesda.excluded <- which(nesda.fullsamp.coefs == 0)
NESDA_excludedSTATS<- cbind(nesda.fullsamp.coefs[nesda.excluded ],  nesda.rand.cv[nesda.excluded ],  nesda.fix.mean[nesda.excluded ],  nesda.fix.sd[nesda.excluded ],  nesda.rand.mean[nesda.excluded ],  nesda.rand.sd[nesda.excluded ],  nesda.fix.lbs[nesda.excluded ],  nesda.fix.ubs[nesda.excluded ],  nesda.rand.lbs[nesda.excluded ],  nesda.rand.ubs[nesda.excluded ])
rownames(NESDA_excludedSTATS) <- colnames(ntr.pred.M)[nesda.excluded]

# reorder excluded variables by random SD
NESDA_excludedSTATS <- NESDA_excludedSTATS[ order(NESDA_excludedSTATS[, 6]), ]
colnames(NESDA_excludedSTATS) <- output.names
write.table(NESDA_excludedSTATS, "NESDA_excludedSTATS.txt", quote = FALSE, sep = "\t")


## COMBINED
combined.excluded <- which(combined.fullsamp.coefs == 0)
COMBINED_excludedSTATS<- cbind(combined.fullsamp.coefs[combined.excluded ],  combined.rand.cv[combined.excluded ],  combined.fix.mean[combined.excluded ],  combined.fix.sd[combined.excluded ],  combined.rand.mean[combined.excluded ],  combined.rand.sd[combined.excluded ],  combined.fix.lbs[combined.excluded ],  combined.fix.ubs[combined.excluded ],  combined.rand.lbs[combined.excluded ],  combined.rand.ubs[combined.excluded ])
rownames(COMBINED_excludedSTATS) <- colnames(ntr.pred.M)[combined.excluded]

# reorder excluded variables by random SD
COMBINED_excludedSTATS <- COMBINED_excludedSTATS[ order(COMBINED_excludedSTATS[, 6]), ]
colnames(COMBINED_excludedSTATS) <- output.names
write.table(COMBINED_excludedSTATS, "COMBINED_excludedSTATS.txt", quote = FALSE, sep = "\t")

# graphical output

# generate appropriate data frame for the plot
print(date())
# combine results for fixed lambda
fixed.results <- rbind( cbind(combined.fix.mean, combined.fix.sd),
                       cbind(ntr.fix.mean, ntr.fix.sd),
                       cbind(nesda.fix.mean, nesda.fix.sd))

# combine results for random lambda
random.results <- rbind( cbind(combined.rand.mean, combined.rand.sd),
                       cbind(ntr.rand.mean, ntr.rand.sd),
                       cbind(nesda.rand.mean, nesda.rand.sd))


full.results <- rbind(fixed.results, random.results)

nobs <- length(combined.fix.mean)

facet.factors <- cbind( c( rep("fixed", ( nobs * 3 )), rep("random", ( nobs * 3 )) ),
                       rep( c( rep("combined", nobs), rep("ntr", nobs), rep("nesda", nobs) ), 2)
                       )

bga_interaction_analysis_bootstrap_moments <- cbind(full.results, facet.factors)

colnames(bga_interaction_analysis_bootstrap_moments) <- c("boot.mn", "boot.sd", "lambda", "study")

bga_interaction_analysis_bootstrap_moments <- as.data.frame(bga_interaction_analysis_bootstrap_moments)

head(bga_interaction_analysis_bootstrap_moments)

bga_interaction_analysis_bootstrap_moments$boot.mn <- as.numeric(bga_interaction_analysis_bootstrap_moments$boot.mn)
bga_interaction_analysis_bootstrap_moments$boot.sd <- as.numeric(bga_interaction_analysis_bootstrap_moments$boot.sd)

qplot( boot.mn, boot.sd, data = bga_interaction_analysis_bootstrap_moments, geom = "point" ) + facet_grid(lambda ~ study, scales = "fixed")
print(date())
ggsave(file = "bootstrap_moments.png")

print(date())

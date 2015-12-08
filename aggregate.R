# aggregation script for dissertation simulations

# Joins all tables resulting from the run script, produces the 108 x 13 total table of results
# And analyses on each margin

# 15 sets of results, each with 10800 rows
# merge together into a giant data frame
# calculate summary and marginal stats
# present summary and marginal stats in nice .tex form  

Ns <- c(6, 18, 54, 100)
ps <- c(5, 10, 20)
rs <- c(0.2, 0.5)
wts <- 1:2
model <- 1:3
method <- 1:5
cells_table <- expand.grid(Ns, ps, rs, wts, model)
analysis_table <- expand.grid(Ns, ps, rs, wts, model, method)

  resultsmatrix <- matrix(ncol = 20, nrow = 43260)
    tempnamelist <- system("ls -1 analyses_matrix*txt", intern = TRUE)
    startcounter <- 1
    endcounter <- 0
    for (name in tempnamelist){
      tempresults <- data.matrix(read.table(name, as.is = TRUE, comment.char = "", header = TRUE))
      endcounter <- endcounter + nrow(tempresults)
      resultsmatrix[startcounter:endcounter,] <- unlist(tempresults)
      startcounter <- (1 + endcounter)

    }
     

  colnames(resultsmatrix) <- c("N", "p", "R2", "wts", "model", "meth", "achR", "vsp", "tpr", "fpr", "hvsp", "htpr", "hfpr", "sachR", "svsp", "stpr", "sfpr", "shvsp", "shtpr", "shfpr")

  analyses_matrix <- matrix(nrow = nrow(analysis_table), ncol = 20)
  for (q in 1:nrow(analysis_table)){
    W1 <- which(resultsmatrix[,'N'] == analysis_table[q,1])
    W2 <- which(resultsmatrix[,'p'] == analysis_table[q,2])
    W3 <- which(resultsmatrix[,'R2'] == analysis_table[q,3])
    W4 <- which(resultsmatrix[,'wts'] == as.numeric(analysis_table[q,4]))
    W5 <- which(resultsmatrix[,'model'] == as.numeric(analysis_table[q,5]))
    W6 <- which(resultsmatrix[,'meth'] == analysis_table[q,6])

    indices <- intersect(W1, intersect( W2, intersect( W3, intersect(W4, intersect(W5, W6)))))
    meanvec <- apply(resultsmatrix[indices, 7:20],2, mean, na.rm = TRUE)

    analyses_matrix[q, 1:6] <- unlist(analysis_table[q,])
    analyses_matrix[q, 7:20] <- unlist(meanvec)

  }
  colnames(analyses_matrix) <- c("N", "p", "R2", "wts", "model","meth", "achR","vsp", "tpr", "fpr", "hvsp", "htpr", "hfpr","sachR","svsp", "stpr", "sfpr", "shvsp", "shtpr", "shfpr")
  write.table(signif(analyses_matrix,4),"hi_analyses_matrix.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

  



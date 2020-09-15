extract_lmekin <- function(x, ...) {
  fcoef <- x$coefficients$fixed
  if (length(fcoef) > 0)  { # Not a ~1 model
    se <- sqrt(diag(x$var))
    tmp <- cbind(fcoef, se, round(fcoef/se,2),
                 signif(1 - pchisq((fcoef/ se)^2, 1), 2))
    dimnames(tmp) <- list(names(fcoef), c("Value", "Std Error", "z", "p")) 
  }
  results = data.frame(n = x$n, tmp)
}

results_miRNA_lmekin <- function(miRNA,data, kmat, cov){
  # miRNA - vector of miRNA names
  # data - all data, miRNA, top POS
  # kmat - kinship matrix
  # cov - vector of covariates
  results <- NULL
  for (j in 1:length(miRNA)){
    print(miRNA[j])
    results1 <- NULL
    for (i in 1:length(data$POS)){
      print(i)
      x <- Sys.time()
      lmmkin1 <- lmekin(formula = as.formula(paste(miRNA[j],"~gt_DS*BB+", cov, "+ (1|colnames(kmat))")), data = data$data[[i]],
                        varlist = list(kmat))
      new_results <- data.frame(extract_lmekin(lmmkin1), POS = data$POS[[i]], Gene = data$Gene[[i]], miRNA = miRNA[j])
      new_results$var <- row.names(new_results)
      results1 <- rbind(results1, new_results)
      y<-Sys.time()
      print(y-x)
    }
    results <- rbind(results, results1)
  }
  return(results)
}
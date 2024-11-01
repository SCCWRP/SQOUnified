SoCalRivpacs <- function(Pcutoff = 0.5,
                         reference.groups = socal.reference.groups,
                         observed.predictors = socal.example.habitat,
                         reference.taxa = socal.reference.taxa,
                         group.means = socal.reference.group.means,
                         reference.cov = socal.reference.covariance,
                         observed.taxa = socal.example.taxa) {  
  
  # Pcutoff is the probability cutoff  
  
  # Names of predictor variables.
  predictor.variables <- c("Latitude", "Longitude", "SampleDepth")
  
  # ----- Format Observed Data -----
  
  FormatObservedData <- function() {
    
    # Align observed (user) data columns with reference data columns. Columns in same
    # order. Observed data may have a different number of taxa (columns) than
    # reference data. 
    
    # Convert observed.taxa to presence/absence (0/1)
    tmp.pa <- observed.taxa
    tmp.pa[tmp.pa > 0] <- 1
    
    # Align rows using predictor variables.
    tmp.pa <- tmp.pa[row.names(observed.predictors), ]    # !!! is this required???
    
    # Container matrix.
    n.observed.sites <- dim(tmp.pa)[[1]]
    n.reference.taxa <- dim(reference.taxa)[[2]]
    
    observed.taxa.pa <- matrix(rep(0, times = n.observed.sites * n.reference.taxa), 
                               nrow = n.observed.sites, 
                               dimnames = list(rownames(tmp.pa), names(reference.taxa)))
    
    # Fill container with observed data.
    col.match <- match(dimnames(observed.taxa.pa)[[2]], dimnames(tmp.pa)[[2]])
    
    for(i in 1:n.reference.taxa) {
      
      if(!is.na(col.match[i])) observed.taxa.pa[, i] <- tmp.pa[, col.match[i]]
      
    }
    
    # The matrix observed.taxa.pa contains the observed.scores used for O/E.
    return(observed.taxa.pa)
    
  }
  
  observed.data <- FormatObservedData()  
  
  # ----- Calculate Expected Data -----
  
  CalculateExpectedData <- function() {
    
    # Calculate probability of sites belonging to groups. Follow RIVPACS assumption
    # of weighting the group probabilities by reference group size. Flags outlier
    # sites, using the chi-squared statistic.
    
    # Definitions.
    n.predictor.variables <- length(predictor.variables)
    group.size <- table(reference.groups)
    n.groups <- length(group.size)
    
    # Chi-squared values for flagging outlier samples.
    degrees.freedom <- min(c(n.predictor.variables, (n.groups - 1)))
    crit.01 <- qchisq(0.99, df = degrees.freedom)
    crit.05 <- qchisq(0.95, df = degrees.freedom)
    
    # Container for probabilities.
    n.observed.sites.filtered <- dim(observed.predictors)[[1]]
    
    group.probabilities <- matrix(rep(0, n.observed.sites.filtered * n.groups), 
                                  nrow = n.observed.sites.filtered,
                                  dimnames = list(dimnames(observed.predictors)[[1]], 
                                                  dimnames(group.means)[[1]]))
    
    # Container for outlier flags and minimum distance.
    outlier.flag <- data.frame(outlier.05 = rep(0, n.observed.sites.filtered), 
                               outlier.01 = rep(0, n.observed.sites.filtered), 
                               min.distance = rep(0, n.observed.sites.filtered), 
                               row.names = dimnames(observed.predictors)[[1]])
    
    # calculate group membership probabilities for each sample and find outliers.
    for(i in 1:n.observed.sites.filtered) {
      
      # Squared Mahalanobis distance from current sample to each group mean.
      distances <- mahalanobis(group.means, 
                               observed.predictors[i,], 
                               reference.cov, 
                               inverted = TRUE)  
      
      group.probabilities[i,] <- group.size * exp(-0.5 * distances) # see Clarke et al. (2000) 
      group.probabilities[i,] <- group.probabilities[i, ] / sum(group.probabilities[i, ])
            
      # Outlier criteria is minimum distance.
      outlier.flag$min.distance[i] <- min(distances) 
      
      # Check for outliers. Each sample is either a pass (0) or fail (1).
      if(outlier.flag$min.distance[i] > crit.05) outlier.flag[i, "outlier.05"] <- 1
      if(outlier.flag$min.distance[i] > crit.01) outlier.flag[i, "outlier.01"] <- 1
      
    } 
    
    # Occurrence frequencies of all taxa in the reference groups.
    freq.in.group <- apply(reference.taxa, 2, 
                           function(x){tapply(x, reference.groups, function(y){sum(y) / length(y)})})
    
    # Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).
    predicted.prob.all <- group.probabilities %*% freq.in.group
    
    # predicted.prob.all are the predicted (expected) probabilites.   
    expected.data <- list(predicted = predicted.prob.all, outliers = outlier.flag, n = n.observed.sites.filtered)
    
  }
  
  expected.data <- CalculateExpectedData()
  
  # ----- Calculate Scores ----- 
  
  CalculateScores <- function() {
    
    observed.score <- vector(mode = "numeric", length = expected.data$n)
    expected.score <- vector(mode = "numeric", length = expected.data$n)
    BC <- vector(mode = "numeric", length = expected.data$n) # Bray-Curtis dissimilarity 
    
    
    
    for(i in 1:expected.data$n) {
      
      predicted.prob <- expected.data$predicted[i, ] # predicted probabilities for current sample
           
      taxa.subset <- names(predicted.prob)[predicted.prob >= Pcutoff]  # subset of taxa with probabilities >= Pcutoff
         
      expected.prob <- predicted.prob[taxa.subset] # probabilites for subset of included taxa
      
      observed.pa <- observed.data[i, taxa.subset] # observed presence/absence for those taxa
            
      observed.score[i] <- sum(observed.pa) # observed richness (O)
      expected.score[i] <- sum(expected.prob) # expected richness (E)
      BC[i] <- sum(abs(observed.pa - expected.prob)) / 
        (observed.score[i] + expected.score[i]) # BC value
      
    }
    
    O.over.E <- observed.score/expected.score 
    
    stats <- data.frame(stations = row.names(observed.predictors), 
                        O = observed.score, 
                        E = round(expected.score, digits = 4), 
                        O.over.E = round(O.over.E, digits = 4))
    
#     stats <- data.frame(stations = row.names(observed.predictors), 
#                         O = observed.score, 
#                         E = expected.score, 
#                         O.over.E = O.over.E)
    
    stats$outlier.05 <- expected.data$outliers$outlier.05
    stats$outlier.01 <- expected.data$outliers$outlier.01
        
    # Convert to "PASS" or "FAIL"
    stats$outlier.05[stats$outlier.05 == 0] <- "PASS"
    stats$outlier.05[stats$outlier.05 == 1] <- "FAIL"
    
    stats$outlier.01[stats$outlier.01 == 0] <- "PASS"
    stats$outlier.01[stats$outlier.01 == 1] <- "FAIL"
    
    #   mean.O.over.E <- mean(OE.stats$O.over.E)
    #   stdev.O.over.E <- sqrt(var(OE.stats$O.over.E))
    
    return(stats)
    
  }

  results <- list(oe.table = CalculateScores(), 
                  observed = observed.data, 
                  predicted = expected.data$predicted,
                  Pcutoff = Pcutoff,
                  region = "scb")
  
}
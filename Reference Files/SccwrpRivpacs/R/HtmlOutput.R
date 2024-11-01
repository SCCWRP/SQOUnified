HtmlOutput <- function(rivpacs,
                       timestamp,
                       user.filename,
                       path = "/var/www/sqo/files/") {
  
  # ----- Rivpacs Output HTML table for results -----
  
  RivpacsOutput <- function(x = rivpacs$oe.table, 
                            Pcutoff = rivpacs$Pcutoff, 
                            region = rivpacs$region,
                            filename = paste(path, timestamp, ".ro", rivpacs$region, ".html", sep = "")) {
      
    # Without setting column width
    # cat("<!DOCTYPE html>\n<html>\n<head>\n<style type=\"text/css\">\ntable {border-collapse:collapse;}\ntable, td, th {border: 1px solid black;}\nbody {font-family: monospace;}\n</style>\n</head>\n<body>\n",
    #     file = filename)
    
    cat("<!DOCTYPE html>\n<html>\n<head>\n<style type=\"text/css\">\ntable {border-collapse:collapse; width: 45%;}\ntable, td, th {border: 1px solid black;}\nbody {font-family: monospace;}\n</style>\n</head>\n<body>\n",
        file = filename)
    
    cat("<p> RIVPACS Output", format(Sys.time(), "%a %b %d %X %Y"), "</p>\n", file = filename, append = TRUE)
    
    cat("<p> SccwrpRivpacs Version", as.character(packageVersion("SccwrpRivpacs")), "</p>\n", file = filename, append = TRUE)
    
    cat("<h2> O/E Table </h2>\n", file = filename, append = TRUE)
    
    if(region == "scb") cat("<h3>Southern California Bays and Estuaries</h3>\n", file = filename, append = TRUE)
    if(region == "sfb") cat("<h3>San Francisco Bay Polyhaline</h3>\n", file = filename, append = TRUE)
    
    cat("<h3> Input File:", user.filename, "</h3>\n", file = filename, append = TRUE)
    
    cat("<h3> Probability Cutoff =", Pcutoff, "</h3>\n", file = filename, append = TRUE)
    
    output <- xtable(x,
                     digits = 4,
                     display = c("d", "s", "d", "f", "f", "d", "d"), 
                     align = c("l", "l", "c", "c", "c", "c", "c"))
    
    # Without modifuing column names
    # print(output, type = "html", file = filename, append = TRUE, include.rownames = FALSE))
    
    captured <- capture.output(print(output, type = "html", include.rownames = FALSE))
    
    captured[4] <- gsub("stations", "Stations", x = captured[4], fixed = TRUE)
    captured[4] <- gsub("O.over.E", "O/E", x = captured[4], fixed = TRUE)
    captured[4] <- gsub("outlier.05", "Data Meets Model Assumptions (P = 0.95)", x = captured[4], fixed = TRUE)
    captured[4] <- gsub("outlier.01", "Data Meets Model Assumptions (P = 0.99)", x = captured[4], fixed = TRUE)
    
    cat(captured, file = filename, sep = "\n", append = TRUE)
    
    cat("</body>\n</html>", file = filename, append = TRUE)
    
  }
  
  # ----- Probability Comparison HTML table for species detection -----
  
  SpeciesDetection <- function(x = rivpacs$observed, 
                               y = rivpacs$predicted, 
                               region = rivpacs$region,
                               filename = paste(path, timestamp, ".pc", rivpacs$region, ".html", sep = "")) {
    
    # Convert observed.data to tall format
    
    obs <- as.data.frame(x)
    obs.long <- reshape(data = obs, varying = list(names(obs)), v.names = "Observed", timevar = "Taxa", 
                        idvar = "Station", ids = row.names(obs), times = names(obs), direction = "long")
    
    # Convert predicted data to tall format
    
    pred <- as.data.frame(y)
    pred.long <- reshape(data = pred, varying = list(names(pred)), v.names = "P", timevar = "Taxa", 
                         idvar = "Station", ids = row.names(pred), times = names(pred), direction = "long")
    
    # Compare based on rules
    
    comparison <- merge(x = obs.long, y = pred.long, by = c("Taxa", "Station"), all = TRUE)
    
    comparison$Taxa <- sub("__", " (", comparison$Taxa, fixed = TRUE)
    comparison$Taxa <- sub("__", ") ", comparison$Taxa, fixed = TRUE)
    comparison$Taxa <- gsub("_", " ", comparison$Taxa, fixed = TRUE)
    
    pred.not.obs <- comparison[comparison$Observed == 0 & comparison$P >= 0.5, ]
    
    obs.not.pred <- comparison[comparison$Observed == 1 & comparison$P < 0.5, ]
    
    pred.and.obs <- comparison[comparison$Observed == 1 & comparison$P >= 0.5, ]
    
    # HTML tables
    
    cat("<!DOCTYPE html>\n<html>\n<head>\n<style type=\"text/css\">\ntable {border-collapse:collapse; width: 45%;}\ntable, td, th {border: 1px solid black;}\nbody {font-family: monospace;}\n</style>\n</head>\n<body>\n",
        file = filename)
    
    cat("<p> RIVPACS Output", format(Sys.time(), "%a %b %d %X %Y"), "</p>\n", file = filename, append = TRUE)

    cat("<p> SccwrpRivpacs Version", as.character(packageVersion("SccwrpRivpacs")), "</p>\n", file = filename, append = TRUE)
    
    cat("<h2> Probability Comparison </h2>\n", file = filename, append = TRUE)
    
    if(region == "scb") cat("<h3>Southern California Bays and Estuaries</h3>\n", file = filename, append = TRUE)
    if(region == "sfb") cat("<h3>San Francisco Bay Polyhaline</h3>\n", file = filename, append = TRUE)
    
    cat("<h3> Input File:", user.filename, "</h3>\n", file = filename, append = TRUE)
    
    cat("<h3> Three tables: </h3>\n", file = filename, append = TRUE)
    
    # Predicted not observed
    
    cat("<h3> 1. Taxa predicted (P >= 0.5) but not observed </h3>\n", file = filename, append = TRUE)
    
    if(nrow(pred.not.obs) == 0) {
      
      cat("None",file = filename, append = TRUE)
      
      cat("<br>\n", file = filename, append = TRUE)
      
    } else {
      
    pred.not.obs$P <- formatC(round(pred.not.obs$P, digits = 4), format = "f", digits = 4)
    
    pred.not.obs$Observed[pred.not.obs$Observed == 0] <- "FALSE"
    
    comparison1 <- xtable(pred.not.obs, digits = 4, align = c("l", "l", "c", "c", "c"))
    
    print(comparison1, type = "html", file = filename, append = TRUE, include.rownames = FALSE)
    
    cat("<br>\n", file = filename, append = TRUE)
    
    }
    
    # Not predicted but observed
    
    cat("<h3> 2. Taxa not predicted (P < 0.5) but observed </h3>\n", file = filename, append = TRUE)
    
    if(nrow(obs.not.pred) == 0) {
      
      cat("None",file = filename, append = TRUE)
      
      cat("<br>\n", file = filename, append = TRUE)
      
    } else {
    
    obs.not.pred$P <- formatC(round(obs.not.pred$P, digits = 4), format = "f", digits = 4)
    
    obs.not.pred$Observed[obs.not.pred$Observed == 1] <- "TRUE"
    
    comparison2 <- xtable(obs.not.pred, digits = 4, align = c("l", "l", "c", "c", "c"))
    
    print(comparison2, type = "html", file = filename, append = TRUE, include.rownames = FALSE)
    
    cat("<br>\n", file = filename, append = TRUE)
    
    }
    
    # Predicted and observed
    
    cat("<h3> 3. Taxa predicted (P >= 0.5) and observed </h3>\n", file = filename, append = TRUE)
    
    if(nrow(pred.and.obs) == 0) {
      
      cat("None",file = filename, append = TRUE)
      
      cat("<br>\n", file = filename, append = TRUE)
      
    } else {
    
    pred.and.obs$P <- formatC(round(pred.and.obs$P, digits = 4), format = "f", digits = 4)
    
    pred.and.obs$Observed[pred.and.obs$Observed == 1] <- "TRUE"
    
    comparison3 <- xtable(pred.and.obs, digits = 4, align = c("l", "l", "c", "c", "c"))
    
    print(comparison3, type = "html", file = filename, append = TRUE, include.rownames = FALSE)
    
    cat("</body>\n</html>", file = filename, append = TRUE)
    
    }
    
  }
  
  
  # ----- Probability Matrix HTML table for probability matrix -----
  
  ProbabilityMatrix <- function(x = rivpacs$predicted, 
                                region = rivpacs$region,
                                filename = paste(path, timestamp, ".pm", rivpacs$region, ".html", sep = "")) {
    
    attr(x, "dimnames")[[2]] <- sub("__", " (", attr(x, "dimnames")[[2]], fixed = TRUE)
    attr(x, "dimnames")[[2]] <- sub("__", ") ", attr(x, "dimnames")[[2]], fixed = TRUE)
    attr(x, "dimnames")[[2]] <- gsub("_", " ", attr(x, "dimnames")[[2]], fixed = TRUE)
    
    matrix.output <- formatC(round(x, digits = 4), format = "f", digits = 4)
    
    Station <- row.names(matrix.output)
    
    matrix.output <- cbind(Station, matrix.output)
    
    cat("<!DOCTYPE html>\n<html>\n<head>\n<style type=\"text/css\">\ntable {border-collapse:collapse;}\ntable, td, th {border: 1px solid black;}\nbody {font-family: monospace;}\n</style>\n</head>\n<body>\n",
        file = filename)
    
    cat("<p> RIVPACS Output", format(Sys.time(), "%a %b %d %X %Y"), "</p>\n", file = filename, append = TRUE)
    
    cat("<p> SccwrpRivpacs Version", as.character(packageVersion("SccwrpRivpacs")), "</p>\n", file = filename, append = TRUE)
    
    cat("<h2> Probability Matrix </h2>\n", file = filename, append = TRUE)
    
    if(region == "scb") cat("<h3>Southern California Bays and Estuaries</h3>\n", file = filename, append = TRUE)
    if(region == "sfb") cat("<h3>San Francisco Bay Polyhaline</h3>\n", file = filename, append = TRUE)
    
    cat("<h3> Input File:", user.filename, "</h3>\n", file = filename, append = TRUE)
    
    probability.matrix <- xtable(matrix.output, digits = 6, align = c("l", "l", rep("c", times = ncol(matrix.output) - 1)))
    
    print(probability.matrix, type = "html", file = filename, append = TRUE, include.rownames = FALSE)
    
    # Not sure this approach for coloring cells will work. The rules are not just > or < a value.
    # captured <- capture.output(print(probability.matrix, type = "html")) #, file = filename, append = TRUE))
    # 
    # captured.color <- vector(mode = "character", length = length(captured))
    # for(i in 1:length(captured)) {
    #  captured.color[i] <- gsub("> 0.210670", "bgcolor=\"red\"> 0.478148", x = captured[i], fixed = TRUE)
    # }
    # cat(captured.color, file = filename, sep = "\n", append = TRUE)
    
    cat("</body>\n</html>", file = filename, append = TRUE)
    
  }
  
  RivpacsOutput()
  SpeciesDetection()
  ProbabilityMatrix()
   
}

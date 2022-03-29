# This function is private, only used by MAMBI
EQR <- function(data) {
  segm <- data[nrow(data),] - data[(nrow(data)-1),]
  vett <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
  for (k in 1: ncol(data)) {vett[, k] <- data[(nrow(data)-1), k]}
  vett <- data - vett
  ris <- round((vett %*% segm / sqrt(sum(segm*segm))) / sqrt(sum(segm*segm)), 3)
  return(ris)
}

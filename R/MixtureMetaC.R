#' Mixture population model approach for meta analysis of continuous outcome
#'
#' @param DSN Input data frame
#' @param WEIGHT Weights used to combine studies: '"size"', '"equal"' or a numeric vector
#' @param ALPHA Two-sided significance level
#' @param SE 'TRUE' for standard error included in input data frame instead of standard deviation and 'FALSE' (the default) otherwise
#' @param DEMOG A character vector for baseline characteristics included in input data frame
#' @returns A data frame
#' @examples
#' MixtureMetaC(DSN = ELDL, WEIGHT = "size", DEMOG = c("Age", "Male", "LDL"))
#' MixtureMetaC(DSN = ELDL, WEIGHT = "equal")
#' MixtureMetaC(DSN = ELDL, WEIGHT = rep(1, 16))
#' @export
MixtureMetaC <- function(DSN = NULL, WEIGHT = NULL, ALPHA = 0.05, SE = FALSE, DEMOG = NULL){
  if (tolower(WEIGHT)[1] == "equal") {
    W <- rep(1 / nrow(DSN), nrow(DSN))
  } else {
    if (tolower(WEIGHT)[1] == "size") {
      W <- with(DSN, (N0 + N1) / sum(N0 + N1))
    } else {
      if (length(WEIGHT) == nrow(DSN) & is.numeric(WEIGHT)) {
        W <- WEIGHT / sum(WEIGHT)
      } else {
        stop(paste("The mixture weight needs to be 'equal', 'size' or a numeric vector with the same length as the number of studies."))
      }
    }
  }

  if (SE == FALSE) {
    V0 <- with(DSN, SD0^2 / N0)
    V1 <- with(DSN, SD1^2 / N1)
  } else {
    V0 <- with(DSN, SD0^2)
    V1 <- with(DSN, SD1^2)
  }

  IRES <- with(DSN, matrix(c(E1 - E0,
                             E1 - E0 - qnorm(1 - ALPHA / 2) * sqrt(V0 + V1),
                             E1 - E0 + qnorm(1 - ALPHA / 2) * sqrt(V0 + V1)), ncol = 3))
  IRES <- cbind(IRES, W)
  colnames(IRES) <- c("DIFF", paste0(round((1 - ALPHA) * 100, 0), "% CI_lower"), paste0(round((1 - ALPHA) * 100, 0), "% CI_lower"), "Weights")

  RESULT <- cbind(DSN, IRES)

  ME0 <- sum(W * DSN$E0)
  ME1 <- sum(W * DSN$E1)

  SE0 <- sqrt(sum(W^2 * V0))
  SE1 <- sqrt(sum(W^2 * V1))

  MD  <- ME1 - ME0
  MDl <- MD - qnorm(1 - ALPHA / 2) * sqrt(SE0^2 + SE1^2)
  MDu <- MD + qnorm(1 - ALPHA / 2) * sqrt(SE0^2 + SE1^2)

  RESULT <- rbind(RESULT[, -1], rep(NA, ncol(DSN) + 3))

  RESULT$E1[nrow(RESULT)]  <- ME1
  RESULT$SD1[nrow(RESULT)] <- SE1
  RESULT$E0[nrow(RESULT)]  <- ME0
  RESULT$SD0[nrow(RESULT)] <- SE0
  RESULT[nrow(RESULT), c((ncol(RESULT) - 3) : (ncol(RESULT) - 1))] <- c(MD, MDl, MDu)

  if (length(DEMOG) > 0){
    DEMOV <- paste("RESULT", DEMOG, sep = "$")
    RESV  <- paste0(paste("RESULT", DEMOG, sep = "$"), "[nrow(RESULT)] <- sum(", DEMOV, " * c(W, NA), na.rm = TRUE)")

    for (i in 1 : length(DEMOG)){
      eval(parse(text = RESV[i]))
    }
  }

  NAMES <- colnames(RESULT)
  RESULT$STUDY <- c(DSN[, 1], paste("Mixture weights =", paste(WEIGHT, collapse = ",")))
  RESULT <- RESULT[, c("STUDY", NAMES)]

  return(RESULT)
}

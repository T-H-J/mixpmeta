#' Mixture population model approach for meta analysis of binary outcome
#'
#' @param DSN Input data frame
#' @param WEIGHT Weights used to combine studies: '"size"', '"equal"' or a numeric vector
#' @param ALPHA Two-sided significance level
#' @param NBOOT Number of iterations. When default (NBOOT = 2000), 2,000 sets of realizations are generated for calculating confidence intervals of risk ratio estimated by the mixture population model approach.
#' @param DEMOG A character vector for baseline characteristics included in input data frame
#' @returns A data frame
#' @examples
#' MixtureMetaB(DSN = EHFE, WEIGHT = "size", DEMOG = c("AGE", "MALE", "BMI", "DD", "FU"))
#' MixtureMetaB(DSN = EHFE, WEIGHT = "equal")
#' MixtureMetaB(DSN = EHFE, WEIGHT = rep(1, 6))
#' @export
MixtureMetaB <- function(DSN = NULL, WEIGHT = NULL, ALPHA = 0.05, NBOOT = 2000, DEMOG = NULL){
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

  PR1   <- with(DSN, R1 / N1)
  PR0   <- with(DSN, R0 / N0)

  eRR   <- PR1 / PR0
  vleRR <- with(DSN, (N1 - R1) / R1 / N1 + (N0 - R0) / R0 / N0)
  eRRl  <- exp(log(eRR) - qnorm(1 - ALPHA / 2) * sqrt(vleRR))
  eRRu  <- exp(log(eRR) + qnorm(1 - ALPHA / 2) * sqrt(vleRR))

  eRD   <- PR1 - PR0
  veRD  <- with(DSN, R1 * (N1 - R1) / N1^3 + R0 * (N0 - R0) / N0^3)
  eRDl  <- eRD - qnorm(1 - ALPHA / 2) * sqrt(veRD)
  eRDu  <- eRD + qnorm(1 - ALPHA / 2) * sqrt(veRD)

  IRES <- with(DSN, matrix(c(PR1, PR0, eRR, eRRl, eRRu, eRD, eRDl, eRDu), ncol = 8))
  IRES <- cbind(IRES, W)
  colnames(IRES) <- c("PR1", "PR0",
                      "Risk ratio", paste0(round((1 - ALPHA) * 100, 0), "% CI_RR_lower"), paste0(round((1 - ALPHA) * 100, 0), "% CI_RR_upper"),
                      "Risk difference", paste0(round((1 - ALPHA) * 100, 0), "% CI_RD_lower"), paste0(round((1 - ALPHA) * 100, 0), "% CI_RD_upper"),
                      "Weights")

  RESULT <- cbind(DSN, IRES)

  MP0 <- sum(PR0 * W)
  MP1 <- sum(PR1 * W)

  MRR <- MP1 / MP0

  R1B <- matrix(rbinom(length(DSN$R1) * NBOOT, size = DSN$N1, prob = PR1), nrow = length(DSN$R1))
  R0B <- matrix(rbinom(length(DSN$R0) * NBOOT, size = DSN$N0, prob = PR0), nrow = length(DSN$R0))

  mgRRB <- numeric(NBOOT)

  for (i in 1 : NBOOT){
    mgRRB[i] <- sum(R1B[, i] / DSN$N1 * W) / sum(R0B[, i] / DSN$N0 * W)
  }

  MRRlu <- quantile(mgRRB, probs = c(ALPHA / 2, 1 - ALPHA / 2))

  MRD <- MP1 - MP0
  vMRD <- sum(W^2 * with(DSN, R1 * (N1 - R1) / N1^3)) + sum(W^2 * with(DSN, R0 * (N0 - R0) / N0^3))
  MRDl <- MRD - qnorm(1 - ALPHA / 2) * sqrt(vMRD)
  MRDu <- MRD + qnorm(1 - ALPHA / 2) * sqrt(vMRD)

  RESULT <- rbind(RESULT[, -1], rep(NA, ncol(DSN) + 8))

  RESULT$PR1[nrow(RESULT)]  <- MP1
  RESULT$PR0[nrow(RESULT)]  <- MP0
  RESULT[nrow(RESULT), c((ncol(RESULT) - 6) : (ncol(RESULT) - 1))] <- c(MRR, MRRlu, MRD, MRDl, MRDu)

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

  RESULT
}


#' Mixture population model approach for meta analysis of binary outcome
#'
#' Mixture population model approach to combine information across multiple studies for binary outcomes. This simple and robust procedure based on a mixture population concept can provide a clinically meaningful between-group contrast summary for a well-defined study population.
#' @param dsn Input data frame with study names in the first column, where variable N0 and N1 are the number of subjects (study size) in control group and treatment group, respectively, and variable R0 and R1 are the number of events in control group and treatment group, respectively.
#' @param weight Weights used to combine studies: '"size"', '"equal"' or a numeric vector
#' @param alpha Two-sided significance level
#' @param nboot Number of iterations. When default (nboot = 2000), 2,000 sets of realizations are generated for calculating confidence intervals of risk ratio estimated by the mixture population model approach.
#' @param demog A character vector for baseline characteristics included in the input data frame
#' @returns A data frame
#' @examples
#' MixtureMetaB(dsn = EHFE, weight = "size", demog = c("AGE", "MALE", "BMI", "DD", "FU"))
#' MixtureMetaB(dsn = EHFE, weight = "equal")
#' MixtureMetaB(dsn = EHFE, weight = rep(1, 6))
#' @export
MixtureMetaB <- function(
  dsn    = NULL, 
  weight = NULL, 
  alpha  = 0.05, 
  nboot  = 2000, 
  demog  = NULL
){
  if (tolower(weight)[1] == "equal") {
    W <- rep(1 / nrow(dsn), nrow(dsn))
  } else {
    if (tolower(weight)[1] == "size") {
      W <- with(dsn, (N0 + N1) / sum(N0 + N1))
    } else {
      if (length(weight) == nrow(dsn) & is.numeric(weight)) {
        W <- weight / sum(weight)
      } else {
        stop(paste("The mixture weight needs to be 'equal', 'size' or a numeric vector with the same length as the number of studies."))
      }
    }
  }

  PR1   <- with(dsn, R1 / N1)
  PR0   <- with(dsn, R0 / N0)

  eRR   <- PR1 / PR0
  vleRR <- with(dsn, (N1 - R1) / R1 / N1 + (N0 - R0) / R0 / N0)
  eRRl  <- exp(log(eRR) - qnorm(1 - alpha / 2) * sqrt(vleRR))
  eRRu  <- exp(log(eRR) + qnorm(1 - alpha / 2) * sqrt(vleRR))

  eRD   <- PR1 - PR0
  veRD  <- with(dsn, R1 * (N1 - R1) / N1^3 + R0 * (N0 - R0) / N0^3)
  eRDl  <- eRD - qnorm(1 - alpha / 2) * sqrt(veRD)
  eRDu  <- eRD + qnorm(1 - alpha / 2) * sqrt(veRD)

  IRES <- with(dsn, matrix(c(PR1, PR0, eRR, eRRl, eRRu, eRD, eRDl, eRDu), ncol = 8))
  IRES <- cbind(IRES, W)
  colnames(IRES) <- c("PR1", "PR0",
                      "Risk ratio", paste0(round((1 - alpha) * 100, 0), "% CI_RR_lower"), paste0(round((1 - alpha) * 100, 0), "% CI_RR_upper"),
                      "Risk difference", paste0(round((1 - alpha) * 100, 0), "% CI_RD_lower"), paste0(round((1 - alpha) * 100, 0), "% CI_RD_upper"),
                      "weights")

  RESULT <- cbind(dsn, IRES)

  MP0 <- sum(PR0 * W)
  MP1 <- sum(PR1 * W)

  MRR <- MP1 / MP0

  R1B <- matrix(rbinom(length(dsn$R1) * nboot, size = dsn$N1, prob = PR1), nrow = length(dsn$R1))
  R0B <- matrix(rbinom(length(dsn$R0) * nboot, size = dsn$N0, prob = PR0), nrow = length(dsn$R0))

  mgRRB <- numeric(nboot)

  for (i in 1:nboot){
    mgRRB[i] <- sum(R1B[, i] / dsn$N1 * W) / sum(R0B[, i] / dsn$N0 * W)
  }

  MRRlu <- quantile(mgRRB, probs = c(alpha / 2, 1 - alpha / 2))

  MRD <- MP1 - MP0
  vMRD <- sum(W^2 * with(dsn, R1 * (N1 - R1) / N1^3)) + sum(W^2 * with(dsn, R0 * (N0 - R0) / N0^3))
  MRDl <- MRD - qnorm(1 - alpha / 2) * sqrt(vMRD)
  MRDu <- MRD + qnorm(1 - alpha / 2) * sqrt(vMRD)

  RESULT <- rbind(RESULT[, -1], rep(NA, ncol(dsn) + 8))

  RESULT$PR1[nrow(RESULT)]  <- MP1
  RESULT$PR0[nrow(RESULT)]  <- MP0
  RESULT[nrow(RESULT), c((ncol(RESULT) - 6) : (ncol(RESULT) - 1))] <- c(MRR, MRRlu, MRD, MRDl, MRDu)

  if (length(demog) > 0){
    DEMOV <- paste("RESULT", demog, sep = "$")
    RESV  <- paste0(paste("RESULT", demog, sep = "$"), "[nrow(RESULT)] <- sum(", DEMOV, " * c(W, NA), na.rm = TRUE)")

    for (i in 1 : length(demog)){
      eval(parse(text = RESV[i]))
    }
  }

  NAMES <- colnames(RESULT)
  RESULT$STUDY <- c(dsn[, 1], paste("Mixture weights =", paste(weight, collapse = ",")))
  RESULT <- RESULT[, c("STUDY", NAMES)]

  RESULT
}


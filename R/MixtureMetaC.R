#' Mixture population model approach for meta analysis of continuous outcome
#'
#' Mixture population model approach to combine information across multiple studies for continuous outcomes. This simple and robust procedure based on a mixture population concept can provide a clinically meaningful between-group contrast summary for a well-defined study population.
#' @param dsn Input data frame with study names in the first column, where variable N0 and N1 are the number of subjects (study size) in control group and treatment group, respectively, variable E0 and E1 are the mean of continuous outcome in control group and treatment group, respectively, and variable SD0 and SD1 are the standard deviation of continuous outcome (se = 'FALSE') or the standard error of mean (se = 'TRUE') in control group and treatment group, respectively.
#' @param weight Weights used to combine studies: '"size"', '"equal"' or a numeric vector
#' @param alpha Two-sided significance level
#' @param se 'TRUE' for standard error included in input data frame instead of standard deviation and 'FALSE' (the default) otherwise
#' @param demog A character vector for baseline characteristics included in input data frame
#' @returns A data frame
#' @examples
#' MixtureMetaC(dsn = ELDL, weight = "size", demog = c("Age", "Male", "LDL"))
#' MixtureMetaC(dsn = ELDL, weight = "equal")
#' MixtureMetaC(dsn = ELDL, weight = rep(1, 16))
#' @export
MixtureMetaC <- function(
  dsn    = NULL, 
  weight = NULL, 
  alpha  = 0.05, 
  se     = FALSE, 
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

  if (se == FALSE) {
    V0 <- with(dsn, SD0^2 / N0)
    V1 <- with(dsn, SD1^2 / N1)
  } else {
    V0 <- with(dsn, SD0^2)
    V1 <- with(dsn, SD1^2)
  }

  IRES <- with(dsn, matrix(c(E1 - E0,
                             E1 - E0 - qnorm(1 - alpha / 2) * sqrt(V0 + V1),
                             E1 - E0 + qnorm(1 - alpha / 2) * sqrt(V0 + V1)), ncol = 3))
  IRES <- cbind(IRES, W)
  colnames(IRES) <- c("DIFF", paste0(round((1 - alpha) * 100, 0), "% CI_lower"), paste0(round((1 - alpha) * 100, 0), "% CI_lower"), "weights")

  RESULT <- cbind(dsn, IRES)

  ME0 <- sum(W * dsn$E0)
  ME1 <- sum(W * dsn$E1)

  SE0 <- sqrt(sum(W^2 * V0))
  SE1 <- sqrt(sum(W^2 * V1))

  MD  <- ME1 - ME0
  MDl <- MD - qnorm(1 - alpha / 2) * sqrt(SE0^2 + SE1^2)
  MDu <- MD + qnorm(1 - alpha / 2) * sqrt(SE0^2 + SE1^2)

  RESULT <- rbind(RESULT[, -1], rep(NA, ncol(dsn) + 3))

  RESULT$E1[nrow(RESULT)]  <- ME1
  RESULT$SD1[nrow(RESULT)] <- SE1
  RESULT$E0[nrow(RESULT)]  <- ME0
  RESULT$SD0[nrow(RESULT)] <- SE0
  RESULT[nrow(RESULT), c((ncol(RESULT) - 3) : (ncol(RESULT) - 1))] <- c(MD, MDl, MDu)

  if (length(demog) > 0){
    DEMOV <- paste("RESULT", demog, sep = "$")
    RESV  <- paste0(paste("RESULT", demog, sep = "$"), "[nrow(RESULT)] <- sum(", DEMOV, " * c(W, NA), na.rm = TRUE)")

    for (i in 1:length(demog)){
      eval(parse(text = RESV[i]))
    }
  }

  NAMES <- colnames(RESULT)
  RESULT$STUDY <- c(dsn[, 1], paste("Mixture weights =", paste(weight, collapse = ",")))
  RESULT <- RESULT[, c("STUDY", NAMES)]

  return(RESULT)
}

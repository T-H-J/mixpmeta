% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MixtureMetaC.R
\name{MixtureMetaC}
\alias{MixtureMetaC}
\title{Mixture population model approach for meta analysis of continuous outcome}
\usage{
MixtureMetaC(dsn = NULL, weight = NULL, alpha = 0.05, se = FALSE, demog = NULL)
}
\arguments{
\item{dsn}{Input data frame with study names in the first column, where variable N0 and N1 are the number of subjects (study size) in control group and treatment group, respectively, variable E0 and E1 are the mean of continuous outcome in control group and treatment group, respectively, and variable SD0 and SD1 are the standard deviation of continuous outcome (se = 'FALSE') or the standard error of mean (se = 'TRUE') in control group and treatment group, respectively.}

\item{weight}{Weights used to combine studies: '"size"', '"equal"' or a numeric vector}

\item{alpha}{Two-sided significance level}

\item{se}{'TRUE' for standard error included in input data frame instead of standard deviation and 'FALSE' (the default) otherwise}

\item{demog}{A character vector for baseline characteristics included in input data frame}
}
\value{
A data frame
}
\description{
Mixture population model approach to combine information across multiple studies for continuous outcomes. This simple and robust procedure based on a mixture population concept can provide a clinically meaningful between-group contrast summary for a well-defined study population.
}
\examples{
MixtureMetaC(dsn = ELDL, weight = "size", demog = c("Age", "Male", "LDL"))
MixtureMetaC(dsn = ELDL, weight = "equal")
MixtureMetaC(dsn = ELDL, weight = rep(1, 16))
}

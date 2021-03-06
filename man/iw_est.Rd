% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/iw_est.R
\name{iw_est}
\alias{iw_est}
\title{The inverse weighting estimator (nonparametric method)}
\usage{
iw_est(Y,
       treat,
       treat_formula,
       data,
       grid_val,
       bandw,
       treat_mod,
       link_function,
       ...)
}
\arguments{
\item{Y}{is the the name of the outcome variable contained in \code{data}.}

\item{treat}{is the name of the treatment variable contained in
\code{data}.}

\item{treat_formula}{an object of class "formula" (or one that can be
coerced to that class) that regresses \code{treat} on a linear combination
of \code{X}: a symbolic description of the model to be fitted.}

\item{data}{is a dataframe containing \code{Y}, \code{treat}, and
\code{X}.}

\item{grid_val}{contains the treatment values to be evaluated.}

\item{bandw}{is the bandwidth.  Default is 1.}

\item{treat_mod}{a description of the error distribution to be used in the
model for treatment. Options include: \code{"Normal"} for normal model,
\code{"LogNormal"} for lognormal model, \code{"Sqrt"} for square-root transformation
to a normal treatment, \code{"Poisson"} for Poisson model,
\code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
model.}

\item{link_function}{is either "log", "inverse", or "identity" for the
"Gamma" \code{treat_mod}.}

\item{...}{additional arguments to be passed to the treatment regression function.}
}
\value{
\code{iw_est} returns an object of class "causaldrf",
a list that contains the following components:
\item{param}{parameter estimates for a iw fit.}
\item{t_mod}{the result of the treatment model fit.}
\item{call}{the matched call.}
}
\description{
This is a nonparametric method that estimates the ADRF by using a local linear
regression of \code{Y} on \code{treat} with weighted kernel function.  For
details, see Flores et. al. (2012).
}
\details{
The ADRF is estimated by
\deqn{(D_{0}(t) S_{2}(t) - D_{1}(t) S_{1}(t)) / (S_{0}(t) S_{2}(t) - S_{1}^{2}(t)) }
where \deqn{D_{j}(t) = \sum_{i = 1}^{N} \tilde{K}_{h, X} (T_i - t)  (T_i - t)^j Y_i} and
\eqn{S_{j}(t) = \sum_{i = 1}^{N} \tilde{K}_{h, X} (T_i - t)  (T_i - t)^j}
\eqn{\tilde{K}_{h, X}(t) = K_{h}(t) / \hat{R}_i(t)}
which is a local linear regression.
More details are given in Flores (2012).
}
\examples{
## Example from Schafer (2015).

example_data <- sim_data

iw_list <- iw_est(Y = Y,
                treat = T,
                treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
                data = example_data,
                grid_val = seq(8, 16, by = 1),
                bandw = bw.SJ(example_data$T),
                treat_mod = "Normal")

sample_index <- sample(1:1000, 100)

plot(example_data$T[sample_index],
      example_data$Y[sample_index],
      xlab = "T",
      ylab = "Y",
      main = "iw estimate")

lines(seq(8, 16, by = 1),
        iw_list$param,
        lty = 2,
        lwd = 2,
        col = "blue")

legend('bottomright',
        "iw estimate",
        lty=2,
        lwd = 2,
        col = "blue",
        bty='Y',
        cex=1)

rm(example_data, iw_list, sample_index)

## Example from Imai & van Dyk (2004).

data("nmes_data")
head(nmes_data)
# look at only people with medical expenditures greater than 0
nmes_nonzero <- nmes_data[which(nmes_data$TOTALEXP > 0), ]


iw_list <- iw_est(Y = TOTALEXP,
                  treat = packyears,
                  treat_formula = packyears ~ LASTAGE + I(LASTAGE^2) +
                    AGESMOKE + I(AGESMOKE^2) + MALE + RACE3 + beltuse +
                    educate + marital + SREGION + POVSTALB,
                  data = nmes_nonzero,
                  grid_val = seq(5, 100, by = 5),
                  bandw = bw.SJ(nmes_nonzero$packyears),
                  treat_mod = "LogNormal")

set.seed(307)
sample_index <- sample(1:nrow(nmes_nonzero), 1000)

plot(nmes_nonzero$packyears[sample_index],
     nmes_nonzero$TOTALEXP[sample_index],
     xlab = "packyears",
     ylab = "TOTALEXP",
     main = "iw estimate",
     ylim = c(0, 10000),
     xlim = c(0, 100))

lines(seq(5, 100, by = 5),
      iw_list$param,
      lty = 2,
      lwd = 2,
      col = "blue")

legend('topright',
       "iw estimate",
       lty=2,
       lwd = 2,
       col = "blue",
       bty='Y',
       cex = 1)
abline(0, 0)
}
\references{
Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
continuous treatment and outcome: alternative estimators for parametric
dose-response models. \emph{Manuscript in preparation}.

Flores, Carlos A., et al. "Estimating the effects of length of exposure to
instruction in a training program: the case of job corps."
\emph{Review of Economics and Statistics} \bold{94.1} (2012): 153-171.
}
\seealso{
\code{\link{nw_est}}, \code{\link{iw_est}}, \code{\link{hi_est}}, \code{\link{gam_est}},
\code{\link{add_spl_est}},
\code{\link{bart_est}}, etc. for other estimates.
}


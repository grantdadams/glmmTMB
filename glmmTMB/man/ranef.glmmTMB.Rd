% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{ranef.glmmTMB}
\alias{ranef.glmmTMB}
\alias{ranef}
\alias{as.data.frame.ranef.glmmTMB}
\title{Extract Random Effects}
\usage{
\method{ranef}{glmmTMB}(object, condVar = TRUE, ...)

\method{as.data.frame}{ranef.glmmTMB}(x, ...,
  stringsAsFactors = default.stringsAsFactors())
}
\arguments{
\item{object}{a \code{glmmTMB} model.}

\item{condVar}{include conditional variances in result?}

\item{...}{some methods for this generic function require additional
arguments.}

\item{x}{a \code{ranef.glmmTMB} object (i.e., the result of running \code{ranef} on a fitted \code{glmmTMB} model)}

\item{stringsAsFactors}{see \code{\link{data.frame}}}
}
\value{
\itemize{
\item For \code{ranef}, an object of class \code{ranef.glmmTMB} with two components:
\describe{
  \item{cond}{a list of data frames, containing random effects
    for the conditional model.}
  \item{zi}{a list of data frames, containing random effects for
    the zero inflation.}
}
If \code{condVar=TRUE}, the individual list elements within the
\code{cond} and \code{zi} components (corresponding to individual
random effects terms) will have associated \code{condVar} attributes
giving the conditional variances of the random effects values.
These are in the form of three-dimensional arrays: see
\code{\link{ranef.merMod}} for details (the only difference between
the packages is that the attributes are called \sQuote{"postVar"}
in \pkg{lme4}, vs. \sQuote{"condVar"} in \pkg{glmmTMB}.
\item For \code{as.data.frame}: a data frame with components
\describe{
\item{component}{part of the model to which the random effects apply (conditional or zero-inflation)}
\item{grpvar}{grouping variable}
\item{term}{random-effects term (e.g., intercept or slope}
\item{grp}{group, or level of the grouping variable}
\item{condval}{value of the conditional mode}
\item{condsd}{conditional standard deviation}
}
}
}
\description{
Generic function to extract random effects from \code{glmmTMB} models, both
for the conditional model and zero inflation.
}
\note{
When a model has no zero inflation, the default behavior of
  \code{ranef} is to simplify the printed format of the random effects. To
  show the full list structure, run \code{print(ranef(model),
  simplify=FALSE)}. In all cases, the full list structure is used to access
  the data frames (see example).
}
\examples{
if (requireNamespace("lme4")) {
   data(sleepstudy, package="lme4")
   model <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy)
   rr <- ranef(model)
   print(rr, simplify=FALSE)
   ## extract Subject conditional modes for conditional model
   rr$cond$Subject
   as.data.frame(rr)
}
}
\seealso{
\code{\link{fixef.glmmTMB}}.
}

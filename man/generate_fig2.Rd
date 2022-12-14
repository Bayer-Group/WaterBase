% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAS_LOQ.R
\name{generate_fig2}
\alias{generate_fig2}
\title{Generate figure 2 (Figure 3 in Manuscript)}
\usage{
generate_fig2(
  Res,
  cutoff = 0,
  df = data.frame(x = 1e-06, y = 8),
  test = "Chronic",
  used = c("Mean"),
  exclude.cas = NULL,
  chemclass,
  ...
)
}
\arguments{
\item{Res}{processed waterbase data}

\item{cutoff}{0, not excluding sites}

\item{df}{indication of plotting legends position in the figure}

\item{test}{"Chronic" or "Acute"}

\item{used}{"Mean" or "Max"}

\item{exclude.cas}{whether to exclude some chemicals}

\item{chemclass}{chemical information}

\item{...}{other to be passed arguments}
}
\value{

}
\description{
Generate figure 2 (Figure 3 in Manuscript)
}

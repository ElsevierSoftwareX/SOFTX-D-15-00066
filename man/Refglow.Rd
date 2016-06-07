\name{Refglow}
\alias{Refglow}
\docType{data}
\title{
Reference glow curves
}

\description{
Two thermoluminescence glow curves.
}

\usage{data(Refglow)}

\format{
A list that contains two thermoluminescence glow curves:
\describe{
    \item{Refglow002}{a synthetic glow curve from the GLOCANIN project (Bos et al., 1993)}
    \item{Refglow009}{an experimental glow-curve from the GLOCANIN project (Bos et al., 1994)}
}
}

\details{
\bold{Refglow} contains two thermoluminescence glow curves (Refglow002 and Refglow009) from the GLOCANIN project.
The author is grateful to Adrie J.J. Bos for providing the reference glow curves of the GLOCANIN project to test function \link{tgcd}.
}

\references{
Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1993. An intercomparison of glow curve analysis computer programs: I. Synthetic glow curves. Radiation Protection Dosimetry, 47(1-4), 473-477.

Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1994. An intercomparison of glow curve analysis computer programs: II. Measured glow curves. Radiation Protection Dosimetry, 51(4): 257-264.
}

\examples{
 # Load package "tgcd".
  require(tgcd)

  data(Refglow)
  names(Refglow)
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve simulation}
\keyword{Glow curve deconvolution}

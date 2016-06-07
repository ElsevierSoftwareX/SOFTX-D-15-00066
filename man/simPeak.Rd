\name{simPeak}
\alias{simPeak}
\alias{simPeak.default}
\title{
Thermoluminescence glow peak simulation
}
\description{
Simulating first-order, second-order, or general-order glow peaks.
}
\usage{
simPeak(temps, n0, Nn = NULL, bv = NULL, ff, 
        ae, hr, typ = c("f", "s", "g"), 
        outfile = NULL, plot = TRUE)
}

\arguments{
  \item{temps}{\link{vector}(\bold{required}): temperature values (K) where the values of the thermoluminescence 
intensity will be computed. It needs to be sorted increasingly. A vector of temperature values may be generated 
using the internal function \link{seq}}

  \item{n0}{\link{numeric}(\bold{required}): initial concentration of trapped electrons (1/cm^3)}

  \item{Nn}{\link{numeric}(\bold{required}): total concentration of the traps in the crystal (1/cm^3)}

  \item{bv}{\link{numeric}(\bold{required}): order number for the general order glow peak}

  \item{ff}{\link{numeric}(\bold{required}): the frequency factor (1/s)}

  \item{ae}{\link{numeric}(\bold{required}): the activation energy (eV)}

  \item{hr}{\link{numeric}(with default): the linear heating rate (K/s)}

  \item{typ}{\link{character}(with default): the type of a glow peak, \code{typ="f"} means first-order, \code{typ="s"} means second-order, \code{typ="g"} means general-order,
default \code{typ="f"}}

  \item{outfile}{\link{character}(optional): if specified, simulated intensities of glow peaks will 
be written to a file named \code{"outfile"} in \code{CSV} format and saved to the current work directory}

  \item{plot}{\link{logical}(with default): draw a plot according to the simulated glow peak or not}
}

\details{
Function \emph{simPeak} simulates glow peaks of various orders. The first-, second-, and general-order glow peak 
can be simulated using the following three ordinary equations, respectively (Pagonis et al., 2006):

\eqn{\frac{dn}{dT}=-n*S*e^{-E/(k*T)}/hr}

\eqn{\frac{dn}{dT}=-n^{2}*S*e^{-E/(k*T)}/(Nn*hr)}

\eqn{\frac{dn}{dT}=-n^{bv}*S*e^{-E/(k*T)}/(Nn*hr)}

where \eqn{n} is the concentration of trapped electrons, \eqn{\frac{dn}{dT}} the rate of change of the concentration 
of trapped electrons, \eqn{S} the frequency factor, \eqn{E} the activation energy, \eqn{T} the absolute temperature, 
\eqn{k} the Boltzmann constant, \eqn{Nn} the total concentration of the traps in the crystal, \eqn{bv} the b value 
(kinetic order), and \eqn{hr} the linear heating rate.

The ordinary equation is solved by the Fortran 77 subroutine \emph{lsoda} (original version written by Linda R. Petzold and 
Alan C. Hindmarsh available at Netlib: \url{http://www.netlib.org/odepack/}, modified version by R. Woodrow Setzer from the R package 
\emph{deSolve} (Soetaert et al., 2010) available at CRAN: \url{http://CRAN.R-project.org/package=deSolve}).
}

\value{
Return an invisible list containing the following elements:
    \item{temps}{a vector of temperature values}
    \item{tl}{values of the thermoluminescence intensity}
    \item{n}{variation of concentration of trapped electrons with temperature}
    \item{sp}{parameters used for describing the shape of a glow peak (Pagonis et al., 2006):\cr 
              the temperature corresponding to half intensity on the left side of the peak (\code{T1});\cr 
              the temperature corresponding to half intensity on the right side of the peak (\code{T2});\cr
              the temperature corresponding to maximum intensity (\code{Tm});\cr
              the half-width at the left side of the peak (\code{d1=Tm-T1});\cr
              the half-width at the right side of the peak (\code{d2=T2-Tm});\cr
              the total half-width (\code{thw=d1+d2});\cr
              the symmetry factor (\code{sf=d2/thw})}
}

\references{
Pagonis V, Kitis G, Furetta C, 2006. Numerical and practical exercises in thermoluminescence. Springer Science & Business Media.

Soetaert K, Petzoldt T, Setzer RW, 2010. Solving Differential Equations in R: Package deSolve. Journal of Statistical Software, 33(9): 1-25.
}

\seealso{
\link{tgcd}; \link{simqOTOR}
}

\examples{
   ### Effect of initial electron trap concentration (n0) 
   ### on second-order glow peaks.
    temps <- seq(400, 600, by=0.5)
    peak1 <- simPeak(temps, n0=0.2e10, Nn=1e10, 
      ff=1e19, ae=2.0, hr=1, typ="s")
    peak2 <- simPeak(temps, n0=0.4e10, Nn=1e10, 
      ff=1e19, ae=2.0, hr=1, typ="s")
    peak3 <- simPeak(temps, n0=0.6e10, Nn=1e10, 
      ff=1e19, ae=2.0, hr=1, typ="s")
    peak4 <- simPeak(temps, n0=0.8e10, Nn=1e10, 
      ff=1e19, ae=2.0, hr=1, typ="s")
    peak5 <- simPeak(temps, n0=1.0e10, Nn=1e10, 
      ff=1e19, ae=2.0, hr=1, typ="s")
    peaks <- cbind(peak1$tl, peak2$tl, peak3$tl, peak4$tl, peak5$tl)
    matplot(temps, peaks, type="l", lwd=2, lty="solid", 
      xlab="Temperature (K)", ylab="TL intensity (counts)")

}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow peak simulation}

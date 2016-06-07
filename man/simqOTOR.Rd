\name{simqOTOR}
\alias{simqOTOR}
\alias{simqOTOR.default}
\title{
Thermoluminescence glow peak simulation
}
\description{
Simulating glow peaks according to the one trap-one recombination 
center (OTOR) model using the quasi-equilibrium approximation.
}
\usage{
simqOTOR(temps, n0, Nn, Ah, An, ff, ae, 
         hr, outfile = NULL, plot = TRUE)
}

\arguments{
  \item{temps}{\link{vector}(\bold{required}): temperature values (K) where the values of the thermoluminescence 
intensity will be computed, it needs to be sorted increasingly}

  \item{n0}{\link{numeric}(\bold{required}): initial concentration of trapped electrons (1/cm^3)}

  \item{Nn}{\link{numeric}(\bold{required}): total concentration of the traps in the crystal (1/cm^3)}

  \item{Ah}{\link{numeric}(optional): probability coefficient of electron recombining with holes in the recombination center (cm^3/s)}

  \item{An}{\link{numeric}(optional): probability coefficient of electron retrapping in the traps (cm^3/s)}

  \item{ff}{\link{numeric}(\bold{required}): the frequency factor (1/s)}

  \item{ae}{\link{numeric}(\bold{required}): the activation energy (eV)}

  \item{hr}{\link{numeric}(with default): the linear heating rate (K/s)}

  \item{outfile}{\link{character}(optional): if specified, simulated intensities of glow peaks will 
be written to a file named \code{"outfile"} in \code{CSV} format and saved to the current work directory}

  \item{plot}{\link{logical}(with default): draw a plot according to the simulated glow peak or not}
}

\details{
Function \emph{simqOTOR} simulates a synthetic glow peak according to the OTOR model using the quasi-equilibrium 
approximation. This function may be used to simulating glow peaks of first-, second-, and general-order, depending 
on the given kinetic parameters. The approximate equation of the OTOR model derived using the quasi-equilibrium 
approximation can be described by (Pagonis et al., 2006):

\eqn{\frac{dn}{dT}=-Ah*n^2*S*e^{-E/(kT)}/((n*Ah+(Nn-n)*An)*hr)}

where \eqn{n} is the concentration of trapped electrons, \eqn{\frac{dn}{dT}} the rate of change of the concentration 
of trapped electrons, \eqn{S} the frequency factor, \eqn{E} the activation energy, \eqn{T} the absolute temperature, 
\eqn{k} the Boltzmann constant, \eqn{Nn} the total concentration of the traps in the crystal, \eqn{Ah} the probability 
coefficient of electron recombining with holes in the recombination center, \eqn{An} the probability coefficient of 
electron retrapping in the traps, and \eqn{hr} the linear heating rate.

The ordinary equation is solved by the Fortran 77 subroutine \emph{lsoda} (original version written by Linda R. Petzold and 
Alan C. Hindmarsh available at Netlib: \url{http://www.netlib.org/odepack/}, modified version by R. Woodrow Setzer from the R package 
\emph{deSolve} (Soetaert et al., 2010) available at CRAN: \url{http://CRAN.R-project.org/package=deSolve}).
}

\value{
Return an invisible list containing the following elements:
    \item{temps}{a vector of temperature values}
    \item{tl}{values of the thermoluminescence intensity}
    \item{n}{variation of concentration of trapped electrons with temperature}
    \item{sp}{parameters used for describing the shape of a glow peak, see function \link{simPeak} for details}
}

\references{
Pagonis V, Kitis G, Furetta C, 2006. Numerical and practical exercises in thermoluminescence. Springer Science & Business Media.

Soetaert K, Petzoldt T, Setzer RW, 2010. Solving Differential Equations in R: Package deSolve. Journal of Statistical Software, 33(9): 1-25. 
}

\seealso{
\link{tgcd}; \link{simPeak}
}

\examples{
    # Synthesizing a glow curve consisting of five glow peaks.
    temps <- seq(330, 730, by=0.5)
    peak1 <- simqOTOR(temps, n0=0.7e10, Nn=1e10, Ah=1e-3, An=1e-7,
      ff=1e14, ae=1.5, hr=1, outfile = NULL, plot = TRUE)
    peak2 <- simqOTOR(temps, n0=0.5e10, Nn=1e10, Ah=1e-7, An=1e-7,
      ff=1e17, ae=1.9, hr=1, outfile = NULL, plot = TRUE)
    peak3 <- simqOTOR(temps, n0=0.2e10, Nn=1e10, Ah=1e-5, An=1e-7,
      ff=1e15, ae=1.45, hr=1, outfile = NULL, plot = TRUE)
    peak4 <- simqOTOR(temps, n0=0.2e10, Nn=1e10, Ah=1e-5, An=1e-7,
      ff=1e9, ae=0.85, hr=1, outfile = NULL, plot = TRUE)
    peak5 <- simqOTOR(temps, n0=0.3e10, Nn=1e10, Ah=1e-7, An=1e-7,
      ff=1e11, ae=1.4, hr=1, outfile = NULL, plot = TRUE)
    peaks <- cbind(peak1$tl, peak2$tl, peak3$tl, peak4$tl, peak5$tl, 
               peak1$tl+peak2$tl+peak3$tl+peak4$tl+peak5$tl)
    matplot(temps, y=peaks, type="l", lwd=2, lty="solid", 
      xlab="Temperature (K)", ylab="TL intensity (counts)")
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow peak simulation}

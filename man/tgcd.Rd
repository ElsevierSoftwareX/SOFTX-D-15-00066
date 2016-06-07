\name{tgcd}
\alias{tgcd}
\alias{tgcd.default}
\title{
Thermoluminescence glow curve deconvolution (TGCD)
}
\description{
Thermoluminescence glow curve deconvolution according to the general-order empirical expression or the semi-analytical expression derived from the one trap-one recombination (OTOR) model based on the Lambert's W function.
}
\usage{
tgcd(Sigdata, npeak, inis = NULL, mwt = 90, 
     mdt = 3, nstart = 30, model=c("g", "lw"), 
     elim = NULL, logy = FALSE, hr = NULL, 
     outfile = NULL, plot = TRUE)
}

\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a \code{two-column} matrix, temperature values and thermoluminescence signal values are stored in the first and second column, respectively}

  \item{npeak}{\link{integer}(\bold{required}): number of glow peaks, the allowed maximum number of glow peaks is set equal to 13}

  \item{inis}{\link{matrix}(optional): a \code{npeak-row} \code{4-column} matrix used for storing starting parameters \code{Im}, \code{E}, \code{Tm}, \code{b} (or \code{R}). 
If \code{inis=NULL}, the user has to click with a mouse on a plot to locate each glow peak maximum}

\item{mwt}{\link{numeric}(with default): allowed maximum total half-width of deconvoluted glow peaks. A smaller \code{mwt} prevents 
the appearance of glow peaks with large total half-width. A conservative \code{mwt} is adopted by default}

  \item{mdt}{\link{numeric}(with default): allowed minimum distance between each optimized temperature at maximum thermoluminescence intensity. A larger \code{mdt} prevents 
the appearance of strongly overlapping peaks. A conservative \code{mdt} is adopted by default}

  \item{nstart}{\link{integer}(with default): number of trials, if \code{nstart>1}, a "try-and-error" protocol will be performed, the upper limit for
\code{nstart} is set equal to 10000}

  \item{model}{\link{character}(with default): \code{"g"} means fitting a general-order model, \code{"lw"} means fitting the Lambert W-function, default \code{model="g"}}

  \item{elim}{\link{vector}(with default): lower and upper limits for the activation energy, default \code{elim=c(0.5, 5)} }

  \item{logy}{\link{logical}(with default): draw the y-axis of the plot used for locating peak maxima with a logarithmic scale or not}

  \item{hr}{\link{numeric}(optional): the linear heating rate used for calculating the frequency factor}

  \item{outfile}{\link{character}(optional): if specified, fitted signal values for each glow peak will be written to a file named \code{"outfile"} in \code{CSV} format and 
saved to the current work directory}

  \item{plot}{\link{logical}(with default): draw a plot according to the fitting result or not}
}

\details{
Function \emph{tgcd} is used for deconvolving thermoluminescence glow curves according to the general-order empirical expression 
(Kitis et al., 1998; Pagonis et al., 2006) or the semi-analytical expression derived from the one trap-one recombination (OTOR) 
model based on the Lambert's W function (Kitis and Vlachos, 2013; Sadek et al., 2015; Kitis et al., 2016) using the Levenberg-Marquardt 
algorithm (plus supports for constraining and fixing parameters).

The general-order empirical expression for a glow peak is:

\code{I(T)=Im*b^(b/(b-1))*expv*((b-1)*(1-xa)*(T/Tm)^2*expv+Zm)^(-b/(b-1))}\cr
\code{xa=2*k*T/E}\cr
\code{xb=2*k*Tm/E}\cr
\code{expv=exp(E/(k*T)*(T-Tm)/Tm)}\cr
\code{Zm=1+(b-1)*xb}

where \code{b} is the kinetic parameter (lies between 1 and 2), \code{I} is the glow peak intensity, 
\code{E} the activation energyin ev, \code{k} the Boltzmann constant in eV/k, T the temperature in K with constant heating rate K/s, 
\code{Tm} the temperature at maximum thermoluminescence intensity in K, and \code{Im} the maximum intensity.
The four parameters for this model are: \code{Im}, \code{E}, \code{Tm}, and \code{b}.\cr

The semi-analytical expression derived from the one trap-one recombination (OTOR) model based on the Lambert's W function is:

\code{I(T)=Im*exp(-E/(k*T)*(Tm-T)/Tm)*(W(Zm)+W(Zm)^2)/(W(Z)+W(Z)^2)}\cr
\code{Zm=R/(1-R)-log((1-R)/R)+E*exp(E/(k*Tm))/(k*Tm^2*(1-1.05*R^1.26))*F(Tm,E)}\cr
\code{Z=R/(1-R)-log((1-R)/R)+E*exp(E/(k*Tm))/(k*Tm^2*(1-1.05*R^1.26))*F(T,E)}\cr
\code{F(Tm,E)=Tm*exp(-E/(k*Tm))+E/k*Ei(-E/(k*Tm))}\cr
\code{F(T,E)=T*exp(-E/(k*T))+E/k*Ei(-E/(k*T))}

where \code{W(x)} is the wright Omega function for variable x, \code{Ei(x)} is the exponential integral function for variable x, 
\code{I} is the glow peak intensity, \code{E} the activation energy in eV, \code{k} the Boltzmann constant in eV/K, 
\code{T} the temperature in K with constant heating rate in K/s, \code{Tm} the temperature at maximum thermoluminescence intensity in K, 
and \code{Im} the maximum intensity. The four parameters for this model are: \code{Im}, \code{E}, \code{Tm}, and \code{R}.
The Fortran 90 subroutine used for evaluating the Wright Omega function is transformed from the Matlab code 
provided by Andrew Horchler available at \url{https://github.com/horchler/wrightOmegaq}.\cr

The procedure minimizes the objective:

\eqn{fcn=\sum_{i=1}^n |y_i^o-y_i^f|, i=1,...,n}

where \eqn{y_i^o} and \eqn{y_i^f} denote the \emph{i-th} observed and fitted signal value, respectively, 
and \eqn{n} indicates the number of data points.

Starting parameters (\code{inis}) can be specified by the user through argument \code{inis} or by clicking with a mouse on 
the plot of the thermoluminescence glow curve showing peak maxima if \code{inis=NULL}.The Levenberg-Marquardt algorithm 
(More, 1978) (\emph{minpack}: Original Fortran 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. Fortran 90 version 
by John Burkardt freely available at \url{http://people.sc.fsu.edu/~jburkardt/f_src/minpack/}) was modified so as to supports 
constraints and fixes of parameters. If argument \code{nstart>1}, a "try-and-error" protocol with starting values generated 
uniformly around the given starting values \code{inis} will be performed repeatedly to search the optimal parameters that 
give a minimum Figure Of Merit (FOM) value.

Parameters can be interactively constrained and fixed by modifying the following elements in a automatically 
generated \bold{Dialog Table} if \code{inis=NULL}:\cr 
\bold{(1)} \emph{INTENS(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{Im} \cr
\bold{(2)} \emph{ENERGY(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{E} \cr
\bold{(3)} \emph{TEMPER(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{Tm} \cr
\bold{(4)} \emph{bValue(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{b} \cr
}

\value{
Return a list containing the following elements:
   \item{pars}{optimized parameters stored in a matrix}
   \item{ff}{calculated frequency factor if \code{hr!=NULL}}
   \item{sp}{parameters used for describing the shape of a glow peak, see function \link{simPeak} for details}
   \item{FOM}{minimized Figure Of Merit}
}

\references{
Kitis G, Gomes-Ros JM, Tuyn JWN, 1998. Thermoluminescence glow curve deconvolution functions for first, second and general orders of kinetics. Journal of Physics D: Applied Physics, 31(19): 2636-2641.

Kitis G, Polymeris GS, Sfampa IK, Prokic M, Meric N, Pagonis V, 2016. Prompt isothermal decay of thermoluminescence in MgB4O7:Dy, Na and LiB4O7:Cu, In dosimeters. Radiation Measurements, 84: 15-25.

Kitis G, Vlachos ND, 2013. General semi-analytical expressions for TL, OSL and other luminescence stimulation modes derived from the OTOR model using the Lambert W-function. Radiation Measurements, 48: 47-54.

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory," in Lecture Notes in Mathematics: Numerical Analysis, Springer-Verlag: Berlin. 105-116.

Pagonis V, Kitis G, Furetta C, 2006. Numerical and practical exercises in thermoluminescence. Springer Science & Business Media. 

Sadek AM, Eissa HM, Basha AM, Carinou E, Askounis P, Kitis G, 2015. The deconvolution of thermoluminescence glow-curves using general expressions derived from the one trap-one recombination (OTOR) level
model. Applied Radiation and Isotopes, 95: 214-221.

\bold{Further reading} 

Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1993. An intercomparison of glow curve analysis computer programs: I. Synthetic glow curves. Radiation Protection Dosimetry, 47(1-4), 473-477.

Chung KS, Choe HS, Lee JI, Kim JL, Chang SY, 2005. A computer program for the deconvolution of thermoluminescence glow curves. Radiation Protection Dosimetry, 115(1-4): 345-349.
Software is freely available at \url{http://physica.gsnu.ac.kr/TLanal}.

Harvey JA, Rodrigues ML, Kearfott JK, 2011. A computerized glow curve analysis (GCA) method for WinREMS thermoluminescent dosimeter data using MATLAB. Applied Radiation and Isotopes, 69(9):1282-1286.
Source codes are freely available at \url{http://www.sciencedirect.com/science/article/pii/S0969804311002685}.

Kiisk V, 2013. Deconvolution and simulation of thermoluminescence glow curves with Mathcad. Radiation Protection Dosimetry, 156(3): 261-267.
Software is freely available at \url{http://www.physic.ut.ee/~kiisk/mcadapps.htm}.

Puchalska M, Bilski P, 2006. GlowFit-a new tool for thermoluminescence glow-curve deconvolution. Radiation Measurements, 41(6): 659-664. 
Software is freely available at \url{http://www.ifj.edu.pl/dept/no5/nz58/deconvolution.htm}.
}

\note{
The model to be optimized should not be underdetermined. This means that the number of data points should exceed the number of parameters. 
A lack of background counts in the analyzed data is assumed. To obtain reliable estimate, the presented background may be accounted for 
by subtracting from measured data before analysis.
}

\seealso{
\link{simPeak}; \link{simqOTOR}
}

\examples{
# Load the data.
  data(Refglow)

# Deconvolve Refglow002 with 4 peaks using the Lambert W-function.
  startingPars <- 
  cbind(c(400, 550, 850, 1600), # Im
        c(1.4, 1.5, 1.6, 2),    # E
        c(420, 460, 480, 510),  # Tm
        c(0.1, 0.1, 0.1, 0.1))  # R
  tgcd(Refglow$Refglow002, npeak=4, model="lw",
       inis=startingPars, nstart=10)

# Do not run.
# Deconvolve Refglow009 with 9 peaks using the general-order equation.
# startingPars <- 
# cbind(c(9824, 21009, 27792, 50520, 7153, 5496, 6080, 1641, 2316), # Im
#       c(1.24, 1.36, 2.10, 2.65, 1.43, 1.16, 2.48, 2.98, 2.25),    # E
#       c(387, 428, 462, 488, 493, 528, 559, 585, 602),             # Tm
#       c(1.02, 1.15, 1.99, 1.20, 1.28, 1.19, 1.40, 1.01, 1.18))    # b
# tgcd(Refglow$Refglow009, npeak=9, model="g",
#      inis=startingPars, nstart=10)
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve deconvolution}

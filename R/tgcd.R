###
tgcd <-
function(Sigdata, npeak, inis=NULL, mwt=90,
         mdt=3, nstart=30, model=c("g","lw"),
         elim=NULL, logy=FALSE, hr=NULL, 
         outfile=NULL, plot=TRUE)  {
    UseMethod("tgcd")
} #
### 2016.05.22.
tgcd.default <- 
function(Sigdata, npeak, inis=NULL, mwt=90, 
         mdt=3, nstart=30, model=c("g","lw"),
         elim=NULL, logy=FALSE, hr=NULL, 
         outfile=NULL, plot=TRUE)  {
        ### Stop if not.
        stopifnot(ncol(Sigdata)==2L, 
                  all(Sigdata[,1L,drop=TRUE]>0),
                  all(Sigdata[,2L,drop=TRUE]>=0),
                  length(npeak)==1L, is.numeric(npeak), 
                  npeak %in% seq(13L), 3L*npeak<nrow(Sigdata), 
                  is.null(inis) || is.matrix(inis),
                  length(mwt)==1L, is.numeric(mwt), mwt>0,
                  length(mdt)==1L, is.numeric(mdt), mdt>0,
                  length(nstart)==1L, is.numeric(nstart), nstart>=1L, nstart<=10000L,
                  length(model) %in% c(1L, 2L), all(model %in% c("g", "lw")),
                  is.null(elim) || is.numeric(elim),
                  is.logical(logy), length(logy)==1L,
                  is.null(hr) || is.numeric(hr),
                  is.null(outfile) || is.character(outfile),
                  is.logical(plot), length(plot)==1L)
        ###
        if (!is.null(inis))  {
            if(!is.numeric(inis))  stop("Error: inis should be a numeric matrix!")
            if(dim(inis)[1L]!=npeak) stop("Error: incorrect dimensions of inis!")
            if(model[1L]=="lw" && dim(inis)[2L]!=4L) stop("Error: incorrect dimensions of inis!")
            if(model[1L]=="g" && dim(inis)[2L]!=4L) stop("Error: incorrect dimensions of inis!")
        } # end if.
        if (!is.null(elim)) {
            if(length(elim)!=2L) stop("Error: elim should be a two-element vector!")
            if(elim[1L]>=elim[2L]) stop("Error: invalid elim!")
            if(elim[1L]>=1.8) stop("Error: lower limit of elim is too large!")
            if(elim[2L]<=2.2) stop("Error: upper limit of elim is too small!")
        } # end if.
        if (!is.null(hr)) {
            if(length(hr)!=1L) stop("Error: hr should be a one-element vector!")
            if(hr<=0) stop("Error: hr should be larger than zero!")
        } # end if.
        if (!is.null(outfile)) {
            if(length(outfile)!=1L) stop("Error: outfile should be an one-element vector!")
        } # end if.
        ###
        ### Temperature and signal values.
        temp <- as.numeric(Sigdata[,1L,drop=TRUE])
        signal <- as.numeric(Sigdata[,2L,drop=TRUE])
        ###
        if(is.null(inis))  {
            abzero <- which(signal>(.Machine$double.eps)^0.3)
            par(mar=c(4.5,4.5,4,2)+0.1)
            plot(temp[abzero], signal[abzero], type="l", col="skyblue3", 
                 lwd=5, xlab="Temperature (K)", log=ifelse(logy,"y",""), 
                 ylab="TL intensity (counts)", main=paste("Click the mouse to select ", 
                 npeak, " peak maxima:", sep=""), cex.lab=1.3)
            grid(col="darkviolet", lwd=1)
            par(mar=c(5,4,4,2)+0.1)
            sldxy <- try(locator(n=npeak), silent=TRUE)
            if(class(sldxy)=="try-error") {
                stop("Error: no available starting values!")
            } else {
                sldxy_index <- order(sldxy$x, decreasing=FALSE)
                sldxy$x <- sldxy$x[sldxy_index]
                sldxy$y <- sldxy$y[sldxy_index]
            } # end if.
        } # end if.
        ###
        minTEMPER <- min(temp)
        maxTEMPER <- max(temp)
        minINTENS <- min(signal)
        maxINTENS <- max(signal)
        ###
        ### Function used for setting initial
        ### parameters manually (interactively).
        setpars <-function(npeak, sldx, sldy)  {
            mat1 <- mat2 <- mat3 <- mat4 <-
            as.data.frame(matrix(nrow=npeak+1L, ncol=5L))
            ###
            ### Default TL growth peak intensity. 
            mat1[1L,] <-  c("Peak", "INTENS(min)", "INTENS(max)",  
                            "INTEN(ini)", "INTENS(fix)")
            mat1[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat1[-1L,2L] <-  round(minINTENS*0.8, 3L)
            mat1[-1L,3L] <-  round(maxINTENS*1.2, 3L) 
            mat1[-1L,4L] <- if (is.null(inis)) {
                round(sldy, 3L)
            } else {
                round(inis[,1L,drop=TRUE], 3L)
            } # end if.
            mat1[-1L,5L] <- FALSE
            ###
            ### Default activation energy.
            mat2[1L,]<- c("Peak", "ENERGY(min)", "ENERGY(max)", 
                          "ENERGY(ini)", "ENERGY(fix)")
            mat2[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat2[-1L,2L] <-  if (is.null(elim)) 0.5 else round(elim[1L], 3L)
            mat2[-1L,3L] <-  if (is.null(elim)) 5.0 else round(elim[2L], 3L)
            mat2[-1L,4L] <- if (is.null(inis)) {
                round(runif(n=npeak,min=1.8, max=2.2), 3L)
            } else {
                round(inis[,2L,drop=TRUE], 3L)
            } # end if.
            mat2[-1L,5L] <- FALSE
            ###
            ### Default temperature at the peak maximum.
            mat3[1L,] <- c("Peak", "TEMPER(min)", "TEMPER(max)",  
                           "TEMPER(ini)", "TEMPER(fix)")
            mat3[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat3[-1L,2L] <-  round(minTEMPER, 3L)
            mat3[-1L,3L] <-  round(maxTEMPER, 3L)
            mat3[-1L,4L] <- if (is.null(inis)) {
                round(sldx, 3L)
            } else {
                round(inis[,3L,drop=TRUE], 3L)
            } # end if.
            mat3[-1L,5L] <- FALSE
            ###
            ### Default bValue or rValue for a glow peak.
            if (model[1L]=="g")  {
                mat4[1L,] <- c("Peak", "bValue(min)", "bValue(max)",  
                               "bValue(ini)", "bValue(fix)")
                mat4[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
                mat4[-1L,2L] <- 1.0
                mat4[-1L,3L] <-  2.0
                mat4[-1L,4L] <- if (is.null(inis)) {
                    round(runif(n=npeak,min=1.1, max=1.9), 3L)
                } else {
                    round(inis[,4L,drop=TRUE], 3L)
                } # end if.
                mat4[-1L,5L] <- FALSE
            } else if(model[1L]=="lw") {
                mat4[1L,] <- c("Peak", "rValue(min)", "rValue(max)",  
                               "rValue(ini)", "rValue(fix)")
                mat4[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
                mat4[-1L,2L] <- 1.0e-13
                mat4[-1L,3L] <-  1.0-1.0e-13
                mat4[-1L,4L] <- if (is.null(inis)) {
                    round(runif(n=npeak,min=0.1, max=0.9), 3L)
                } else {
                    inis[,4L,drop=TRUE]
                } # end if.
                mat4[-1L,5L] <- FALSE
            } # end if.
            ###
            ###
            mat <- rbind(mat1, rep("    ", 5L),
                         mat2, rep("    ", 5L), 
                         mat3, rep("    ", 5L),
                         mat4)       
            ###
            if(is.null(inis)) {
                pars <- try(edit(name=mat), silent=TRUE)
                if (class(pars)=="try-error") {
                    stop("Error: incorrect parameter modification!")
                } # end if.
            } else {
                pars <- mat
            } # end if.
            ###
            return(pars)
        } # end function setpars.
       ###
       ###
       ### cat("Set parameter constraints:\n")
       pars <- setpars(npeak=npeak, sldx=sldxy$x, sldy=sldxy$y)
       indx <- seq(from=2L, to=npeak+1L, by=1L)
       ###
       ###############################################################
       ###############################################################
       ### 1. TL growth peak intensity.
       ### Check non-finite vlaues.
       intensity1 <- as.numeric(pars[indx,2L,drop=TRUE])
       whichloc <- which(!is.finite(intensity1))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite lower bound of INTENS")
       } # end if.
       intensity2 <- as.numeric(pars[indx,3L,drop=TRUE])
       whichloc <- which(!is.finite(intensity2))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite upper bound of INTENS")
       } # end if.
       intensity3 <- as.numeric(pars[indx,4L,drop=TRUE])
       whichloc <- which(!is.finite(intensity3))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite initial of INTENS")
       } # end if.
       ###
       ### Check bounds.
       whichloc <- which(intensity3<intensity1 | 
                         intensity3>intensity2)
       if (length(whichloc)>=1L)  {
           stop("Error: unbounded initial of INTENS")
       } # end if.
       ### 
       ### Check logical values.
       fix_intensity <- pars[indx,5L,drop=TRUE]
       if (!all(fix_intensity %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of INTENS!")
       } # end if.
       fix_intensity <- as.logical(fix_intensity)
       intensity1[fix_intensity==TRUE] <- 
       intensity3[fix_intensity==TRUE]
       intensity2[fix_intensity==TRUE] <- 
       intensity3[fix_intensity==TRUE]
       ###
       ###
       ### 2. Activation energy. 
       ### Check non-finite values.
       energy1 <- as.numeric(pars[indx+(npeak+2L),2L,drop=TRUE])
       whichloc <- which(!is.finite(energy1)) 
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite lower bound of ENERGY!")
       } # end if.
       energy2 <- as.numeric(pars[indx+(npeak+2L),3L,drop=TRUE])
       whichloc <- which(!is.finite(energy2)) 
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite upper bound of ENERGY!")
       } # end if.
       energy3 <- as.numeric(pars[indx+(npeak+2L),4L,drop=TRUE])
       whichloc <- which(!is.finite(energy3)) 
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite initial of ENERGY!")
       } # end if.
       ###
       ### Check bounds.
       whichloc <- which(energy3<energy1 | 
                         energy3>energy2)
       if (length(whichloc)>=1L)  {
           stop("Error: unbounded initial of ENERGY")
       } # end if.
       ### Check logical values.
       fix_energy <- pars[indx+(npeak+2L),5L,drop=TRUE]
       if (!all(fix_energy %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of ENERGY!")
       } # end if.
       fix_energy <- as.logical(fix_energy)
       energy1[fix_energy==TRUE] <- 
       energy3[fix_energy==TRUE]
       energy2[fix_energy==TRUE] <- 
       energy3[fix_energy==TRUE]
       ###
       ###
       ### 3. Temperature at the peak maximum.
       ### Check non-finite values.
       temperature1 <- as.numeric(pars[indx+2L*(npeak+2L),2L,drop=TRUE])
       whichloc <- which(!is.finite(temperature1))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite lower bound of TEMPER")
       } # end if.
       temperature2 <- as.numeric(pars[indx+2L*(npeak+2L),3L,drop=TRUE])
       whichloc <- which(!is.finite(temperature2))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite upper bound of TEMPER")
       } # end if.
       temperature3 <- as.numeric(pars[indx+2L*(npeak+2L),4L,drop=TRUE])
       whichloc <- which(!is.finite(temperature3))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite initial of TEMPER")
       } # end if.
       ### Check bounds.
       whichloc <- which(temperature3<temperature1 | 
                         temperature3>temperature2)
       if (length(whichloc)>=1L)  {
           stop("Error: unbounded initial of TEMPER")
       } # end if.
       ### Check logical values.
       fix_temperature <- pars[indx+2L*(npeak+2L),5L,drop=TRUE]
       if (!all(fix_temperature %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of TEMPER!")
       } # end if.
       fix_temperature <- as.logical(fix_temperature)
       temperature1[fix_temperature==TRUE] <- 
       temperature3[fix_temperature==TRUE]
       temperature2[fix_temperature==TRUE] <- 
       temperature3[fix_temperature==TRUE]
       ###
       ###
       ### 4. bValue or rValue for the the glow peak.
       ### Check non-finite values.
       label <- ifelse(model[1L]=="g", "bValue", "rValue")
       ###if (model[1L]=="g") {
       bValue1 <- as.numeric(pars[indx+3L*(npeak+2L),2L,drop=TRUE])
       whichloc <- which(!is.finite(bValue1))
       if (length(whichloc)>=1L)  {
           stop(paste("Error: non-finite lower bound of ", label, "!", sep=""))
       } # end if.
       bValue2 <- as.numeric(pars[indx+3L*(npeak+2L),3L,drop=TRUE])
       whichloc <- which(!is.finite(bValue2))
       if (length(whichloc)>=1L)  {
           stop(paste("Error: non-finite upper bound of ", label, "!", sep=""))
       } # end if.
       bValue3 <- as.numeric(pars[indx+3L*(npeak+2L),4L,drop=TRUE])
       whichloc <- which(!is.finite(bValue3))
       if (length(whichloc)>=1L)  {
           stop(paste("Error: non-finite initial of ", label, "!", sep=""))
       } # end if.
       ### Check bounds.
       whichloc <- which(bValue3<bValue1 | 
                         bValue3>bValue2)
       if (length(whichloc)>=1L)  {
           stop(paste("Error: unbounded initial of ", label, "!", sep=""))
       } # end if.
       ### Check logical values.
       fix_bValue <- pars[indx+3L*(npeak+2L),5L,drop=TRUE]
       if (!all(fix_bValue %in% c("TRUE", "FALSE", "T", "F")))  {
           stop(paste("Error: non-logical variable in the 5th column of ", label, "!", sep=""))
       } # end if.
       fix_bValue <- as.logical(fix_bValue)
       bValue1[fix_bValue==TRUE] <- 
       bValue3[fix_bValue==TRUE]
       bValue2[fix_bValue==TRUE] <- 
       bValue3[fix_bValue==TRUE]
       ###} # end if.
       ###############################################################
       ###############################################################
       ###
       nd <- length(temp)
       n2 <- 4L*npeak
       fmin <- message <- 0
       ###
       lower <- c(intensity1, energy1, temperature1, bValue1)
       upper <- c(intensity2, energy2, temperature2, bValue2)
       pars <- c(intensity3, energy3, temperature3, bValue3)
       ###
       subroutine <- ifelse(model[1L]=="lw", "tgcd", "tgcd1")
       res <- .Fortran(subroutine, as.double(temp), as.double(signal),
                       as.integer(nd), pars=as.double(pars), as.integer(n2), 
                       fmin=as.double(fmin), message=as.integer(message), 
                       as.double(lower), as.double(upper), as.integer(nstart),
                       as.double(mdt), as.double(mwt), PACKAGE="tgcd")
       if (res$message!=0) {
           stop("Error: fail in glow curve deconvolution!")
       } # end if.
       ###
       pars <- matrix(res$pars, ncol=4L)
       index <- order(pars[,3L,drop=TRUE], decreasing=FALSE)
       pars <- pars[index,,drop=FALSE]
       colnames(pars) <- if (model[1L]=="lw") {
           c("INTENS", "ENERGY", "TEMPER", "rValue")
       } else {
           c("INTENS", "ENERGY", "TEMPER", "bValue")
       } # end if.
       rownames(pars) <- paste(seq(npeak),"th-Peak",sep="")
       ###
       kbz <- 8.617385e-5
       ###
       ####################################################################
       ### R interfaces for fortran subroutine calcei() and wrightOmega().
       ####################################################################
       wrightOmega <- function(Z) {
           W <- double(1L)
           res <- .Fortran("wrightOmega", as.double(Z), 
                           W=as.double(W), PACKAGE="tgcd")
           return(res$W)
       } # end function wrightOmega.
       ###
       calcEi <- function(x) {
           int <- as.integer(1L)
           val <- double(1L)
           res <- .Fortran("calcei", as.double(x), val=as.double(val),
                           as.integer(int), PACKAGE="tgcd")
           return(res$val)
       } # end function calcEi.
       ####################################################################
       ###
       ### Calculate signal values for each components.
       CompSig <- matrix(nrow=nd, ncol=npeak)
       if (model[1L]=="lw")  {
           for(i in seq(npeak)) {
               maxi <- pars[i,1L]
               engy <- pars[i,2L]
               maxt <- pars[i,3L]
               rv <- pars[i, 4L]
               ###
               ftev <- wz1v <- vector(length=nd)
               ###
               xi <- min(temp)
               ###
               eivi <- calcEi(-engy/kbz/xi)
               Feivi <- xi*exp(-engy/kbz/xi) + engy/kbz*eivi
               ###
               ### Calculate part1: vector wz1v.
               for (j in seq(nd)) {
                   eiv <- calcEi(-engy/kbz/temp[j])
                   ftev[j] <- (temp[j]*exp(-engy/kbz/temp[j]) + engy/kbz*eiv) - Feivi
               } # end for.
               z1v <- rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/
                      kbz/maxt^2/(1.0-1.05*rv^1.26)*ftev
               for (j in seq(nd)) wz1v[j] <- wrightOmega(z1v[j])
               ###
               ### Calculate part2: scalar wz1m.
               eiv <- calcEi(-engy/kbz/maxt)
               ftem <- (maxt*exp(-engy/kbz/maxt) +engy/kbz*eiv) - Feivi
               z1m <- rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/
                      kbz/maxt^2/(1.0-1.05*rv^1.26)*ftem
               wz1m <- wrightOmega(z1m)
               ### 
               ### Calculate signals.
               CompSig[,i] <- maxi*(wz1m+wz1m^2)/(wz1v+wz1v^2)*
                              exp(-engy/kbz*(1.0/temp-1.0/maxt))
           } # end for. 
       } else if (model[1L]=="g") {
           for(i in seq(npeak)) {
               maxi <- pars[i,1L]
               engy <- pars[i,2L]
               maxt <- pars[i,3L]
               bv <- pars[i,4L]
               xa <- 2.0*kbz*temp/engy
               xb <- 2.0*kbz*maxt/engy 
               expv <- exp(engy/kbz/temp*(temp-maxt)/maxt)
               CompSig[,i] <- maxi*(bv^(bv/(bv-1.0)))*expv*
                 ((bv-1.0)*(1.0-xa)*((temp/maxt)^2L)*expv+
                   1.0+(bv-1.0)*xb)^(-bv/(bv-1.0))
           } # end for. 
       } # end if.
       rowsumSig <- rowSums(CompSig)
       FOM <- res$fmin/sum(rowsumSig)*100
       ###
       ###
       ### Calculate shape parameters for glow peaks.
       calShape <- function(y, x)  {
           ny <- length(y)
           maxloc <- which.max(y)
           hmaxval <- max(y)/2.0
           Tm <- x[maxloc]
           T1 <- try(approx(x=y[1L:maxloc], y=x[1L:maxloc], 
                     xout=hmaxval)$y, silent=TRUE)
           T2 <- try(approx(x=y[maxloc:ny], y=x[maxloc:ny], 
                     xout=hmaxval)$y, silent=TRUE)
           ###
           if (class(T1)!="try-error") {
               d1 <- Tm-T1
           } else {
               T1 <- d1 <- NA
           } # end if.
           ###
           if (class(T2)!="try-error") {
               d2 <- T2-Tm
           } else {
               T2 <- d2 <- NA
           } # end if.  
           ###          
           thw <- T2-T1
           sf <- d2/thw
           ###         
           return( c("T1"=T1, "T2"=T2, "Tm"=Tm, "d1"=d1, 
                     "d2"=d2, "thw"=thw, "sf"=sf) )
       } # end function calShape.
       ###
       sp <- t(apply(CompSig, MARGIN=2L, calShape, temp))
       rownames(sp) <- paste(seq(npeak),"th-Peak",sep="")
       ###
       ###
       ### Calculate frequency factors for glow peaks.
       if (!is.null(hr))  {
           calff1 <- function(et)  {
               energy <- et[1L]
               temper <- et[2L]
               rv <- et[3L]
               ###
               xi <- min(temp)
               eivi <- calcEi(-energy/kbz/xi)
               Feivi <- xi*exp(-energy/kbz/xi) + energy/kbz*eivi
               eiv <- calcEi(-energy/kbz/temper)
               ftem <- (temper*exp(-energy/kbz/temper) + energy/kbz*eiv) - Feivi
               z1m <- rv/(1.0-rv) - log((1.0-rv)/rv) + energy*exp(energy/kbz/temper)/
                      kbz/temper^2/(1.0-1.05*rv^1.26)*ftem
               wz1m <- wrightOmega(z1m)
               ###
               return(hr*energy/kbz/temper^2L/exp(-energy/kbz/temper)/
                      (1.0/(1.0-rv)*(1.0+2.0*wz1m)/(1.0+wz1m)^2))             
           } # end function calff1.
           ###
           calff2 <- function(et) {
               energy <- et[1L]
               temper <- et[2L]
               bv <- et[3L]
               return(hr*energy/kbz/temper^2L*exp(energy/kbz/temper)/
                      (1.0+(bv-1.0)*(2.0*kbz*temper)/energy))
           } # end function calff2.
           ###
           if (model[1L]=="lw")  {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff1)
           } else if (model[1L]=="g") {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff2)
           } # end if.
           ###
           output <-list("pars"=pars, "ff"=ff, "sp"=sp, "FOM"=FOM)
       } else {
           output <-list("pars"=pars, "sp"=sp, "FOM"=FOM)
       } # end if
       ###
       ###
       ### Plot the results.
       if (plot==TRUE) {
           layout(cbind(c(rep(1,13), 2, rep(3,6)),
                        c(rep(1,13), 2, rep(3,6))))
           ###
           ### The first plot.
           par(mar=c(0,5.1,3.1,1.1))
           lineCol <- c("deepskyblue", "orangered", "purple", 
                        "violetred", "yellowgreen", "lightblue", 
                        "goldenrod", "forestgreen", "blue", 
                        "plum", "tan", "violet", "grey50")
           plot(temp, signal, type="p", pch=21, bg="white", cex=1.0,
                ylab="TL intensity (counts)", las=0, lab=c(7,7,9), xaxt="n", 
                xaxs="r", yaxs="i", cex.lab=2.0, cex.axis=1.5)
           box(lwd=2L)
           XaxisCentral <- median(axTicks(side=1L))
           for (i in seq(npeak)) {
               points(temp,CompSig[,i,drop=TRUE], type="l", 
                      lwd=2, col=lineCol[i])
           } # end for.
           points(temp, rowsumSig, 
                  type="l", lwd=2, col="black")
           legend(ifelse(temp[which.max(signal)]>XaxisCentral,"topleft","topright"),
                  legend=c("Fitted.Curve", paste(seq(npeak),"th-Peak",sep="")), 
                  col=c("black", lineCol[seq(npeak)]), pch=c(21, rep(NA,npeak)),
                  lty=rep("solid",npeak), yjust=2, ncol=1, cex=2.0, bty="o", 
                  lwd=2, pt.bg="white")
           ###
           ### The second plot.
           par(mar=c(0,5.1,0,1.1))
           plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
           ###
           ### The third plot.
           par(mar=c(5.1,5.1,0,1.1))
           plot(temp, signal-rowsumSig, type="p", 
                xlab="Temperature (K)", ylab="Residuals",
                las=0, lab=c(7,7,9), xaxs="r", yaxs="i", 
                pch=21, bg="grey", cex=0.8, cex.lab=2, 
                cex.axis=1.5)
           ### abline(h=0)
           box(lwd=2L)
           ###
           par(mar=c(5,4,4,2)+0.1)
           layout(1L)
       } # end if.
       ###
       if (!is.null(outfile)) {
           CompSig <- cbind(temp, signal, rowsumSig,  CompSig)
           colnames(CompSig) <- c("Temperature", "Obs.Signal", 
               "Fit.Signal", paste("Comp.", seq(npeak), sep = ""))
           write.csv(CompSig, file=paste(outfile, ".csv", sep = ""))
       } # end if.
       ###
       return(output)                 
} # end fucntion tgcd.
###

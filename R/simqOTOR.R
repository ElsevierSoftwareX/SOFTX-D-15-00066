###
simqOTOR <- 
function(temps, n0, Nn, Ah, An, ff, ae,
         hr, outfile=NULL, plot=TRUE) {
    UseMethod("simqOTOR")
} #
### 2016.05.22.
simqOTOR.default <- 
function(temps, n0, Nn, Ah, An, ff, ae,
         hr, outfile=NULL, plot=TRUE) {
    ### Stop if not.
    stopifnot(is.numeric(temps), all(temps>0),
              length(n0)==1L, n0>0,
              length(Nn)==1L, is.numeric(Nn), Nn>0, Nn>=n0,
              length(Ah)==1L, is.numeric(Ah), Ah>0,
              length(An)==1L, is.numeric(An), An>0,
              length(ff)==1L, is.numeric(ff), ff>0,
              length(ae)==1L, is.numeric(ae), ae>0,
              length(hr)==1L, is.numeric(hr), hr>0,
              is.null(outfile) || is.character(outfile),
              length(plot)==1L, is.logical(plot))
    ###
    if (!is.null(outfile)) {
        if (length(outfile)!=1L) 
            stop("Error: outfile should be an one-element vector!")
    } # end if.
    ###
    nt <- length(temps)
    vecy <- double(nt)
    info <- -100
    ###
    res <- .Fortran("qeOTOR", as.integer(nt), as.double(temps),
                    as.double(n0), as.double(Nn), as.double(Ah), 
                    as.double(An), as.double(ff), as.double(ae),
                    as.double(hr), vecy=as.double(vecy), 
                    info=as.integer(info), PACKAGE="tgcd")
    if(res$info!=2L) {
        stop("Error: fail in solving ordinary equation!")
    } # end if.
    ###
    kbz = 8.617385e-05
    tout <- temps
    yout <- ff*(res$vecy)^2L*exp(-ae/kbz/temps)*Ah/
            (Ah*res$vecy+An*(Nn-res$vecy))/hr
    ###
    calShape <- function(x, y)  {
        ny <- length(y)
        maxloc <- which.max(y)
        hmaxval <- max(y)/2.0
        Tm <- x[maxloc]
        T1 <- approx(x=y[1L:maxloc], y=x[1L:maxloc], xout=hmaxval)$y
        T2 <- approx(x=y[maxloc:ny], y=x[maxloc:ny], xout=hmaxval)$y
        d1 <- Tm-T1
        d2 <- T2-Tm
        thw <- T2-T1
        sf <- d2/thw
        return( c("T1"=T1, "T2"=T2, "Tm"=Tm, "d1"=d1, 
                  "d2"=d2, "thw"=thw, "sf"=sf) )
    } # end function calShape.
    ###
    sp <- try(calShape(tout, yout), silent=TRUE)
    ###
    if (plot==TRUE)  {
        layout(cbind(c(rep(1,13), 2, rep(3,6)), 
                     c(rep(1,13), 2, rep(3,6))))
        ###
        ### The first plot.
        par(mar=c(0,5.1,3.1,1.1))
        plot(tout, yout, type="l", lwd=5, col="skyblue3",
             ylab="TL Intensity (counts)", las=0, lab=c(7,7,9), 
             xaxt="n", xaxs="r", yaxs="r", cex.lab=2.0, cex.axis=1.5)
        box(lwd=2L)
        ###
        XaxisCentral <- median(axTicks(side=1L))
        ###
        if (class(sp)=="try-error")  {
            legend(ifelse(tout[which.max(yout)] > XaxisCentral, "topleft", 
                   "topright"), legend=c(paste("n0: ",
                   format(n0,digits=3,scientific=TRUE)," (1/cm^3)",sep=""),
                   paste("ae: ",round(ae,2L)," (eV)",sep=""),
                   paste("ff: ",format(ff,digits=3,scientific=TRUE)," (1/s)",sep=""), 
                   paste("hr: ",round(hr,2L)," (K/s)",sep=""), 
                   paste("Ah: ",format(Ah,digits=3,scientific=TRUE)," (cm^3/s)",sep=""),
                   paste("An: ",format(An,digits=3,scientific=TRUE)," (cm^3/s)",sep="")),
                   yjust=2, ncol=1, cex=2.0, bty="n", pt.bg="white")
        } else {
            legend(ifelse(tout[which.max(yout)] > XaxisCentral, "topleft", 
                   "topright"), legend=c(paste("n0: ",
                   format(n0,digits=3,scientific=TRUE)," (1/cm^3)",sep=""),
                   paste("ae: ",round(ae,2L)," (eV)",sep=""),
                   paste("ff: ",format(ff,digits=3,scientific=TRUE)," (1/s)",sep=""), 
                   paste("hr: ",round(hr,2L)," (K/s)",sep=""), 
                   paste("Ah: ",format(Ah,digits=3,scientific=TRUE)," (cm^3/s)",sep=""),
                   paste("An: ",format(An,digits=3,scientific=TRUE)," (cm^3/s)",sep=""),
                   paste("sf: ",round(sp["sf"],2L),sep="")),
                   yjust=2, ncol=1, cex=2.0, bty="n", pt.bg="white")
        } # end if.
        ###
        ### The second plot.
        par(mar=c(0,5.1,0,1.1))
        plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="n", ylab="")
        ###
        ### The third plot.
        par(mar = c(5.1, 5.1, 0, 1.1))
        plot(tout, res$vecy, type="l", xlab = "Temperature(K)", 
             ylab="Concentration", las=0, lab=c(7,7,9), xaxs="r", yaxs="i",
             ylim=c(-0.1*n0, 1.1*n0), col="grey", lwd=5, cex.lab=2.0, cex.axis=1.5)
        box(lwd=2L) 
        par(mar=c(5,4,4,2)+0.1)
        layout(1L)
    } # end if.
    ###
    if (class(sp)=="try-error")  sp <- NULL
    ###
    if (!is.null(outfile))  {
        PeakSig <- cbind(tout, yout, res$vecy)
        colnames(PeakSig) <- c("temps","tl", "n")
        write.csv(PeakSig, file=paste(outfile, ".csv", sep=""))
    } # end if.
    return(invisible(list("temps"=tout, "tl"=yout, "n"=res$vecy, "sp"=sp)))
} # end function simqOTOR.
###

subroutine qeOTOR(nt, vect, y0, Nn, Ah,& 
                  An, ff, ae, hr, vecy, info)
!--------------------------------------------------------------------------------
! Subroutine qeOTOR is used for solving QE approximation equation 
! for the one-trap-one recombination center (OTOR) model.
!--------------------------------------------------------------------------------
!        nt, input:: integer, the number of data points.
!  vect(nt), input:: real values, temperature values.
!        y0, input:: real value, starting value of n.
!        Nn, input:: total concentration of traps in the crystal.
!        Ah, input:: real value, probability of electron recombine with holes.
!        An, input:: real value, probability of electron retrapping.
!        ff, input:: real value, the frequency factor.
!        ae, input:: real value, the activation energy.
!        hr, input:: real value, the heating rate.
! vecy(nt), output:: real values, y values at various temperatures.
!     info, output:: integer, error message, info=2 means a successful work.
!--------------------------------------------------------------------------------
! Author:: Peng Jun, 2015.09.13; revised in 2016.01.19.
!--------------------------------------------------------------------------------
! Dependence:: subroutine dlsoda.
!--------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER nt, info
    DOUBLE PRECISION vect(nt), vecy(nt), y0,& 
                     Nn, Ah, An, ff, ae, hr
    !
    INTEGER NEQ(1), ITOL, ITASK, ISTATE, IOPT, & 
            JT, LRW, LIW, IWORK(21), i
    DOUBLE PRECISION RWORK(36), Y(1), T, TOUT, & 
                     RTOL(1), ATOL(1), kbz, bv
    EXTERNAL FUN, JAC
    !
    kbz = 8.617385d-05
    !  
    NEQ = 1
    ITOL = 1
    ITASK = 1
    ISTATE = 1
    IOPT = 1
    LRW = 36
    LIW = 21
    JT = 2

    RWORK = 0.0d+00
    RWORK(6) = maxval(vect(2:nt)-vect(1:nt-1))

    IWORK = 0
    IWORK(1) = 1
    IWORK(2) = 1
    IWORK(6) = 50000
    !
    vecy(1) = y0
    RTOL = 1.0d-6 
    ATOL = 1.0d-6
    !
    bv = 0.0D+00
    !
    DO i=1, nt-1
        Y = vecy(i)
        T = vect(i)
        TOUT = vect(i+1)
        !
        CALL DLSODA (FUN, NEQ, Y, T, TOUT, ITOL, RTOL,& 
                     ATOL, ITASK, ISTATE, IOPT, RWORK,&
                     LRW, IWORK, LIW, JAC, JT, ff, ae,& 
                     Ah, An, Nn, hr, bv)
        info = ISTATE
        IF (ISTATE .LT. 0) RETURN   
        vecy(i+1) = Y(1)
    END DO
    !
    RETURN
END SUBROUTINE qeOTOR
!
SUBROUTINE FUN(NEQ, T, Y, YDOT, & 
               ff, ae, Ah, An, Nn, hr, bv)
    IMPLICIT NONE
    INTEGER NEQ
    DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), & 
                     kbz, ff, ae, Ah, An, Nn, hr, bv
    kbz = 8.617385d-05 
    YDOT(1) = -ff*(Y(1))**2*dexp(-ae/kbz/T)*&
               Ah/(Ah*Y(1)+(Nn-Y(1))*An)/hr + bv*bv        
    RETURN
END SUBROUTINE FUN
!
SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
    IMPLICIT NONE
    INTEGER NEQ, ML, MU, NROWPD
    DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ)
    RETURN
END SUBROUTINE JAC

subroutine simPeak(nt, vect, y0, Nn, bv, model,&
                   ff, ae, hr, vecy, info)
!--------------------------------------------------------------------------------
! Subroutine qeOTOR is used for solving QE approximation equation 
! for the one-trap-one recombination center (OTOR) model.
!--------------------------------------------------------------------------------
!        nt, input:: integer, the number of data points.
!  vect(nt), input:: real values, temperature values.
!        y0, input:: real value, starting value of n.
!        Nn, input:: total concentration of traps in the crystal.
!     model, input:: integer, 1=first-order, 2=second-order.
!        ff, input:: real value, the frequency factor.
!        ae, input:: real value, the activation energy.
!        hr, input:: real value, the heating rate.
! vecy(nt), output:: real values, y values at various temperatures.
!     info, output:: integer, error message, info=2 means a successful work.
!--------------------------------------------------------------------------------
! Author:: Peng Jun, 2015.09.14; revised in 2016.01.19.
!--------------------------------------------------------------------------------
! Dependence:: subroutine dlsoda.
!--------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER nt, model, info
    DOUBLE PRECISION vect(nt), vecy(nt), y0,& 
                     Nn, ff, ae, hr, bv
    !
    INTEGER NEQ(1), ITOL, ITASK, ISTATE, IOPT, & 
            JT, LRW, LIW, IWORK(21), i
    DOUBLE PRECISION RWORK(36), Y(1), T, TOUT, & 
                     RTOL(1), ATOL(1), kbz, Ah, An
    EXTERNAL FUN1, FUN2, FUN3, JAC
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
    RTOL = 1.0d-06
    ATOL = 1.0d-06
    !
    DO i=1, nt-1
        Y = vecy(i)
        T = vect(i)
        TOUT = vect(i+1)
        !
        IF (model==1)  then
            Ah = 0.0d+00
            An = 0.0d+00
            Nn = 0.0d+00
            bv = 0.0d+00
            CALL DLSODA (FUN1, NEQ, Y, T, TOUT, ITOL, RTOL,& 
                         ATOL, ITASK, ISTATE, IOPT, RWORK,&
                         LRW, IWORK, LIW, JAC, JT, ff, ae,& 
                         Ah, An, Nn, hr, bv)
        ELSE IF (model==2) then
            Ah = 0.0d+00
            An = 0.0d+00
            bv = 0.0d+00
            CALL DLSODA (FUN2, NEQ, Y, T, TOUT, ITOL, RTOL,& 
                         ATOL, ITASK, ISTATE, IOPT, RWORK,&
                         LRW, IWORK, LIW, JAC, JT, ff, ae,& 
                         Ah, An, Nn, hr, bv)
        ELSE IF (model==3) then
            Ah = 0.0d+00
            An = 0.0d+00
            CALL DLSODA (FUN3, NEQ, Y, T, TOUT, ITOL, RTOL,& 
                         ATOL, ITASK, ISTATE, IOPT, RWORK,&
                         LRW, IWORK, LIW, JAC, JT, ff, ae,& 
                         Ah, An, Nn, hr, bv)
        END IF
        info = ISTATE
        IF (ISTATE .LT. 0) RETURN   
        vecy(i+1) = Y(1)
    END DO
    !
    RETURN
END SUBROUTINE simPeak
!
SUBROUTINE FUN1(NEQ, T, Y, YDOT, & 
                ff, ae, Ah, An, Nn, hr, bv)
    IMPLICIT NONE
    INTEGER NEQ
    DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), & 
                     kbz, ff, ae, hr
    DOUBLE PRECISION Ah, An, Nn, bv
    kbz = 8.617385d-05 
    YDOT(1) = -ff*Y(1)*dexp(-ae/kbz/T)/hr + (Ah*An*Nn*bv)      
    RETURN
END SUBROUTINE FUN1
!
SUBROUTINE FUN2(NEQ, T, Y, YDOT, & 
                ff, ae, Ah, An, Nn, hr, bv)
    IMPLICIT NONE
    INTEGER NEQ
    DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), & 
                     kbz, ff, ae, hr, Nn
    DOUBLE PRECISION Ah, An, bv
    kbz = 8.617385d-05 
    YDOT(1) = -ff*(Y(1))**2*dexp(-ae/kbz/T)/Nn/hr + (Ah*An*bv)      
    RETURN
END SUBROUTINE FUN2
!
SUBROUTINE FUN3(NEQ, T, Y, YDOT, & 
                ff, ae, Ah, An, Nn, hr, bv)
    IMPLICIT NONE
    INTEGER NEQ
    DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), & 
                     kbz, ff, ae, hr, Nn, bv
    DOUBLE PRECISION Ah, An
    kbz = 8.617385d-05 
    YDOT(1) = -ff*(Y(1))**bv*dexp(-ae/kbz/T)/Nn/hr + (Ah*An)      
    RETURN
END SUBROUTINE FUN3

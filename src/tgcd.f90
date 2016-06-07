subroutine tgcd(xd,yd,nd,pars,n2,fmin,message,&
                lower,upper,nstart,mdt,mwt)
!--------------------------------------------------------------
! Subroutine tgcd() is used for thermoluminescence glow curve
! deconvolution according to the Lambert’s W function using
! the Levenberg–Marquardt algorithm.
!--------------------------------------------------------------
!    xd(nd):: input, real values, observation X.
!    yd(nd):: input, real vlaues, observations Y.
!        nd:: input, integer, number of points.
!  pars(n2):: input/output, paraneters.
!        n2:: input, integer, number of pars (<=30).
!      fmin:: output, real value, minimized objective.
!   message:: output, integer, error message:
!                     0=success, 1=fail.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real values, upper bounds.
!    nstart:: input, integer, number of trials.
!       mdt:: input, real value, minimum distance
!             between each temperature.
!       mwt:: input, real value, maximum total 
!             half-width for glow peaks.
!-------------------------------------------------------------
! Author:: Peng Jun, 2016.05.21.
!-------------------------------------------------------------
! Dependence:: subroutine lmtl; 
!              subroutine hpSort;
!              subroutine calcD;
!              subroutine calcMaty.
!-------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nd, n2, nstart
    real   (kind=8), intent(in):: xd(nd), yd(nd), mdt, mwt
    real   (kind=8), intent(in):: lower(n2), upper(n2)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    real   (kind=8):: ran(n2), ranpars(n2), ranfmin, &
                      minfmin, unifa(n2), unifb(n2),&
                      orderTemp(n2/4), mindist,&
                      maty(nd,n2/4), maxwidth
    integer(kind=4):: i, info, indx(n2/4), flag
    !
    fmin = -99.0
    message = 1
    !
    ! Do optimization for the first time.
    ranpars = pars
    call lmtl(xd,yd,nd,ranpars,n2,ranfmin,&
              info,lower,upper)
    !
    if (info==0) then
        ! Calculate minimum distance between peaks.
        orderTemp = ranpars(2*n2/4+1:3*n2/4)
        call hpSort(orderTemp, n2/4, indx) 
        mindist = minval(orderTemp(2:n2/4)-&
                         orderTemp(1:n2/4-1))
        !
        ! Calculate maximum total half-width of peaks.
        call calcMaty(nd,n2,ranpars,xd,maty)
        call calcD(nd,n2/4,xd,maty,maxwidth,flag)
        !
        if (mindist>=mdt .and. flag==0 .and. maxwidth<=mwt) then
            message = 0
            pars = ranpars
            fmin = ranfmin
        end if
    end if
    !
    if (nstart==1)  return
    !
    if (message/=0)  then
        minfmin = 1.0e+20
    else 
        minfmin = fmin          
    end if
    !
    ! INTENS.
    unifa(1:n2/4) = pars(1:n2/4)*0.98
    unifb(1:n2/4) = pars(1:n2/4)*1.02
    ! ENERGY.
    unifa(n2/4+1:2*n2/4) = 1.8
    unifb(n2/4+1:2*n2/4) = 2.2
    ! TEMPER.
    unifa(2*n2/4+1:3*n2/4) = pars(2*n2/4+1:3*n2/4)*0.98
    unifb(2*n2/4+1:3*n2/4) = pars(2*n2/4+1:3*n2/4)*1.02
    ! rValue. 
    unifa(3*n2/4+1:n2) = 0.1
    unifb(3*n2/4+1:n2) = 0.9
    !
    ! Try-and-error.
    do i=1,  nstart
        call random_number(ran)
        ranpars = unifa + ran*(unifb-unifa)
        !
        call lmtl(xd,yd,nd,ranpars,n2,ranfmin,&
                  info,lower,upper)
        if (info==0) then
            ! Calculate minimum distance between peaks.
            orderTemp = ranpars(2*n2/4+1:3*n2/4)
            call hpSort(orderTemp, n2/4, indx) 
            mindist = minval(orderTemp(2:n2/4)-orderTemp(1:n2/4-1))
            !
            ! Calculate maximum total half-width of peaks.
            call calcMaty(nd,n2,ranpars,xd,maty)
            call calcD(nd,n2/4,xd,maty,maxwidth,flag)
            !!!print*, maxwidth, ranfmin
            !
            if (ranfmin<minfmin .and. mindist>=mdt .and. &
                flag==0 .and. maxwidth<=mwt)  then
                pars = ranpars
                fmin = ranfmin
                minfmin = ranfmin
                message = 0
            end if
        end if
    end do
    !
    return
end subroutine tgcd
!
subroutine lmtl(xd,yd,nd,pars,n2,fmin,&
                message,lower,upper)
!----------------------------------------------------
! Subroutine lmtl() is used for fitting a TL glow 
! curve using the Levenberg-Marquardt algorithm.
!----------------------------------------------------
!    xd(nd):: input, real values, observation X.
!    yd(nd):: input, real vlaues, observations Y.
!        nd:: input, integer, number of points.
!  pars(n2):: input/output, paraneters.
!        n2:: input, integer, number of pars (<=30).
!      fmin:: output, real value, minimized objective.
!   message:: output, integer, error message:
!                     0=success, 1=fail.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real values, upper bounds.
!------------------------------------------------------
! Author:: Peng Jun, 2016.05.20.
!------------------------------------------------------
! Dependence:: subroutine lmdif1; 
!              subroutine tgcfunc.
!------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nd, n2
    real   (kind=8), intent(in):: xd(nd), yd(nd)
    real   (kind=8), intent(in):: lower(n2), upper(n2)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: fmin
    integer(kind=4), intent(out):: message
    !
    integer(kind=4):: info
    real   (kind=8), parameter:: tol=1.0e-07
    real   (kind=8):: fvec(nd)
    external:: tgcfunc
    !
    fmin = -99.0
    !
    call lmdif1(tgcfunc,nd,n2,pars,fvec,tol,&
                info,xd,yd,lower,upper)
    !
    if (info==1 .or. info==2 .or. info==3) then
        message = 0
    else 
        message = 1
        return
    end if
    !
    fmin = sum(fvec**2)
    !
    return
end subroutine lmtl
!
subroutine tgcfunc(nd,n2,pars,fvec,iflag,&
                   xd,yd,lower,upper)
!---------------------------------------------------
! Subroutine tgcfunc() is used for calculating
! the residual vector of a TL glow curve.
!---------------------------------------------------
!        nd:: input, integer, number of data points.
!        n2:: input, integer, number of pars (<=52).
!  pars(n2):: input, real values, pars.
!  fvec(nd):: output, real values, residuals.
!     iflag:: input, integer, working variable.
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real vlaues, upper bounds.
!----------------------------------------------------
! Author:: Peng Jun, 2016.05.20.
!----------------------------------------------------
! Dependence:: subroutine calcei; 
!              subroutine wrightOmega.
!----------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)               
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(52), maxi, engy, maxt, &
                      rv, eiv, eivi, Feivi, ftev(nd), ftem, &
                      z1v(nd), z1m, wz1v(nd), wz1m, xi, wv
    integer(kind=4):: i, j
    !
    ! Bound constraints.
    do i=1, n2
        if (pars(i)<lower(i))  then
            pars(i) = lower(i)
        else if (pars(i)>upper(i)) then
            pars(i) = upper(i)
        end if
    end do
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    xi = minval(xd)
    !
    fvec = 0.0
    do i=1, n2/4
        maxi = xx(i)
        engy = xx(i+n2/4)
        maxt = xx(i+2*n2/4)
        rv = xx(i+3*n2/4)
        !
        call calcei(-engy/kbz/xi, eivi, 1)
        Feivi = xi*exp(-engy/kbz/xi) + engy/kbz*eivi
        !
        ! Calculate part1: vector wz1v.
        do j=1, nd
            call calcei(-engy/kbz/xd(j), eiv, 1)
            ftev(j) = (xd(j)*exp(-engy/kbz/xd(j)) + engy/kbz*eiv) - Feivi
        end do
        !
        z1v = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
              kbz/maxt**2/(1.0-1.05*rv**1.26)*ftev
        !
        do j=1, nd
            call wrightOmega(z1v(j), wv)
            wz1v(j) = wv
        end do 
        !
        !
        ! Calculate part2: scalar wz1m.
        call calcei(-engy/kbz/maxt, eiv, 1)
        ftem = (maxt*exp(-engy/kbz/maxt) + engy/kbz*eiv) - Feivi
        !
        z1m = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
              kbz/maxt**2/(1.0-1.05*rv**1.26)*ftem
        !
        call wrightOmega(z1m, wv)
        wz1m = wv
        !
        ! Calculate residuals.
        fvec = fvec + maxi*(wz1m+wz1m**2)/(wz1v+wz1v**2)*&
               exp(-engy/kbz*(1.0/xd-1.0/maxt))
        !
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc
!
subroutine calcMaty(nd,n2,pars,xd,maty)
!-----------------------------------------------------------
! Subroutine calcMaty() is used for calculating
! the signal matrix of optimized TL glow peaks.
!----------------------------------------------------------
!            nd:: input, integer, number of data points.
!            n2:: input, integer, number of pars (<=52).
!      pars(n2):: input, real values, pars.
!        xd(nd):: input, real values, observations X.
! maty(nd,n2/4):: output, real values, calculated signals.
!----------------------------------------------------------
! Author:: Peng Jun, 2016.05.21.
!----------------------------------------------------------
! Dependence:: subroutine calcei; 
!              subroutine wrightOmega.
!----------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4), intent(in):: nd, n2
    real   (kind=8), intent(in):: pars(n2),xd(nd)
    real   (kind=8), intent(out):: maty(nd,n2/4)                 
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(52), maxi, engy, maxt, &
                      rv, eiv, eivi, Feivi, ftev(nd), ftem, &
                      z1v(nd), z1m, wz1v(nd), wz1m, xi, wv
    integer(kind=4):: i, j
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    xi = minval(xd)
    !
    do i=1, n2/4
        maxi = xx(i)
        engy = xx(i+n2/4)
        maxt = xx(i+2*n2/4)
        rv = xx(i+3*n2/4)
        !
        call calcei(-engy/kbz/xi, eivi, 1)
        Feivi = xi*exp(-engy/kbz/xi) + engy/kbz*eivi
        !
        ! Calculate part1: vector wz1v.
        do j=1, nd
            call calcei(-engy/kbz/xd(j), eiv, 1)
            ftev(j) = (xd(j)*exp(-engy/kbz/xd(j)) + engy/kbz*eiv) - Feivi
        end do
        !
        z1v = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
              kbz/maxt**2/(1.0-1.05*rv**1.26)*ftev
        !
        do j=1, nd
            call wrightOmega(z1v(j), wv)
            wz1v(j) = wv
        end do 
        !
        !
        ! Calculate part2: scalar wz1m.
        call calcei(-engy/kbz/maxt, eiv, 1)
        ftem = (maxt*exp(-engy/kbz/maxt) + engy/kbz*eiv) - Feivi
        !
        z1m = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
              kbz/maxt**2/(1.0-1.05*rv**1.26)*ftem
        !
        call wrightOmega(z1m, wv)
        wz1m = wv
        !
        ! Calculate residuals.
        maty(:,i) = maxi*(wz1m+wz1m**2)/(wz1v+wz1v**2)*&
                    exp(-engy/kbz*(1.0/xd-1.0/maxt))
        !
    end do
    !
    return
    !
end subroutine calcMaty

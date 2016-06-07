subroutine calcD(nd,np,xd,maty,maxd,flag)
!------------------------------------------------------------------------
! Subroutine calcD() is used for calculating the
! maximum total half-width from a series of glow peaks.
!------------------------------------------------------------------------
!          nd:: integer, input, number of data points.
!          np:: integer, input, number of glow peaks.
!      xd(nd):: real values, input, temperature values.
! maty(nd,np):: real values, input, TL values for each component.
!        maxd:: real value, output, calculated maximum total half-width.
!        flag:: integer, output, error flag, 0=success, 1=fail.
!------------------------------------------------------------------------
! Author:: Peng jun, 2016.05.21.
!------------------------------------------------------------------------
! Dependence:: No.
!------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nd, np
    real   (kind=8), intent(in):: xd(nd), maty(nd,np)
    real   (kind=8), intent(out):: maxd
    integer(kind=4), intent(out):: flag
    !
    integer(kind=4):: i, j, indx, loc1
    real   (kind=8):: yd(nd), hymax, T1, T2,& 
                      x0, y0, x1, y1, width
    !  
    maxd=-1.0e20
    flag = 0
    !
    do i=1, np
        yd=maty(:,i)
        hymax=maxval(yd)/2.0
        indx=maxloc(yd,dim=1)
        !
        ! Calculate T1 for each component.
        loc1 = -99
        do j=1, indx-1
            if (hymax>yd(j) .and. hymax<yd(j+1)) then
                loc1=j
            end if 
        end do
        !
        if (loc1==-99) then
            flag = 1
            return
        end if
        !
        x0=yd(loc1)
        y0=xd(loc1)
        x1=yd(loc1+1)
        y1=xd(loc1+1)
        T1=y0*(hymax-x1)/(x0-x1) + y1*(hymax-x0)/(x1-x0)
        !
        ! Calculate T2 for each component.
        loc1 = -99
        do j=indx, nd-1
            if (hymax<yd(j) .and. hymax>yd(j+1)) then
                loc1=j
            end if 
        end do
        !
        if (loc1==-99) then
            flag = 1
            return
        end if
        !
        x0=yd(loc1)
        y0=xd(loc1)
        x1=yd(loc1+1)
        y1=xd(loc1+1)
        T2=y0*(hymax-x1)/(x0-x1) + y1*(hymax-x0)/(x1-x0)
        !
        ! Calculate width for each component.
        width=T2-T1
        ! 
        ! Select a maximum width.
        if(width>maxd) maxd=width
    end do
    !
    return
    !
end subroutine calcD

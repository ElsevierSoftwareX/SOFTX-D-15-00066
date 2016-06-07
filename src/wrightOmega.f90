subroutine wrightOmega(Z, W)
!----------------------------------------------------------
! Subroutine wrightOmega is used for calculating 
! the wright Omega function for real values.
!----------------------------------------------------------
! Z:: input, real value, a value from which the
!     wright Omega function will be calculated.
! W:: output, real value, the calculated 
!     wright Omega function for Z value.
!----------------------------------------------------------
! Author:: Peng Jun, 2016.05.19.
!----------------------------------------------------------
! Dependence:: No.
!--------------------------------------------------------------------------------
! References: Piers WL, Robert MC, David JJ, Algorithm 917: Complex 
!             Double-Precision Evaluation of the Wright omega Function.
!             ACM Transactions on Mathematical Software, 38(3):1-17.
!-------------------------------------------------------------------------------
!      Note: This subroutine is transformed from the matlab code provided by 
!            Andrew Horchler available at:
!            https://github.com/horchler/wrightOmegaq/blob/master/wrightOmegaq.m
!-------------------------------------------------------------------------------
    implicit none
    real   (kind=8), intent(in):: Z
    real   (kind=8), intent(out):: W
    ! Local variables.
    real   (kind=8):: X, y, lzi, r, w1, w2
    real   (kind=8), parameter:: PI=3.141593D+00,&
                                 PI2=1.570796D+00,&
                                 EXP1=2.718282D+00,&
                                 EXP2=7.389056D+00,&
                                 LN2=0.6931472D+00,&
                                 omega=0.5671432904097838D+00,&
                                 oneThird=0.3333333333333D+00,&
                                 tol=2.2204D-16,&
                                 lg=5.764608D+17
    X=Z
    !
    if(Z>lg) then
        W=Z
    else if(Z==0.0D+00) then
        W=omega
    else if(Z==1.0D+00) then
        W=1.0D+00
    else if(Z==1.0D+00+EXP1) then
        W=EXP1
    else 
        if(Z<-745.1332D+00) then
            W=0.0D+00
        else 
            if(Z<=-2.0D+00) then
                y=exp(X)
                W=y*(1.0D+00-y*(1.0D+00-y*(36.0D+00-&
                  y*(64.0D+00-125.0D+00*y))/24.0D+00))
                if(X<-EXP2) return
            else if(Z>PI+1.0D+00) then
                y=log(X)
                lzi=y/X
                W=X-y+lzi*(1.0D+00+lzi*(0.5D+00*y-1.0D+00+&
                  lzi*((oneThird*y-1.5D+00)*y+1.0D+00)))
            else 
                y=X-1.0D+00
                W=1.0D+00+y*(0.5D+00+y*(1.0D+00/16.0D+00-y*(1.0D+00/192.0D+00+&
                  y*(1.0D+00/3072.0D+00-13.0D+00/61440.0D+00*y))))
            end if
            !
            r=X-(W+log(W))
            if(abs(r)>tol) then
                w1=1.0D+00+W
                w2=w1+2.0D+00*oneThird*r
                W=W*(1.0D+00+r*(w1*w2-0.5D+00*r)/(w1*(w1*w2-r)))
                r=X-(W+log(W))
                !
                if(abs(r)>tol) then
                    w1=1.0D+00+W
                    w2=w1+2.0D+00*oneThird*r
                    W=W*(1.0D+00+r*(w1*w2-0.5D+00*r)/(w1*(w1*w2-r)))
                end if
            end if
            !
        end if
        !
    end if
    !
    return
end subroutine wrightOmega

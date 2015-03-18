



  subroutine eigenbound_geometry2
! include modules

  use param
  use arrays
  use derivatives
  use procinfo
  
  implicit none

  integer k
  logical contains

  real(8) dr_Aa,idr,dr_phi,dr_B,dr_K,B

  idr = 1.0d0/dr
  k=Nr

! ************************
! ***   SANITY CHECK   ***
! ************************

  if (rank==0) then

!    Only allow standard BSSN.

     if (eta/=2.d0) then
        print *
        print *, 'Boundary=eigen only works for eta=2.'
        print *
        call die
     end if

  end if
! **********************************************
! ***   DERIVATIVES OF SOURCES AT BOUNDARY   ***
! **********************************************

! Derivative of salpha.

  if (order=="two") then

!    I do this to fourth order to improve the boundary behavior.

     dr_salpha =0
     dr_K=0

  else if (order=="four") then

!    I do this to fifth order to improve the boundary behavior.

     dr_salpha =0
     dr_K=0

  end if

! Derivative of phi.

  if (order=="two") then

!    I do this to fourth order to improve the boundary behavior.

     dr_phi = idr*(25.d0*phi(Nr) - 48.d0*phi(Nr-1) &
             + 36.d0*phi(Nr-2) - 16.d0*phi(Nr-3) + 3.d0*phi(Nr-4))/12.d0

  else if (order=="four") then

!    I do this to fifth order to improve the boundary behavior.

     dr_phi = idr*(137.0d0/60.0d0*phi(Nr) - 5.0d0*phi(Nr-1) + 5.0d0*phi(Nr-2) &
             - 10.0d0/3.0d0*phi(Nr-3) + 5.0d0/4.0d0*phi(Nr-4) - 0.2d0*phi(Nr-5))

! Derivative of B.
  if (order=="two") then

!    I do this to fourth order to improve the boundary behavior.

     dr_B = idr*(25.d0*B(Nr) - 48.d0*B(Nr-1) &
           + 36.d0*B(Nr-2) - 16.d0*B(Nr-3) + 3.d0*B(Nr-4))/12.d0

  else if (order=="four") then

!    I do this to fifth order to improve the boundary behavior.

     dr_B = idr*(137.0d0/60.0d0*B(Nr) - 5.0d0*B(Nr-1) + 5.0d0*B(Nr-2) &
           - 10.0d0/3.0d0*B(Nr-3) + 5.0d0/4.0d0*B(Nr-4) - 0.2d0*B(Nr-5))

  end if
! Principal equation
  dr_Aa=2/3*dr_K+6*Aa*dr_phi(Nr)+2/3*Aa*(2*1/r+dr_B/B(Nr))
  if (order=='two') then
     Aa = 3.d0*dr*D1_Aa
  else if (order=='four') then
     Aa = (12.d0*dr*d1_Aa + 48.d0*auxarray(Nr-1) - 36.d0*auxarray(Nr-2) &
            + 16.d0*auxarray(Nr-3) - 3.d0*auxarray(Nr-4))/25.d0
  end
																			

  end subroutine eigenbound_geometry 

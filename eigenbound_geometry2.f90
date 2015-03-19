
! Ax

  subroutine eigenbound_geometry2
! include modules

  use param
  use arrays
  use derivatives
  use procinfo
  
  implicit none

  integer k
  logical contains

  idr = 1.0d0/dr

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

     dr_alpha =0
     dr_K=0
     alpha[Nr]=alpha[Nr]+idr*0.5*(3*alpha[Nr]-alpha[Nr-1]+alpha[Nr-2]) 
     K[Nr]=K[Nr]+idr*0.5*(3*K[Nr]-K[Nr-1]+K[Nr-2])

  else if (order=="four") then

     dr_alpha =0
     dr_K=0
     alpha[Nr]=alpha[Nr]+(idr/12)*(25*alpha[Nr]-48*alpha[Nr-1]+36*alpha[Nr-2]-16*alpha[Nr-3]+3*alpha[Nr-4])
     K[Nr]=K[Nr]+(idr/12)*(25*K[Nr]-48*K[Nr-1]+36*K[Nr-2]-16*K[Nr-3]+3*K[Nr-4])


  end if

! Derivative of phi.

  if (order=="two") then

     dr_phi = (idr/12)*(idr*0.5*(3*phi[Nr]-phi[Nr-1]+phi[Nr-2])

  else if (order=="four") then

     dr_phi = idr/12*(25*phi[Nr]-48*phi[Nr-1]+36*phi[Nr-2]-16*phi[Nr-3]+3*phi[Nr-4])

! Derivative of B.
  if (order=="two") then

     dr_B = (idr/12)*(idr*0.5*(3*B[Nr]-B[Nr-1]+B[Nr-2])

  else if (order=="four") then

     dr_B = idr/12*(25*B[Nr]-48*B[Nr-1]+36*B[Nr-2]-16*B[Nr-3]+3*B[Nr-4])

  end if
! Principal equation
  dr_KTA=2/3*dr_K+6*KTA[Nr]*dr_phi(Nr)+2/3*Aa*(2*1/r+dr_B/B(Nr))
  idr_KTA=1./dr_KTA
  if (order=="two") then

     KTA =KTA[Nr]+(idr_KTA/12)*(idr*0.5*(3*KTA[Nr]-KTA[Nr-1]+KTA[Nr-2])

  else if (order=="four") then

     KTA[] =KTA[Nr]+idr_KTA/12*(25*B[Nr]-48*KTA[Nr-1]+36*KTA[Nr-2]-16*KTA[Nr-3]+3*KTA[Nr-4])
  end
																			

  end subroutine eigenbound_geometry2

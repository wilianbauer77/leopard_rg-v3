!> Computes the spline coefficients for both parallel and perpendicular interpolation of the velocity distribution data
!! \param iarb index of species among the particle species with arbitrary velocity distribution which is to be interpolated
!! \param splcoeff1c array of spline coefficients for which distribution was first interpolated in parallel direction and then coefficients were interpolated again in perpendicular direction
!! \param splcoeff2c array of spline coefficients for which distribution was first interpolated in perpendicular direction and then coefficients were interpolated again in parallel direction
subroutine get_splinecoeff(iarb,splcoeff1c,splcoeff2c)
  use param_mod
  implicit none
  integer :: iarb
  integer :: ipara,iperp

  real, dimension (npara(iarb),nperp(iarb),4) :: splcoeff1a
  real, dimension (npara(iarb),nperp(iarb),4,3) :: splcoeff1b

  real, dimension (npara(iarb),nperp(iarb),4) :: splcoeff2a
  real, dimension (npara(iarb),nperp(iarb),4,3) :: splcoeff2b

  real, dimension (npara(iarb)-1,nperp(iarb)-1,4,3) :: splcoeff1c
  real, dimension (npara(iarb)-1,nperp(iarb)-1,4,3) :: splcoeff2c

  !Interpolate distribution over perpendicular velocity while scanning through the parallel velocity grid

  do ipara=1,npara(iarb)

     do iperp=1,nperp(iarb)
        splcoeff1a(ipara,iperp,1)=distribution(ipara,iperp,iarb)
     enddo

     call spline_interpol(vperp(:,iarb),splcoeff1a(ipara,:,1),splcoeff1a(ipara,:,4),&
          & splcoeff1a(ipara,:,3),splcoeff1a(ipara,:,2),nperp(iarb))

  enddo

  !Interpolate coefficients of perpendicular spline interpolation over parallel velocity

  do iperp=1,nperp(iarb)-1

     do ipara=1,npara(iarb)

        splcoeff1b(ipara,iperp,1,1)= 3*splcoeff1a(ipara,iperp,2)
        splcoeff1b(ipara,iperp,1,2)= 2*splcoeff1a(ipara,iperp,3)-6*splcoeff1a(ipara,iperp,2)*vperp(iperp,iarb)
        splcoeff1b(ipara,iperp,1,3)= splcoeff1a(ipara,iperp,4)-2*splcoeff1a(ipara,iperp,3)*vperp(iperp,iarb)+&
             & 3*splcoeff1a(ipara,iperp,2)*vperp(iperp,iarb)*vperp(iperp,iarb)
     enddo

     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,1,1),splcoeff1b(:,iperp,4,1),&
          & splcoeff1b(:,iperp,3,1),splcoeff1b(:,iperp,2,1),npara(iarb))
     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,1,2),splcoeff1b(:,iperp,4,2),&
          & splcoeff1b(:,iperp,3,2),splcoeff1b(:,iperp,2,2),npara(iarb))
     call spline_interpol(vpara(:,iarb),splcoeff1b(:,iperp,1,3),splcoeff1b(:,iperp,4,3),&
          & splcoeff1b(:,iperp,3,3),splcoeff1b(:,iperp,2,3),npara(iarb))

     do ipara=1,npara(iarb)-1

        splcoeff1c(ipara,iperp,1,1)=splcoeff1b(ipara,iperp,2,1)
        splcoeff1c(ipara,iperp,2,1)=splcoeff1b(ipara,iperp,3,1)-3*splcoeff1b(ipara,iperp,2,1)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,3,1)=splcoeff1b(ipara,iperp,4,1)-2*splcoeff1b(ipara,iperp,3,1)*vpara(ipara,iarb)+&
             &  3*splcoeff1b(ipara,iperp,2,1)*vpara(ipara,iarb)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,4,1)=splcoeff1b(ipara,iperp,1,1)-splcoeff1b(ipara,iperp,4,1)*vpara(ipara,iarb)+&
             & splcoeff1b(ipara,iperp,3,1)*vpara(ipara,iarb)*vpara(ipara,iarb)-&
             & splcoeff1b(ipara,iperp,2,1)*vpara(ipara,iarb)*vpara(ipara,iarb)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,1,2)=splcoeff1b(ipara,iperp,2,2)
        splcoeff1c(ipara,iperp,2,2)=splcoeff1b(ipara,iperp,3,2)-3*splcoeff1b(ipara,iperp,2,2)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,3,2)=splcoeff1b(ipara,iperp,4,2)-2*splcoeff1b(ipara,iperp,3,2)*vpara(ipara,iarb)+&
             &  3*splcoeff1b(ipara,iperp,2,2)*vpara(ipara,iarb)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,4,2)=splcoeff1b(ipara,iperp,1,2)-splcoeff1b(ipara,iperp,4,2)*vpara(ipara,iarb)+&
             & splcoeff1b(ipara,iperp,3,2)*vpara(ipara,iarb)*vpara(ipara,iarb)-&
             & splcoeff1b(ipara,iperp,2,2)*vpara(ipara,iarb)*vpara(ipara,iarb)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,1,3)=splcoeff1b(ipara,iperp,2,3)
        splcoeff1c(ipara,iperp,2,3)=splcoeff1b(ipara,iperp,3,3)-3*splcoeff1b(ipara,iperp,2,3)*vpara(ipara,iarb)

        splcoeff1c(ipara,iperp,3,3)=splcoeff1b(ipara,iperp,4,3)-2*splcoeff1b(ipara,iperp,3,3)*vpara(ipara,iarb)+&
             &  3*splcoeff1b(ipara,iperp,2,3)*vpara(ipara,iarb)*vpara(ipara,iarb)
        splcoeff1c(ipara,iperp,4,3)=splcoeff1b(ipara,iperp,1,3)-splcoeff1b(ipara,iperp,4,3)*vpara(ipara,iarb)+&
             & splcoeff1b(ipara,iperp,3,3)*vpara(ipara,iarb)*vpara(ipara,iarb)-&
             & splcoeff1b(ipara,iperp,2,3)*vpara(ipara,iarb)*vpara(ipara,iarb)*vpara(ipara,iarb)

     enddo

  enddo




  !Now, interpolate distribution over parallel velocity while scanning through the perpendicular velocity grid

  do iperp=1,nperp(iarb)

     do ipara=1,npara(iarb)
        splcoeff2a(ipara,iperp,1)=distribution(ipara,iperp,iarb)
     enddo

     call spline_interpol(vpara(:,iarb),splcoeff2a(:,iperp,1),splcoeff2a(:,iperp,4),&
          & splcoeff2a(:,iperp,3),splcoeff2a(:,iperp,2),npara(iarb))

  enddo


  !Interpolate coefficients of parallel spline interpolation over perpendicular velocity

  do ipara=1,npara(iarb)-1

     do iperp=1,nperp(iarb)
        splcoeff2b(ipara,iperp,1,1)=3*splcoeff2a(ipara,iperp,2)
        splcoeff2b(ipara,iperp,1,2)=2*splcoeff2a(ipara,iperp,3)-6*splcoeff2a(ipara,iperp,2)*vpara(ipara,iarb)
        splcoeff2b(ipara,iperp,1,3)=splcoeff2a(ipara,iperp,4)-2*splcoeff2a(ipara,iperp,3)*vpara(ipara,iarb)+&
             & 3*splcoeff2a(ipara,iperp,2)*vpara(ipara,iarb)*vpara(ipara,iarb)

     enddo

     call spline_interpol(vperp(:,iarb),splcoeff2b(ipara,:,1,1),splcoeff2b(ipara,:,4,1),&
          & splcoeff2b(ipara,:,3,1),splcoeff2b(ipara,:,2,1),nperp(iarb))
     call spline_interpol(vperp(:,iarb),splcoeff2b(ipara,:,1,2),splcoeff2b(ipara,:,4,2),&
          & splcoeff2b(ipara,:,3,2),splcoeff2b(ipara,:,2,2),nperp(iarb))
     call spline_interpol(vperp(:,iarb),splcoeff2b(ipara,:,1,3),splcoeff2b(ipara,:,4,3),&
          & splcoeff2b(ipara,:,3,3),splcoeff2b(ipara,:,2,3),nperp(iarb))


     do iperp=1,nperp(iarb)-1

        splcoeff2c(ipara,iperp,1,1)=splcoeff2b(ipara,iperp,2,1)
        splcoeff2c(ipara,iperp,2,1)=splcoeff2b(ipara,iperp,3,1)-3*splcoeff2b(ipara,iperp,2,1)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,3,1)=splcoeff2b(ipara,iperp,4,1)-2*splcoeff2b(ipara,iperp,3,1)*vperp(iperp,iarb)+&
             &  3*splcoeff2b(ipara,iperp,2,1)*vperp(iperp,iarb)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,4,1)=splcoeff2b(ipara,iperp,1,1)-splcoeff2b(ipara,iperp,4,1)*vperp(iperp,iarb)+&
             & splcoeff2b(ipara,iperp,3,1)*vperp(iperp,iarb)*vperp(iperp,iarb)-&
             & splcoeff2b(ipara,iperp,2,1)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)

        splcoeff2c(ipara,iperp,1,2)=splcoeff2b(ipara,iperp,2,2)
        splcoeff2c(ipara,iperp,2,2)=splcoeff2b(ipara,iperp,3,2)-3*splcoeff2b(ipara,iperp,2,2)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,3,2)=splcoeff2b(ipara,iperp,4,2)-2*splcoeff2b(ipara,iperp,3,2)*vperp(iperp,iarb)+&
             &  3*splcoeff2b(ipara,iperp,2,2)*vperp(iperp,iarb)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,4,2)=splcoeff2b(ipara,iperp,1,2)-splcoeff2b(ipara,iperp,4,2)*vperp(iperp,iarb)+&
             & splcoeff2b(ipara,iperp,3,2)*vperp(iperp,iarb)*vperp(iperp,iarb)-&
             & splcoeff2b(ipara,iperp,2,2)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)

        splcoeff2c(ipara,iperp,1,3)=splcoeff2b(ipara,iperp,2,3)
        splcoeff2c(ipara,iperp,2,3)=splcoeff2b(ipara,iperp,3,3)-3*splcoeff2b(ipara,iperp,2,3)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,3,3)=splcoeff2b(ipara,iperp,4,3)-2*splcoeff2b(ipara,iperp,3,3)*vperp(iperp,iarb)+&
             &  3*splcoeff2b(ipara,iperp,2,3)*vperp(iperp,iarb)*vperp(iperp,iarb)
        splcoeff2c(ipara,iperp,4,3)=splcoeff2b(ipara,iperp,1,3)-splcoeff2b(ipara,iperp,4,3)*vperp(iperp,iarb)+&
             & splcoeff2b(ipara,iperp,3,3)*vperp(iperp,iarb)*vperp(iperp,iarb)-&
             & splcoeff2b(ipara,iperp,2,3)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)

     enddo

  enddo


end subroutine get_splinecoeff

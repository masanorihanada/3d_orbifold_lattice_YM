!***********************************************************
!***********************************************************
!    Polaakov loop 
SUBROUTINE Calc_Polyakov(umat,Pol_re,Pol_im)
  
  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex trace,MAT1(1:nmat,1:nmat),MAT2(1:nmat,1:nmat),eig_P(1:nmat)
  double precision Pol_re,Pol_im,pi
  double complex pol_pre
  integer imat,jmat,kmat
  integer isite,ix,iy,it

  pi=2d0*dasin(1d0)

  pol_pre=(0d0,0d0)
  do ix=1,nx
     do iy=1,ny
        do imat=1,nmat
           do jmat=1,nmat
              MAT1(imat,jmat)=umat(imat,jmat,1,ix,iy)
           end do
        end do
        do it=2,nt
           MAT2=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    MAT2(imat,jmat)=MAT2(imat,jmat)&
                         +MAT1(imat,kmat)*umat(kmat,jmat,it,ix,iy)
                 end do
              end do
           end do
           MAT1=MAT2
        end do
        trace=(0d0,0d0)
        do imat=1,nmat
           trace=trace+MAT2(imat,imat)
        end do
        Pol_pre=Pol_pre+trace

     end do
  end do
  Pol_re=dble(Pol_pre)/dble(NMAT*nx*ny)
  Pol_im=dble(Pol_pre*(0d0,-1d0))/dble(NMAT*nx*ny)

  return
  
END SUBROUTINE Calc_Polyakov

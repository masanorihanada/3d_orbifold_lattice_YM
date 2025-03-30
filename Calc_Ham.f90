!**********************
!*** ndim must be 3 ***
!**********************
subroutine Calc_Ham(at,as,mass2,mass2_U1,umat,zmat,P_umat,P_zmat,ham)

  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)
  double complex P_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny)
  double complex P_zmat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim-1)
  double precision ham
  double precision at,as
  double precision mass2,mass2_U1

  integer it,ix,iy
  integer idim
  integer imat,jmat

  call Calc_action(at,as,mass2,mass2_U1,umat,zmat,ham)

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           do imat=1,nmat
              do jmat=1,nmat
                 ham=ham&
                      +dble(P_umat(imat,jmat,it,ix,iy)&
                      *P_umat(jmat,imat,it,ix,iy))*0.5d0
                 
                 do idim=1,ndim-1
                    ham=ham&
                         +dble(P_zmat(imat,jmat,it,ix,iy,idim)&
                         *dconjg(P_zmat(imat,jmat,it,ix,iy,idim)))
                 end do
              end do
           end do
        end do
     end do
  end do
  
  return
  
END subroutine Calc_Ham

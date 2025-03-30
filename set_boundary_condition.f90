! umat -> exp(i*P_umat*dtau_umat)*umat
! xmat -> xmat + P_xmat*dtau_xmat
! P_umat -> P_umat - delh_umat*dtau_umat
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dS/dxmat(jmat,imat)

subroutine set_bc(umat,zmat)

  implicit none

  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)

  integer imat,jmat,it,ix,iy,idim


  do imat=1,nmat
     do jmat=1,nmat
        !do it=1,nt
        do ix=1,nx
           do iy=1,ny
              umat(imat,jmat,0,ix,iy)=umat(imat,jmat,nt,ix,iy)
              umat(imat,jmat,nt+1,ix,iy)=umat(imat,jmat,1,ix,iy)
              do idim=1,ndim-1
                 zmat(imat,jmat,0,ix,iy,idim)=zmat(imat,jmat,nt,ix,iy,idim)
                 zmat(imat,jmat,nt+1,ix,iy,idim)=zmat(imat,jmat,1,ix,iy,idim)
              end do
           end do
        end do
           
        do it=1,nt
           !do ix=1,nx
           do iy=1,ny
              umat(imat,jmat,it,0,iy)=umat(imat,jmat,it,nx,iy)
              umat(imat,jmat,it,nx+1,iy)=umat(imat,jmat,it,1,iy)
              do idim=1,ndim-1
                 zmat(imat,jmat,it,0,iy,idim)=zmat(imat,jmat,it,nx,iy,idim)
                 zmat(imat,jmat,it,nx+1,iy,idim)=zmat(imat,jmat,it,1,iy,idim)
              end do
           end do
        end do
           
        do it=1,nt
           do ix=1,nx
              !do iy=1,ny
              umat(imat,jmat,it,ix,0)=umat(imat,jmat,it,ix,ny)
              umat(imat,jmat,it,ix,ny+1)=umat(imat,jmat,it,ix,1)
              do idim=1,ndim-1
                 zmat(imat,jmat,it,ix,0,idim)=zmat(imat,jmat,it,ix,ny,idim)
                 zmat(imat,jmat,it,ix,ny+1,idim)=zmat(imat,jmat,it,ix,1,idim)
              end do
           end do
        end do

        do it=1,nt
           umat(imat,jmat,it,0,0)=umat(imat,jmat,it,nx,ny)
           umat(imat,jmat,it,0,ny+1)=umat(imat,jmat,it,nx,1)
           umat(imat,jmat,it,nx+1,0)=umat(imat,jmat,it,1,ny)
           umat(imat,jmat,it,nx+1,ny+1)=umat(imat,jmat,it,1,1)
           do idim=1,ndim-1
              zmat(imat,jmat,it,0,0,idim)=zmat(imat,jmat,it,nx,ny,idim)
              zmat(imat,jmat,it,0,ny+1,idim)=zmat(imat,jmat,it,nx,1,idim)
              zmat(imat,jmat,it,nx+1,0,idim)=zmat(imat,jmat,it,1,ny,idim)
              zmat(imat,jmat,it,nx+1,ny+1,idim)=zmat(imat,jmat,it,1,1,idim)
           end do
        end do
     
        do ix=1,nx
           umat(imat,jmat,0,ix,0)=umat(imat,jmat,nt,ix,ny)
           umat(imat,jmat,nt+1,ix,0)=umat(imat,jmat,1,ix,ny)
           umat(imat,jmat,0,ix,ny+1)=umat(imat,jmat,nt,ix,1)
           umat(imat,jmat,nt+1,ix,ny+1)=umat(imat,jmat,1,ix,1)
           do idim=1,ndim-1
              zmat(imat,jmat,0,ix,0,idim)=zmat(imat,jmat,nt,ix,ny,idim)
              zmat(imat,jmat,nt+1,ix,0,idim)=zmat(imat,jmat,1,ix,ny,idim)
              zmat(imat,jmat,0,ix,ny+1,idim)=zmat(imat,jmat,nt,ix,1,idim)
              zmat(imat,jmat,nt+1,ix,ny+1,idim)=zmat(imat,jmat,1,ix,1,idim)
           end do
        end do

        do iy=1,ny
           umat(imat,jmat,0,0,iy)=umat(imat,jmat,nt,nx,iy)
           umat(imat,jmat,nt+1,0,iy)=umat(imat,jmat,1,nx,iy)
           umat(imat,jmat,0,nx+1,iy)=umat(imat,jmat,nt,1,iy)
           umat(imat,jmat,nt+1,nx+1,iy)=umat(imat,jmat,1,1,iy)
           do idim=1,ndim-1
              zmat(imat,jmat,0,0,iy,idim)=zmat(imat,jmat,nt,nx,iy,idim)
              zmat(imat,jmat,nt+1,0,iy,idim)=zmat(imat,jmat,1,nx,iy,idim)
              zmat(imat,jmat,0,nx+1,iy,idim)=zmat(imat,jmat,nt,1,iy,idim)
              zmat(imat,jmat,nt+1,nx+1,iy,idim)=zmat(imat,jmat,1,1,iy,idim)
           end do
        end do

        umat(imat,jmat,0,0,0)=umat(imat,jmat,nt,nx,ny)
        umat(imat,jmat,nt+1,0,0)=umat(imat,jmat,1,nx,ny)
        umat(imat,jmat,0,nx+1,0)=umat(imat,jmat,nt,1,ny)
        umat(imat,jmat,0,0,ny+1)=umat(imat,jmat,nt,nx,1)
        umat(imat,jmat,nt+1,nx+1,0)=umat(imat,jmat,1,1,ny)
        umat(imat,jmat,0,nx+1,ny+1)=umat(imat,jmat,nt,1,1)
        umat(imat,jmat,nt+1,0,ny+1)=umat(imat,jmat,1,nx,1)
        umat(imat,jmat,nt+1,nx+1,ny+1)=umat(imat,jmat,1,1,1)
        do idim=1,ndim-1
           zmat(imat,jmat,0,0,0,idim)=zmat(imat,jmat,nt,nx,ny,idim)
           zmat(imat,jmat,nt+1,0,0,idim)=zmat(imat,jmat,1,nx,ny,idim)
           zmat(imat,jmat,0,nx+1,0,idim)=zmat(imat,jmat,nt,1,ny,idim)
           zmat(imat,jmat,0,0,ny+1,idim)=zmat(imat,jmat,nt,nx,1,idim)
           zmat(imat,jmat,nt+1,nx+1,0,idim)=zmat(imat,jmat,1,1,ny,idim)
           zmat(imat,jmat,0,nx+1,ny+1,idim)=zmat(imat,jmat,nt,1,1,idim)
           zmat(imat,jmat,nt+1,0,ny+1,idim)=zmat(imat,jmat,1,nx,1,idim)
           zmat(imat,jmat,nt+1,nx+1,ny+1,idim)=zmat(imat,jmat,1,1,1,idim)
        end do
           
     end do
  end do
  
  return

END subroutine Set_Bc

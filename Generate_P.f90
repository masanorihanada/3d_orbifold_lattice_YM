! Generate P with Gaussian weight.
SUBROUTINE Generate_P(P_umat,P_zmat)
  
  implicit none
  include 'size.h'
  
  integer it,ix,iy,idim,imat,jmat
  double precision r1,r2,trace

  double complex P_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny)
  double complex P_zmat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim-1)
  !**************
  !*** P_umat ***
  !**************
  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           do imat=1,nmat-1
              do jmat=imat+1,nmat
                 call BoxMuller(r1,r2)
                 P_umat(imat,jmat,it,ix,iy)=&
                      dcmplx(r1/dsqrt(2d0))+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                 P_umat(jmat,imat,it,ix,iy)=&
                      dcmplx(r1/dsqrt(2d0))-dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
              end do
           end do
           trace=0d0
           do imat=1,nmat
              call BoxMuller(r1,r2)
              P_umat(imat,imat,it,ix,iy)=dcmplx(r1)
              trace=trace+r1
           end do
           !traceless condition
           do imat=1,nmat
              call BoxMuller(r1,r2)
              P_umat(imat,imat,it,ix,iy)=P_umat(imat,imat,it,ix,iy)-dcmplx(trace/dble(nmat))
           end do
        end do
     end do
  end do
  !**************
  !*** P_zmat ***
  !**************
  do idim=1,ndim-1
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat
                 do jmat=1,nmat
                    call BoxMuller(r1,r2)
                    P_zmat(imat,jmat,it,ix,iy,idim)=&
                         dcmplx(r1/dsqrt(2d0))+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                 end do
              end do
           end do
        end do
     end do
  end do
   
  return
  
END SUBROUTINE Generate_P

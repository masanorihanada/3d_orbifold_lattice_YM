!**********************
!*** ndim must be 3 ***
!**********************
subroutine measurements(as,zmat,plaq_U,plaq_Z,Wm1,av_det_U)

  implicit none
  include 'size.h'
  
  double complex zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)
  double complex wmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)
  double complex umat_spatial(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)!this umat is dummy, needed for set_bc subroutine. 
  double precision as !not really used for 3d theory

  double precision plaq_U,plaq_Z,Wm1

  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat

  double complex UU(1:nmat,1:nmat)
  double complex UdUd(1:nmat,1:nmat)
  double complex ZxZy(1:nmat,1:nmat)
  double complex ZxbarZybar(1:nmat,1:nmat)
  
  
  double complex MAT(1:nmat,1:nmat)
  double complex MATSQRT(1:nmat,1:nmat)
  double complex temp
  double complex av_det_U,det_U

  
  plaq_U = 0d0
  plaq_Z = 0d0
  Wm1 = 0d0

  wmat = (0d0,0d0)
  umat_spatial = (0d0,0d0)
  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           do idim=1,ndim-1

              MAT=(0d0,0d0)
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       MAT(imat,jmat) = MAT(imat,jmat) + zmat(imat,kmat,it,ix,iy,idim) * dconjg(zmat(jmat,kmat,it,ix,iy,idim))
                    end do
                 end do
              end do
              call MATRIX_SQRT(NMAT,MAT,MATSQRT)
              do imat=1,nmat
                 do jmat=1,nmat
                    wmat(imat,jmat,it,ix,iy,idim) = MATSQRT(imat,jmat) * dcmplx(dsqrt(2d0))
                 end do
              end do
              call MatrixInverse(NMAT,MATSQRT)
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       umat_spatial(imat,jmat,it,ix,iy,idim) = umat_spatial(imat,jmat,it,ix,iy,idim) + MATSQRT(imat,kmat) * zmat(kmat,jmat,it,ix,iy,idim)
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  !call set_bc(umat,wmat)
  call set_bc(umat,umat_spatial)

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           UU = (0d0,0d0)
           UdUd = (0d0,0d0)
           ZxZy = (0d0,0d0)
           ZxbarZybar = (0d0,0d0)

  
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    UU(imat,jmat)=UU(imat,jmat)&
                         &+umat_spatial(imat,kmat,it,ix,iy,1)*umat_spatial(kmat,jmat,it,ix+1,iy,2)
                    UdUd(imat,jmat)=UdUd(imat,jmat)&
                         &+dconjg(umat_spatial(kmat,imat,it,ix,iy+1,1))*dconjg(umat_spatial(jmat,kmat,it,ix,iy,2))

                    ZxZy(imat,jmat)=ZxZy(imat,jmat)&
                         &+zmat(imat,kmat,it,ix,iy,1)*zmat(kmat,jmat,it,ix+1,iy,2)
                    
                    ZxbarZybar(imat,jmat)=ZxbarZybar(imat,jmat)&
                         &+dconjg(zmat(kmat,imat,it,ix,iy+1,1))*dconjg(zmat(jmat,kmat,it,ix,iy,2))

                 end do
              end do
           end do

           do imat=1,nmat
              do jmat=1,nmat
                 plaq_U = plaq_U + dble(UU(imat,jmat)*UdUd(jmat,imat))
                 plaq_Z = plaq_Z + dble(ZxZy(imat,jmat)*ZxbarZybar(jmat,imat))
              end do
           end do

        end do
     end do
  end do
  plaq_U = plaq_U / dble(nx * ny * nt)
  plaq_Z = plaq_Z / dble(nx * ny * nt)

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           do idim=1,ndim-1

              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat) = wmat(imat,jmat,it,ix,iy,idim) 
                 end do
              end do
              do imat=1,nmat
                 MAT(imat,imat) = MAT(imat,imat) - (1d0,0d0)
              end do

              do imat=1,nmat
                 do jmat=1,nmat
                    Wm1 = Wm1 + dble(MAT(imat,jmat)*MAT(jmat,imat))
                 end do
              end do
           end do
        end do
     end do
  end do

  Wm1 = Wm1 / dble( nx * ny * nt * (ndim - 1) )

  av_det_U = (0d0,0d0)
  do idim = 1,ndim-1
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              !MAT=(0d0,0d0)
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat) = umat_spatial(imat,jmat,it,ix,iy,idim)
                 end do
              end do
              call MATRIX_DETERMINANT_COMPLEX(nmat,MAT,det_U)
              av_det_U = av_det_U + det_U
           end do
        end do
     end do
  end do
  av_det_U = av_det_U / dble( nx * ny * nt * (ndim - 1) )
  
  return

END subroutine measurements

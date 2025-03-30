!**********************
!*** ndim must be 3 ***
!**********************
!
! umat -> exp(i*P_umat*dtau)*umat
! P_umat -> P_umat - delh_umat*dtau
!
! zmat -> zmat + P_zmat*dtay
! P_umat -> P_umat - delh_umat*dtau

subroutine Molecular_Dynamics(at,as,mass2,mass2_U1,ntau,dtau_t,dtau_s,umat,zmat,P_umat,P_zmat)

  implicit none

  include 'size.h'

  integer ntau

  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)
  
  double precision at,as
  double precision mass2,mass2_U1
  double precision dtau_t,dtau_s
  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat
  integer itau
  
  double complex P_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny)
  double complex P_zmat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim-1)
  
  double complex delh_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny)
  double complex delh_zmat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim-1)
  
  double complex MAT(1:NMAT,1:NMAT),EXP_iP(1:NMAT,1:NMAT)

  ! first step of leap frog
  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           do imat=1,nmat
              do jmat=1,nmat
                 MAT(imat,jmat)=P_umat(imat,jmat,it,ix,iy)*dcmplx(0.5d0*dtau_t)
              end do
           end do
           call MATRIX_iEXP(NMAT,MAT,EXP_iP)
           do imat=1,nmat
              do jmat=1,nmat
                 MAT(imat,jmat)=umat(imat,jmat,it,ix,iy)
                 umat(imat,jmat,it,ix,iy)=(0d0,0d0)
              end do
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    umat(imat,jmat,it,ix,iy)=umat(imat,jmat,it,ix,iy)&
                         +EXP_iP(imat,kmat)*MAT(kmat,jmat)
                 end do
                 do idim=1,ndim-1
                    zmat(imat,jmat,it,ix,iy,idim)&
                         =zmat(imat,jmat,it,ix,iy,idim)&
                         +P_zmat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau_s*0.5d0)
                 end do
              end do
           end do
        end do
     end do
  end do
  call set_bc(umat,zmat)
  ! second,...,Ntau-th   
  do itau=2,ntau
     call Calc_DELH(at,as,mass2,mass2_U1,umat,zmat,delh_umat,delh_zmat)
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat
                 do jmat=1,nmat
                    P_umat(imat,jmat,it,ix,iy)&
                         &=P_umat(imat,jmat,it,ix,iy)&
                         &-delh_umat(imat,jmat,it,ix,iy)*dcmplx(dtau_t)
                    do idim=1,ndim-1
                       P_zmat(imat,jmat,it,ix,iy,idim)&
                            &=P_zmat(imat,jmat,it,ix,iy,idim)&
                            &-delh_zmat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau_s)
                    end do
                 end do
              end do
           end do
        end do
     end do
     
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat)=P_umat(imat,jmat,it,ix,iy)*dcmplx(dtau_t)
                 end do
              end do
              call MATRIX_iEXP(NMAT,MAT,EXP_iP)
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat)=umat(imat,jmat,it,ix,iy)
                    umat(imat,jmat,it,ix,iy)=(0d0,0d0)
                 end do
              end do
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       umat(imat,jmat,it,ix,iy)=umat(imat,jmat,it,ix,iy)&
                            +EXP_iP(imat,kmat)*MAT(kmat,jmat)
                    end do
                    do idim=1,ndim-1
                       zmat(imat,jmat,it,ix,iy,idim)&
                            =zmat(imat,jmat,it,ix,iy,idim)&
                            +P_zmat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau_s)
                    end do
                 end do
              end do
                 
           end do
        end do
     end do
     call set_bc(umat,zmat)   
  end do
  ! last step
  call Calc_DELH(at,as,mass2,mass2_U1,umat,zmat,delh_umat,delh_zmat)
  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           do imat=1,nmat
              do jmat=1,nmat
                 P_umat(imat,jmat,it,ix,iy)&
                      &=P_umat(imat,jmat,it,ix,iy)&
                      &-delh_umat(imat,jmat,it,ix,iy)*dcmplx(dtau_t)
                 do idim=1,ndim-1
                    P_zmat(imat,jmat,it,ix,iy,idim)&
                         &=P_zmat(imat,jmat,it,ix,iy,idim)&
                         &-delh_zmat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau_s)
                 end do
              end do
           end do
        end do
     end do
  end do
  
  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           do imat=1,nmat
              do jmat=1,nmat
                 MAT(imat,jmat)=P_umat(imat,jmat,it,ix,iy)*dcmplx(0.5d0*dtau_t)
              end do
           end do
           call MATRIX_iEXP(NMAT,MAT,EXP_iP)
           do imat=1,nmat
              do jmat=1,nmat
                 MAT(imat,jmat)=umat(imat,jmat,it,ix,iy)
                 umat(imat,jmat,it,ix,iy)=(0d0,0d0)
              end do
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    umat(imat,jmat,it,ix,iy)=umat(imat,jmat,it,ix,iy)&
                         +EXP_iP(imat,kmat)*MAT(kmat,jmat)
                 end do
                 do idim=1,ndim-1
                    zmat(imat,jmat,it,ix,iy,idim)&
                         =zmat(imat,jmat,it,ix,iy,idim)&
                         +P_zmat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau_s*0.5d0)
                 end do
              end do
           end do
           
        end do
     end do
  end do
  call set_bc(umat,zmat)

  return

END subroutine Molecular_Dynamics

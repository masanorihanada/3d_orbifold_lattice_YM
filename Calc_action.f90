!**********************
!*** ndim must be 3 ***
!**********************
subroutine Calc_action(at,as,mass2,mass2_U1,umat,zmat,action)

  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)

  double precision action,kin,pot,pot_scalar_mass
  double precision at,as
  double precision mass2,mass2_U1

  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat

  double complex UZ(1:nmat,1:nmat,1:ndim-1)
  double complex ZU(1:nmat,1:nmat,1:ndim-1)
  double complex ZxZy(1:nmat,1:nmat)
  double complex ZyZx(1:nmat,1:nmat)
  double complex ZxZxbar(1:nmat,1:nmat)
  double complex ZxbarZx(1:nmat,1:nmat)
  double complex ZyZybar(1:nmat,1:nmat)
  double complex ZybarZy(1:nmat,1:nmat)

  double complex MAT(1:nmat,1:nmat)
  double complex temp

  double complex detZ

  
  kin = 0d0
  pot = 0d0
  pot_scalar_mass = 0d0

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           ZxZy = (0d0,0d0)
           ZyZx = (0d0,0d0)
           UZ = (0d0,0d0)
           ZU = (0d0,0d0)
           ZxZxbar = (0d0,0d0)
           ZyZybar = (0d0,0d0)
           ZxbarZx = (0d0,0d0)
           ZybarZy = (0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    UZ(imat,jmat,1)=UZ(imat,jmat,1)&
                         &+umat(imat,kmat,it,ix,iy)*zmat(kmat,jmat,it+1,ix,iy,1)

                    UZ(imat,jmat,2)=UZ(imat,jmat,2)&
                         &+umat(imat,kmat,it,ix,iy)*zmat(kmat,jmat,it+1,ix,iy,2)
                    
                    ZU(imat,jmat,1)=ZU(imat,jmat,1)&
                         &+zmat(imat,kmat,it,ix,iy,1)*umat(kmat,jmat,it,ix+1,iy)

                    ZU(imat,jmat,2)=ZU(imat,jmat,2)&
                         &+zmat(imat,kmat,it,ix,iy,2)*umat(kmat,jmat,it,ix,iy+1)
                    
                    ZxZy(imat,jmat)=ZxZy(imat,jmat)&
                         &+zmat(imat,kmat,it,ix,iy,1)*zmat(kmat,jmat,it,ix+1,iy,2)
                    
                    ZyZx(imat,jmat)=ZyZx(imat,jmat)&
                         &+zmat(imat,kmat,it,ix,iy,2)*zmat(kmat,jmat,it,ix,iy+1,1)

                    ZxZxbar(imat,jmat)=ZxZxbar(imat,jmat)&
                         &+zmat(imat,kmat,it,ix,iy,1)*dconjg(zmat(jmat,kmat,it,ix,iy,1))
                    
                    ZyZybar(imat,jmat)=ZyZybar(imat,jmat)&
                         &+zmat(imat,kmat,it,ix,iy,2)*dconjg(zmat(jmat,kmat,it,ix,iy,2))
                    
                    ZxbarZx(imat,jmat)=ZxbarZx(imat,jmat)&
                         &+dconjg(zmat(kmat,imat,it,ix-1,iy,1))*zmat(kmat,jmat,it,ix-1,iy,1)

                    ZybarZy(imat,jmat)=ZybarZy(imat,jmat)&
                         &+dconjg(zmat(kmat,imat,it,ix,iy-1,2))*zmat(kmat,jmat,it,ix,iy-1,2)

                    
                 end do
              end do
           end do

           do imat=1,nmat
              do jmat=1,nmat
                 temp = UZ(imat,jmat,1) - ZU(imat,jmat,1)
                 kin = kin+dble(temp * dconjg(temp))
                 temp = UZ(imat,jmat,2) - ZU(imat,jmat,2)
                 kin = kin+dble(temp * dconjg(temp))
 
                 temp = ZxZy(imat,jmat) - ZyZx(imat,jmat)
                 pot = pot + dble(temp * dconjg(temp)) * 2d0

                 temp = ZxZxbar(imat,jmat) - ZxbarZx(imat,jmat)&
                      + ZyZybar(imat,jmat) - ZybarZy(imat,jmat)
                 pot = pot + dble(temp * dconjg(temp)) * 0.5d0

                 temp = ZxZxbar(imat,jmat)
                 if(imat.EQ.jmat)then
                    temp = temp - (0.5d0,0d0)
                 end if
                 pot_scalar_mass = pot_scalar_mass + dble(temp * dconjg(temp))

                 temp = ZyZybar(imat,jmat)
                 if(imat.EQ.jmat)then
                    temp = temp - (0.5d0,0d0)
                 end if
                 pot_scalar_mass = pot_scalar_mass + dble(temp * dconjg(temp))

              end do
           end do

        end do
     end do
  end do
  
  kin = kin / at
  pot = pot * at / ( as * as )
  pot_scalar_mass = pot_scalar_mass * 0.5d0 * mass2 * at
  action = kin + pot + pot_scalar_mass

  !*****************************************
  !*** constraint terms with determinant ***
  !*****************************************

  do idim=1,ndim-1
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat)=zmat(imat,jmat,it,ix,iy,idim)
                 end do
              end do
              call MATRIX_DETERMINANT_COMPLEX(nmat,MAT,detZ)
              temp = detZ * dcmplx(2.0**(dble(nmat)*0.5)) - (1d0,0d0)
              
              action = action + 0.5d0 * mass2_U1 * at * dble( temp * dconjg(temp) )
           end do
        end do
     end do
  end do


  
  return

END subroutine Calc_action

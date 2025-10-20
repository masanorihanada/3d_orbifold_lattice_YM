! umat -> exp(i*P_umat*dtau_umat)*umat
! xmat -> xmat + P_xmat*dtau_xmat
! P_umat -> P_umat - delh_umat*dtau_umat
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dS/dxmat(jmat,imat)

SUBROUTINE Calc_DELH(at,as,mass2,mass2_U1,umat,zmat,delh_umat,delh_zmat)

  implicit none
  include 'size.h'

  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)
  double complex delh_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny)
  double complex delh_zmat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim-1)
  double precision at,as
  double precision mass2,mass2_U1

  integer ix,iy,it,ix2,iy2
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
  double complex temp_mat(1:nmat,1:nmat),Zbar_times_temp_mat(1:nmat,1:nmat),zinv(1:nmat,1:nmat),trace

  double complex coeff
  integer it_p1,it_m1
  integer ix_p1,ix_m1
  integer iy_p1,iy_m1
  double complex detZ

  delh_umat=(0d0,0d0)
  delh_zmat=(0d0,0d0)

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           if(it.EQ.nt)then
              it_p1 = 1
           else
              it_p1 = it + 1
           end if
           if(it.EQ.1)then
              it_m1 = nt 
           else
              it_m1 = it - 1
           end if
           
           if(ix.EQ.nx)then
              ix_p1 = 1
           else
              ix_p1 = ix + 1
           end if
           if(ix.EQ.1)then
              ix_m1 = nx 
           else
              ix_m1 = ix - 1
           end if
            
           if(iy.EQ.ny)then
              iy_p1 = 1
           else
              iy_p1 = iy + 1
           end if
           if(iy.EQ.1)then
              iy_m1 = ny 
           else
              iy_m1 = iy - 1
           end if
           
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
           !***************************
           !*** delh_zmat, 1st line ***
           !***************************
           coeff = dcmplx(1d0/at)
           do idim=1,ndim-1
              do imat=1,nmat
                 do jmat=1,nmat
                    temp_mat(imat,jmat) = UZ(imat,jmat,idim) - ZU(imat,jmat,idim)
                 end do
              end do
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       delh_zmat(imat,jmat,it_p1,ix,iy,idim)=&
                            delh_zmat(imat,jmat,it_p1,ix,iy,idim)&
                            + dconjg(umat(kmat,imat,it,ix,iy))*temp_mat(kmat,jmat) * coeff
                    end do
                 end do
              end do
              
              if(idim.EQ.1)then
                 ix2 = ix_p1
                 iy2 = iy
              elseif(idim.EQ.2)then
                 ix2 = ix
                 iy2 = iy_p1
              end if
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       delh_zmat(imat,jmat,it,ix,iy,idim)=&
                            delh_zmat(imat,jmat,it,ix,iy,idim)&
                            - temp_mat(imat,kmat)*dconjg(umat(jmat,kmat,it,ix2,iy2)) * coeff
                       
                    end do
                 end do
              end do
              
           end do
           !***************************
           !*** delh_zmat, 2nd line ***
           !***************************
           temp_mat = ZxZxbar - ZxbarZx + ZyZybar - ZybarZy
           coeff = dcmplx(at/as/as)
           do idim = 1,ndim-1
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       delh_zmat(imat,jmat,it,ix,iy,idim)=&
                            delh_zmat(imat,jmat,it,ix,iy,idim)&
                            + temp_mat(imat,kmat) * zmat(kmat,jmat,it,ix,iy,idim) * coeff
                    end do
                 end do
              end do
           end do
           
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    delh_zmat(imat,jmat,it,ix_m1,iy,1)=&
                         delh_zmat(imat,jmat,it,ix_m1,iy,1)&
                         - zmat(imat,kmat,it,ix_m1,iy,1) * temp_mat(kmat,jmat) * coeff
                    
                    delh_zmat(imat,jmat,it,ix,iy_m1,2)=&
                         delh_zmat(imat,jmat,it,ix,iy_m1,2)&
                         - zmat(imat,kmat,it,ix,iy_m1,2) * temp_mat(kmat,jmat) * coeff
                 end do
              end do
           end do
           !***************************
           !*** delh_zmat, 3rd line ***
           !***************************
           temp_mat = ZxZy - ZyZx 
           coeff = dcmplx( 2d0 * at / as / as )

           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    delh_zmat(imat,jmat,it,ix,iy,1)=&
                         delh_zmat(imat,jmat,it,ix,iy,1)&
                         + temp_mat(imat,kmat) * dconjg(zmat(jmat,kmat,it,ix_p1,iy,2)) * coeff

                    delh_zmat(imat,jmat,it,ix_p1,iy,2)=&
                         delh_zmat(imat,jmat,it,ix_p1,iy,2)&
                         + dconjg(zmat(kmat,imat,it,ix,iy,1)) * temp_mat(kmat,jmat) * coeff

                    delh_zmat(imat,jmat,it,ix,iy,2)=&
                         delh_zmat(imat,jmat,it,ix,iy,2)&
                         - temp_mat(imat,kmat) * dconjg(zmat(jmat,kmat,it,ix,iy_p1,1)) * coeff

                    delh_zmat(imat,jmat,it,ix,iy_p1,1)=&
                         delh_zmat(imat,jmat,it,ix,iy_p1,1)&
                         - dconjg(zmat(kmat,imat,it,ix,iy,2)) * temp_mat(kmat,jmat) * coeff

                 end do
              end do
           end do
           
           !***********************************
           !*** delh_zmat, scalar mass term ***
           !***********************************
           coeff = dcmplx( mass2 * at )
           
           temp_mat = ZxZxbar
           do imat=1,nmat
              temp_mat(imat,imat) = temp_mat(imat,imat) - (0.5d0,0d0)
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    delh_zmat(imat,jmat,it,ix,iy,1)=&
                         delh_zmat(imat,jmat,it,ix,iy,1)&
                         + temp_mat(imat,kmat) * zmat(kmat,jmat,it,ix,iy,1) * coeff

                 end do
              end do
           end do

           temp_mat = ZyZybar
           do imat=1,nmat
              temp_mat(imat,imat) = temp_mat(imat,imat) - (0.5d0,0d0)
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    delh_zmat(imat,jmat,it,ix,iy,2)=&
                         delh_zmat(imat,jmat,it,ix,iy,2)&
                         + temp_mat(imat,kmat) * zmat(kmat,jmat,it,ix,iy,2) * coeff

                 end do
              end do
           end do
           !*****************
           !*** delh_umat ***
           !*****************
           coeff = dcmplx(1d0/at) * (0d0,1d0)
           
           do idim=1,ndim-1
              
              do imat=1,nmat
                 do jmat=1,nmat
                    temp_mat(imat,jmat) = UZ(imat,jmat,idim) - ZU(imat,jmat,idim)
                 end do
              end do
              
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       delh_umat(imat,jmat,it,ix,iy)=&
                            delh_umat(imat,jmat,it,ix,iy)&
                            + UZ(imat,kmat,idim) * dconjg(temp_mat(jmat,kmat)) * coeff                      
                    end do
                 end do
              end do
              
              !** actually the next one is just the c.c. of the last one **
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       delh_umat(imat,jmat,it,ix,iy)=&
                            delh_umat(imat,jmat,it,ix,iy)&
                            - temp_mat(imat,kmat) * dconjg(UZ(jmat,kmat,idim)) * coeff
                    end do
                 end do
              end do
              
              Zbar_times_temp_mat = (0d0,0d0)
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat   
                       Zbar_times_temp_mat(imat,jmat)=&
                            Zbar_times_temp_mat(imat,jmat)&
                            + dconjg(zmat(kmat,imat,it,ix,iy,idim)) * temp_mat(kmat,jmat)
                    end do
                 end do
              end do
              
              if(idim.EQ.1)then
                 ix2 = ix_p1
                 iy2 = iy
              else if(idim.EQ.2)then
                 ix2 = ix
                 iy2 = iy_p1
              end if
              ! let's re-use temp_mat for another purpose, i.e., for \bar{Z}_{ix,iy}(UZ-ZU)*\bar{U}_{ix2,iy2}
              temp_mat = (0d0,0d0)
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       temp_mat(imat,jmat) = temp_mat(imat,jmat)&
                            + Zbar_times_temp_mat(imat,kmat) * dconjg(umat(jmat,kmat,it,ix2,iy2))
                    end do
                 end do
              end do
              
              do imat=1,nmat
                 do jmat=1,nmat
                    delh_umat(imat,jmat,it,ix2,iy2)=&
                         delh_umat(imat,jmat,it,ix2,iy2)&
                         + ( temp_mat(imat,jmat) - dconjg(temp_mat(jmat,imat)) ) * coeff
                 end do
              end do
              
           end do
           
        end do
     end do
  end do
  !*****************************************
  !*** constraint terms with determinant ***
  !*****************************************
  do idim=1,ndim-1
     do it=1,nt
        do ix=1,nx
           do iy=1,ny

              do imat=1,nmat
                 do jmat=1,nmat
                    temp_mat(imat,jmat) = zmat(imat,jmat,it,ix,iy,idim)
                 end do
              end do
              zinv = temp_mat
              call MATRIX_DETERMINANT_COMPLEX(nmat,temp_mat,detZ)
              call MatrixInverse(nmat,zinv)
              coeff = dcmplx(  0.5d0 * mass2_U1 * at ) &
                   * ( detZ * dcmplx(2.0**(dble(nmat)*0.5)) - (1d0,0d0) )&
                   * dcmplx(2.0**(dble(nmat)*0.5)) * dconjg(detZ)

              !re-use temp_mat
              do imat=1,nmat
                 do jmat=1,nmat
                    temp_mat(imat,jmat) = dconjg( zinv(jmat,imat) ) 
                 end do
              end do
              temp_mat = temp_mat * coeff

              do imat=1,nmat
                 do jmat=1,nmat
                    delh_zmat(imat,jmat,it,ix,iy,idim)=&
                         delh_zmat(imat,jmat,it,ix,iy,idim) + temp_mat(imat,jmat) 
                 end do
              end do
              
           end do
        end do
     end do
  end do
  !*************************************
  !*** Traceless projection ************
  !*** (needed because Z is complex) ***
  !*************************************
  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           trace=(0d0,0d0)
           do imat=1,nmat
              trace = trace + delh_umat(imat,imat,it,ix,iy)
           end do
           do imat=1,nmat
              delh_umat(imat,imat,it,ix,iy) = delh_umat(imat,imat,it,ix,iy) - trace/dcmplx(nmat)
           end do
        end do
     end do
  end do


  
  return

END SUBROUTINE Calc_DELH

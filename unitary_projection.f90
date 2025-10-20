
subroutine SU_N_projection_all(umat)!projection to SU(N) matrix

  implicit none

  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1),u_temp(1:nmat,1:nmat)

  integer imat,jmat,it,ix,iy,idim

  do it=0,nt+1
     do ix=0,nx+1
        do iy=0,ny+1
           
           do imat=1,nmat
              do jmat=1,nmat
                 u_temp(imat,jmat)=umat(imat,jmat,it,ix,iy)
              end do
           end do
           
           call SU_N_projection_each(u_temp)

           do imat=1,nmat
              do jmat=1,nmat
                 umat(imat,jmat,it,ix,iy)=u_temp(imat,jmat)
              end do
           end do
           
        end do
     end do
  end do
  
  return

END subroutine SU_N_Projection_All
!**************************************************
subroutine SU_N_projection_each(u_temp)!projection to SU(N) matrix

  implicit none

  include 'size.h'
  
  double complex u_temp(1:nmat,1:nmat),u_temp_2(1:nmat,1:nmat),det,inner_product,phase
  double precision norm,norm_temp,arg,pi

  integer i_column,j_column,i_row,imat,jmat,kmat

  pi = dasin(1d0)*2d0
  !**************************
  !*** Projection to U(N) ***
  !**************************
  do i_column=1,nmat
     do j_column=1,i_column-1
        inner_product = (0d0,0d0)
        do i_row=1,nmat
           inner_product = inner_product + dconjg(u_temp(j_column,i_row)) * u_temp(i_column,i_row)
        end do
        do i_row=1,nmat
           u_temp(i_column,i_row) = u_temp(i_column,i_row) - inner_product * u_temp(j_column,i_row)
        end do
     end do
     !*******************************
     !** normalize a column vector **
     !*******************************
     norm_temp=0d0
     do i_row=1,nmat
        norm_temp = norm_temp + abs(u_temp(i_column,i_row))**2d0
     end do
     norm = dsqrt(norm_temp)
     do i_row=1,nmat
        u_temp(i_column,i_row) = u_temp(i_column,i_row) / dcmplx(norm)
     end do
  end do

  !***************************
  !*** Projection to SU(N) ***
  !***************************
  u_temp_2=u_temp
  call MATRIX_DETERMINANT_COMPLEX(nmat,u_temp_2,det)
  arg = dasin (dble((0d0,-1d0)*det))
  if(arg.GT.0d0)then
     if(dble(det).LT.0d0)then
        arg = pi - arg
     end if
  else if(arg.LT.0d0)then
     if(dble(det).LT.0d0)then
        arg = -pi - arg
     end if
  end if

  arg = arg / dble(nmat)
  phase = dcmplx(dcos(arg)) - (0d0,1d0) * dcmplx(dsin(arg))
  u_temp = u_temp * phase
  u_temp_2=u_temp
  call MATRIX_DETERMINANT_COMPLEX(nmat,u_temp_2,det)
  !** check **
  !write(*,*)"detU = ",det
  
  return

END subroutine SU_N_Projection_Each


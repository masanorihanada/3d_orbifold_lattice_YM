!#####################################£##############################
!######              3d YM, orbifold action, HMC            #########
!######                                                     #########
!######                 written by Masanori Hanada          #########
!######                                                     #########
!######         Note: g^2 = 1, instead of g^2 * N = 1.      #########
!####################################################################
!Mersenne twister.
include 'mt19937.f90'
program D3YM

  use mtmod !Mersenne twistor
  implicit none
  include 'size.h'
  include 'include.h'

  double precision pol_phase(1:nmat,1:nx,1:ny),Pol_re,Pol_im
  
  open(unit=10,status='OLD',file='input_v1.dat',action='READ')
  read(10,*) input_config
  read(10,*) output_config
  read(10,*) data_output
  read(10,*) init
  read(10,*) at
  read(10,*) as
  read(10,*) mass2 ! m^2
  read(10,*) mass2_U1  ! (m_{U(1)})^2
  read(10,*) ntraj
  read(10,*) nskip
  read(10,*) ntau
  read(10,*) dtau_t
  read(10,*) dtau_s
  read(10,*) mersenne_seed

  close(10)

  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  if(init.EQ.0)then
     !continue from old config
     open(unit=9,status='OLD',file=input_config,action='READ')
     call mtgetu(9)
     read(9,*) itraj
     read(9,*) umat
     read(9,*) zmat
     close(9)
     
  else if(init.EQ.1)then
     !initialize random number generator
     call sgrnd(mersenne_seed)
     !new config, cold start
     itraj=1
     umat=(0d0,0d0)
     zmat=(0d0,0d0)
     do it=0,nt+1
        do ix=0,nx+1
           do iy=0,ny+1
              do imat=1,nmat
                 umat(imat,imat,it,ix,iy)=(1d0,0d0)
                 do idim=1,ndim-1
                    zmat(imat,imat,it,ix,iy,idim)&
                         =(1d0,0d0)*dcmplx(dsqrt(as**(dble(ndim-3))/2))
                 end do
              end do
           end do
        end do
     end do
  end if
  !**************************************************
  !**************************************************
  nacceptance=0 !number of acceptance
  ntrial=0 !number of trial
  !************************************
  !************************************
  !     Make the output file
  !************************************
  !************************************
  
  open(unit=10,status='REPLACE',file=data_output,action='WRITE')
  write(10,*) "#size of the gauge group: nmat=",nmat
  write(10,*) "#lattice spacings, a_t, a_s=", at, as
  write(10,*) "#ntau=",ntau
  write(10,*) "#dtau for U_t=",Dtau_t
  write(10,*) "#dtau for Z_x,Z_y=",Dtau_s
  write(10,*) "#deformation parameters: m^2, m^2_{U(1)}",mass2,mass2_U1
  write(10,*) "# traj, ham_fin - ham-fin, spatial plaquette from Z, spatial plaquette from U, temporal plaquette from U, Tr(W-1)^2, Re(det U) from spatial links,  Re(Pol.), Im(Pol.), acceptance"
  
  write(10,*) "#------------------------------------------------"

  nacceptance=0
  ntrial=0
  do while (itraj.LE.ntraj)
     backup_umat=umat
     backup_zmat=zmat
     
     call Generate_P(P_umat,P_zmat)
     call Calc_Ham(at,as,mass2,mass2_U1,umat,zmat,P_umat,P_zmat,ham_init)
     call Molecular_Dynamics(at,as,mass2,mass2_U1,ntau,dtau_t,dtau_s,umat,zmat,P_umat,P_zmat)
     call Calc_Ham(at,as,mass2,mass2_U1,umat,zmat,P_umat,P_zmat,ham_fin)
    
     metropolis=grnd()
     ntrial=ntrial+1
     If(dexp(ham_init-ham_fin) > metropolis)THEN
        !accept
        nacceptance=nacceptance+1
     else
        !reject
        umat=backup_umat
        zmat=backup_zmat
     end If

     ! measurements
     if(MOD(itraj,nskip).EQ.0)then
        call measurements(as,umat,zmat,plaq_U,plaq_Z,plaq_temp_U,Wm1,av_det_U)
        call Calc_Polyakov(umat,Pol_re,Pol_im)
        
        write(10,*)itraj,-ham_init+ham_fin,plaq_Z,plaq_U,plaq_temp_U,Wm1,dble(av_det_U),&
             &Pol_re,Pol_im,dble(nacceptance)/dble(ntrial)
        write(*,*)itraj,-ham_init+ham_fin,plaq_Z,plaq_U,plaq_temp_U,Wm1,dble(av_det_U),&
             &Pol_re,Pol_im,dble(nacceptance)/dble(ntrial)
        flush(10)
     end if

     itraj=itraj+1
  end do
  !**************************************************
  !**************************************************
  !   End of iteration
  !**************************************************
  !**************************************************

  close(10)
  close(20)

  open(UNIT = 22, File = output_config, STATUS = "REPLACE", ACTION = "WRITE")
  call mtsaveu(22)
  write(22,*) itraj
  write(22,*) umat
  write(22,*) zmat
  close(22)

end program D3YM

include 'BoxMuller.f90'
include 'MATRIX_iEXP.f90'
include 'Generate_P.f90'
include 'Calc_Ham.f90'
include 'Calc_action.f90'
include 'Molecular_Dynamics.f90'
include 'set_boundary_condition.f90'
include 'Calc_DELH.f90'
include 'Calc_Polyakov.f90'
include 'Calc_measurement.f90'
include 'MATRIX_DET_COMPLEX.f90'
include 'MatrixInverse.f90'
include 'MATRIX_SQRT.f90'

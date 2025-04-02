
  integer nbc !boundary condition for fermions; 0 -> pbc, 1 -> apbc
  integer init !initial condition; 0 -> continue, 1 -> new

!matrices
				     
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)

  double complex backup_umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1)
  double complex backup_zmat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim-1)
				     
  double complex P_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny)
  double complex P_zmat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim-1)
				     
  !matrix indices
  integer imat,jmat,kmat,lmat
  integer it,ix,iy,jt,jx,jy
  integer idim,jdim,kdim,ldim

!lattice spacing along temporal and spatial directions
				     
  double precision at,as

!deformation parameter to remove scalar and U(1)
				     
  double precision mass2,mass2_U1

				     !parameters for molecular evolution
				     
				     
  integer ntau
  double precision dtau_t,dtau_s

  !number of trajectories
				     
  integer ntraj     !total number of trajectories at the end of the run
  integer itraj

  double precision ham_init,ham_fin,metropolis

  !measurements
				     
  integer nskip !measurement is performed every nskiptrajectories
  integer nacceptance !number of acceptance
  integer ntrial !number of trial

  !absolute value of Polyakov loop
				     
  double precision Pol_local(1:nx,1:ny),Pol

  !energy, action, etc
				     
				     double precision energy,action, trx2, trf2,trace, plaq_U, plaq_Z, Wm1
				     double complex av_det_U
character(150) input_config,data_output,output_config,output_pol

!For Mersenne Twister
integer mersenne_seed

				     

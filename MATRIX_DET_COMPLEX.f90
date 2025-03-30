
!***********************************************************
!***********************************************************
!   Sqrt of a matrix : P -> sqrt(P)
SUBROUTINE MATRIX_DETERMINANT_COMPLEX(MatSize,MAT,det)

  implicit none 

  integer MatSize
  double complex MAT(1:MatSize,1:MatSize),MATSQRT(1:MatSize,1:MatSize)

  integer i
  double precision eig_abs(1:MatSize)
  double complex eig_phase(1:MatSize)
  double precision log_det_abs
  double complex det,det_phase
  
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork
  double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       work(1:2*MatSize),rwork(1:2*MatSize)

  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  
  do i=1,MatSize
      eig_abs(i)=abs(w(i))
      eig_phase(i) = w(i)/dcmplx(abs(w(i)))
  end do

  log_det_abs=0d0
  det_phase = (1d0,0d0)
  do i=1,MatSize
     log_det_abs = log_det_abs + dlog(eig_abs(i))
     det_phase = det_phase * eig_phase(i)
  end do

  det = dcmplx(dexp(log_det_abs)) * det_phase
  
  return

END SUBROUTINE MATRIX_DETERMINANT_COMPLEX

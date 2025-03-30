!##############################################################
!######    Estimation of auto-correlation                 #####
!######    written by Masanori Hanada                     #####     
!##############################################################
program JackKnife

  implicit none
  !----------------------------------------
  !     Iteration
  integer iteration
  parameter(iteration=500)
  !----------------------------------------
  ! skip unthermalized configs
  !integer ntherm
  !parameter(ntherm=0)
  !---------------------------------------- 
  integer Binsize,BinStart,BinEnd,BinStep
  parameter(BinStart=2)
  parameter(BinEnd=51)
  parameter(BinStep=2)
  !----------------------------------------
  character(150) input,output,aho
  integer i,j

!  doubleprecision readdata(1:6),average(1:6),err(1:6),BinAverage(1:6)
 doubleprecision readdata(1:7),average(1:7),err(1:7),BinAverage(1:7)
  integer read1
  !*************************
  !*** evaluation of err *** 
  !*************************
  write(input,"('gauged-N8L16T06.txt')")
  BinSize=BinStart
  do while(Binsize.LT.BinEnd)
     
     open(UNIT = 12, File = input, STATUS = "OLD", ACTION = "READ")
     do i=1,6
        read(12,*)aho
     end do

     average=0d0
     err=0d0
     do i=1,INT(dble(iteration)/dble(BinSize))
        BinAverage=0d0
        do j=1,Binsize
           read(12,*)readdata
           BinAverage=BinAverage+readdata
        end do
        BinAverage=BinAverage/dble(Binsize)
        average=average+BinAverage
        err=err+BinAverage**2d0

     end do
     close(12)
     average=average/dble(INT(dble(iteration)/dble(BinSize)))
     err=err/dble(INT(dble(iteration)/dble(BinSize)))
     err=err-average**2d0
     err=sqrt(err/dble(INT(dble(iteration)/dble(BinSize))-1))
     write(*,*)dble(binsize),&
          average(3),err(3),average(5),err(5),average(6),err(6)
     Binsize=Binsize+BinStep
  end do
  




 


end program JackKnife

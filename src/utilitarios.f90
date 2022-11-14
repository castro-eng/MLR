!=================================================================================
!
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Garcia,        bidu@lncc.br [4]
!         Tuane Lopes,           tuane@lncc.br [5]
!
!         LNCC/MCT
!         Petropolis, 02.2014
!=================================================================================
!
!
      subroutine ztest(a,n,lzero)
      use mGlobaisEscalares, only : zero
!
!.... program to determine if an array contains only zero entries
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer*4, intent(in)    :: n 
      real*8, intent(in)     :: a(n)
      logical, intent(inout) :: lzero
!
      integer*4:: i
!
      lzero = .true.
!
      do 100 i=1,n
      if (a(i).ne.zero) then
         lzero = .false.
         return
      endif
  100 continue
!
      end     
!
!**** new **********************************************************************
!
      subroutine timing(tempoDeParede)
      implicit none     
!
!.... program to determine elapsed cpu time
!
!
      real*8, intent(inout) :: tempoDeParede
!
      character(LEN=8)  :: date
      character(LEN=10) :: minhaHora
      character(LEN=5)  :: zone
      integer*4,dimension(8) :: values
      integer*4::  horas, minutos, segundos, milesimosSeg
!
      call date_and_time(date,minhaHora,zone,values);
!
      horas=values(5);      minutos=values(6); 
      segundos=values(7); milesimosSeg=values(8);    
!
      tempoDeParede = (60*horas+minutos)*60.0+segundos+milesimosSeg/1000.00
!
      return
      end
!
!:::: new ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine imove(ia,ib,n) 
! 
!.... program to move an integer*4array 
! 
      integer*4:: ia(1),ib(*) 
      integer*4:: n
!
      integer*4:: i
! 
      do 100 i=1,n 
      ia(i)=ib(i) 
  100 continue 
 
      return 
      end subroutine
!:::: new :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
      subroutine move(a,b,n)
! 
!.... program to move a floating-point array 
! 
      implicit real*8(a-h,o-z) 
! 
!.... remove above card for single-precision operation 
! 
      real*8 :: a(*),b(*) 
      integer*4:: n
!
      integer*4 i
! 
      do 100 i=1,n 
      a(i) = b(i) 
  100 continue 
! 
      return 
      end subroutine
!
!**** new **********************************************************************
!
      subroutine gerarLabel(label,tempo)
!
      implicit none


      character(LEN=21), intent(out) :: label
      real*8, intent(in) :: tempo
!
      character(LEN=21) :: labelAux, num     
      integer :: i
            
      write(num,'(f20.4)') tempo
!     labelAux="t="//ADJUSTL(num)
      labelAux="t="//num
      do i = 1, 21
         if(labelAux(i:i) .ne. ' ') then
            label(i:i) = labelAux(i:i)
         else
            label(i:i) = '0'
         end if
      end do
     
      end subroutine                                                     

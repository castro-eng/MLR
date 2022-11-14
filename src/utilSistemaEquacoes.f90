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
module mUtilSistemaEquacoes

   use mEstruturasDadosSistEq

   ! funcoes e subrotinas
   public :: load, dirichletConditions
   public :: btod, kdbc, pivots, colht, coldot, matadd
   public :: montarEstrutDadosSistEqAlq

   contains

   !**** new **********************************************************************
   subroutine montarEstrutDadosSistEqAlq(umaEstSistEq_)
      use mMalha,            only: conecNodaisElem, numnp, numel, nen, nsd
      use mMalha,            only: conecNodaisFratura, numelFratura, nenFratura
      use mMalha,            only: numelBordas, conecNodaisBordas, numelDescontNaBorda, conecNodaisDesconNaBorda
      
      use mMalha,            only: listaDosElemsPorNoCSR, criarListaVizinhosCRS
      use mLeituraEscrita,   only: iecho
      use mSolverHypre

      implicit none
      type(estruturasArmazenamentoSistemaEq) :: umaEstSistEq_

      logical   :: simetria
      integer*4 :: meanbw
      real*8    :: t1, t2, t3
      integer*4 :: i, localSize

      write(*,'(a,i2,a)',advance='NO') " em montarEstrutDadosSistEqAlq, "
    
      call timing(t1)
      allocate(umaEstSistEq_%lm(umaEstSistEq_%ndof,nen,numel))
      call formlm(umaEstSistEq_%id,conecNodaisElem,umaEstSistEq_%lm,umaEstSistEq_%ndof,umaEstSistEq_%ndof,nen,numel)
      allocate(umaEstSistEq_%idiag(umaEstSistEq_%neq+1)); umaEstSistEq_%idiag = 0  ! posicoes dos elementos da diagonal principal no vetor alhs    
      
      allocate(umaEstSistEq_%lmFratura(umaEstSistEq_%ndof,nenFratura,numelFratura))
      call formlm(umaEstSistEq_%id,conecNodaisFratura,umaEstSistEq_%lmFratura,umaEstSistEq_%ndof,&
                            umaEstSistEq_%ndof,nenFratura,numelFratura)
      
      allocate(umaEstSistEq_%lmBorda(umaEstSistEq_%ndof,nenFratura,numelBordas))
      call formlm(umaEstSistEq_%id,conecNodaisBordas,umaEstSistEq_%lmBorda,umaEstSistEq_%ndof,&
                  umaEstSistEq_%ndof,nenFratura,numelBordas)
      if ( numelFratura > 0 ) then
         allocate(umaEstSistEq_%lmFraturaNaBorda(umaEstSistEq_%ndof,(nenFratura-1),numelDescontNaBorda))      
         call formlm(umaEstSistEq_%id,conecNodaisDesconNaBorda(:,:),umaEstSistEq_%lmFraturaNaBorda,umaEstSistEq_%ndof,&
                     umaEstSistEq_%ndof,nenFratura-1,numelDescontNaBorda)                            
      end if
            
      if(umaEstSistEq_%optSolver=="skyline") then
         call colht(umaEstSistEq_%idiag,umaEstSistEq_%lm,umaEstSistEq_%ndof,nen,numel,umaEstSistEq_%neq);
         call diag (umaEstSistEq_%idiag,umaEstSistEq_%neq,umaEstSistEq_%nalhs); 
         if(.not.associated(umaEstSistEq_%alhs))allocate(umaEstSistEq_%alhs(umaEstSistEq_%nalhs));
         umaEstSistEq_%alhs=0.0 
      end if 
      
      if(umaEstSistEq_%optSolver=="pardiso") then

         ! if(.not.associated (listaDosElemsPorNoCSR) ) then
         !    allocate(listaDosElemsPorNoCSR(nen,numnp)); listaDosElemsPorNoCSR=0
         ! end if

         ! call criarListaVizinhos(nen, numnp, numel, conecNodaisElem, listaDosElemsPorNoCSR)
         call criarListaVizinhosCRS(nen,numnp,numel,conecNodaisElem,umaEstSistEq_%nVizinMax)
            
         simetria=.true.
         umaEstSistEq_%Ap=> umaEstSistEq_%idiag; !posicoes dos elementos da diagonal principal no vetor alhs
            
            
         call criarPonteirosMatEsparsa_CSR(umaEstSistEq_, nsd, conecNodaisElem, listaDosElemsPorNoCSR, & 
                  numnp, nen, umaEstSistEq_%nVizinMax,  simetria)       
                                                       
         ! if(.not.associated(umaEstSistEq_%alhs)) 
         allocate(umaEstSistEq_%alhs(umaEstSistEq_%nalhs)); umaEstSistEq_%alhs=0.0
      endif 

      if(umaEstSistEq_%optSolver=="hypre") then
         
         Clower = 1 - 1
         Cupper = umaEstSistEq_%neq-1
         localSize=CUpper-Clower+1
         if(.not.allocated(rows)) allocate(rows(localSize))
         do i = 1, localSize
            rows(i) = i - 1
         end do
         call criarMatriz_HYPRE    (A_HYPRE, Clower, Cupper, mpi_comm )
         call criarVetor_HYPRE     (b_HYPRE, Clower, Cupper, mpi_comm )
         call criarVetor_HYPRE     (u_HYPRE, Clower, Cupper, mpi_comm )

      endif

      if(.not.associated(umaEstSistEq_%brhs)) allocate(umaEstSistEq_%brhs(umaEstSistEq_%neq));
      umaEstSistEq_%brhs=0.0

      call timing(t2)
      write(*,*) "criar estruturas de dados,  tempo de parede = ", t2 - t1
      
      meanbw = umaEstSistEq_%nalhs/umaEstSistEq_%neq
      write(*,    6000) 'Informacao do sistema de equacoes ',  &
               umaEstSistEq_%neq, umaEstSistEq_%nalhs, meanbw, (8.0*umaEstSistEq_%nalhs)/1000.0/1000.0, &
               umaEstSistEq_%optSolver
      write(iecho,6000) 'Informacao do sistema de equacoes ',  &
               umaEstSistEq_%neq, umaEstSistEq_%nalhs, meanbw, (8.0*umaEstSistEq_%nalhs)/1000.0/1000.0, &
               umaEstSistEq_%optSolver

         
      6000 format(a///&
     ' e q u a t i o n    s y s t e m    d a t a              ',  //5x,&
     ' number of equations . . . . . . . . . . . . (neq    ) = ',i8//5x,&
     ' number of terms in left-hand-side matrix  . (nalhs  ) = ',i12//5x,&
     ' mean half bandwidth . . . . . . . . . . . . (meanbw ) = ',i8//5x,&
     ' memoria necessaria para a matriz do sistema (Mbytes)  = ',e10.2//5x, &
     ' Solver escolhido                                      = ',a)

   end subroutine montarEstrutDadosSistEqAlq

   !**** new **********************************************************************
   subroutine load(id,f,brhs,ndof,numnp,nlvect)
      !.... program to accumulate nodal forces and transfer into
      !        right-hand-side vector

      implicit none

      integer*4 :: ndof, numnp, nlvect
      integer*4:: id(ndof,*)
      real*8  :: f(ndof,numnp,*),brhs(*)

      integer*4:: nlv
      integer*4:: i, j, k

      do i=1,ndof

         do j=1,numnp
            k = id(i,j)
            if (k.gt.0) then

               do nlv=1,nlvect
                  brhs(k) = brhs(k) + f(i,j,nlv)
               enddo
            endif
         enddo
      enddo

      do  i=1,ndof  ! BD
         do  j=1,numnp
            k = id(i,j)
            ! write(*,*) 'j=',j,', id=', id(i, j), ', brhs=',  brhs (k)
         end do
      end do

      return
   end subroutine

   !**** new **********************************************************************
   subroutine dirichletConditions(id,d,f,ndof,numnp,nlvect)
      !.... program to compute displacement boundary conditions

      implicit none

      integer*4:: ndof, numnp, nlvect
      integer*4:: id(ndof,*)
      real*8  :: d(ndof,*),f(ndof,numnp,*)

      integer*4:: i, j, k, lv
      real*8  :: val

      do 300 i=1,ndof
         do 200 j=1,numnp
            k = id(i,j)

            if (k.gt.0) go to 200
            val = 0.d0
            do 100 lv=1,nlvect
               val = val + f(i,j,lv)
            100 continue

            d(i,j) = val
         200 continue
      300 continue
      return
   end subroutine

   !**** new **********************************************************************
   subroutine btod(id,d,brhs,ndof,numnp)
      !.... program to perform transfer from r.h.s. to displacement array

      implicit none

      integer*4:: ndof, numnp
      integer*4:: id(ndof,*)
      real*8  :: d(ndof,*),brhs(*)

      integer*4:: i, j, k

      do i=1,ndof

         do j=1,numnp
            k = id(i,j)
            if (k.gt.0) then 
               d(i,j) = brhs(k)
            end if
         enddo
      enddo

      return
   end subroutine

   !**** new **********************************************************************
   subroutine kdbc(eleffm,elresf,dl,nee)
      !.... program to adjust load vector for prescribed displacement
      !     boundary condition

      implicit none

      integer*4:: nee
      real*8  :: eleffm(nee,*),elresf(*),dl(*)

      integer*4:: i,j
      real*8  :: val

      do 200 j=1,nee
         val=dl(j)
      
         if(val.eq.0.0d0) go to 200

         do 100 i=1,nee
            elresf(i)=elresf(i)-eleffm(i,j)*val
         100 continue

      200 continue

      return
   end subroutine

   !**** NEW **********************************************************************
   SUBROUTINE KDBC2(ELEFFM,ELRESF,DL,NEE,LMT,NEL)
      !.... PROGRAM TO ADJUST LOAD VECTOR FOR PRESCRIBED DISPLACEMENT
      !     BOUNDARY CONDITION
      use mGlobaisEscalares, only: zero

      IMPLICIT NONE

      real*8 :: ELEFFM(NEE,*),ELRESF(*),DL(*),VAL
      integer :: nel, nee
      integer*4, pointer :: lmT(:,:,:)

      integer :: i, j, l
      integer*4:: lm(nee)
      
      lm(:) =reshape(lmT(:,:,nel),(/nee/)); 

      DO 200 J=1,NEE

         L=LM(J) 
         VAL=DL(J)
   
         IF(L.GT.0) GO TO 200
         IF(VAL.EQ.ZERO) GO TO 200

         DO 100 I=1,NEE
            ELRESF(I)=ELRESF(I)-ELEFFM(I,J)*VAL
         100 CONTINUE

      200 CONTINUE

      RETURN
   END SUBROUTINE

   !**** new **********************************************************************
   subroutine pivots(a,idiag,neq,nsq,iecho,*)
      !.... program to determine the number of zero and negative terms in
      !        array d of factorization a = u(transpose) * d * u

      implicit none

      real*8  :: a(*)                                                       
      integer*4:: idiag(*)
      integer*4:: neq 
      integer*4:: nsq, iecho

      integer*4:: iz, in, n, i

      iz = 0
      in = 0

      do n=1,neq
         i = idiag(n)
         if (a(i).eq.0.) iz = iz + 1
         if (a(i).lt.0.) in = in + 1
      enddo

      write(iecho,1000) nsq,iz,in

      return 1

      1000 format(' ',&
               ' zero and/or negative pivots encountered                ', ///5x,&
               ' time sequence number   . . . . . . . . . . . (nsq  ) = ',i10//5x,&
               ' number of zeroes . . . . . . . . . . . . . . . . . . = ',i10//5x,&
               ' number of negatives  . . . . . . . . . . . . . . . . = ',i10//5x)

   end subroutine

   !**** new **********************************************************************
   subroutine btdb(elstif,b,db,nee,nrowb,nstr)
      !.... program to multiply b(transpose) * db taking account of symmetry
      !        and accumulate into element stiffness matrix

      implicit none

      integer*4:: nee,nrowb
      real*8  :: elstif(nee,*),b(nrowb,*),db(nrowb,*)
      integer*4:: nstr

      integer*4:: i,j

      do j=1,nee
         do i=1,j
            elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)
         enddo
      enddo

      return
   end subroutine

   !**** new **********************************************************************
   function coldot(a,b,n)
      !.... program to compute the dot product of vectors stored column-wise

      implicit none

      integer*4:: n
      real*8  :: a(n),b(n)

      real*8  :: coldot
      integer*4:: i

      real*8  :: dot_product
      real*8  :: ddot
    
      coldot = 0.0d0

      do i=1,n
         coldot = coldot + a(i)*b(i)
      enddo

      ! coldot = dot_product(a,b) ! intrinsec fortran funtion
      ! coldot = ddot(n,a,1,b,1)    ! external blas function 

      return
   end function

   !**** new **********************************************************************
   subroutine matadd(a,b,c,ma,mb,mc,m,n,iopt)
      !.... program to add rectangular matrices

      implicit none

      integer*4:: ma,mb,mc,m,n,iopt
      real*8  :: a(ma,*),b(mb,*),c(mc,*)

      integer*4:: i,j

      go to (1000,2000,3000),iopt

      !.... iopt = 1, add entire matrices
      1000 do j=1,n
         do i=1,m 
            c(i,j) = a(i,j) + b(i,j)
         enddo
      enddo
      return

      !.... iopt = 2, add lower triangular and diagonal elements
      2000 do j=1,n
         do i=j,m 
            c(i,j) = a(i,j) + b(i,j)
         enddo
      enddo  
      return

      !.... iopt = 3, add upper triangular and diagonal elements
      3000 do j=1,n
         do i=1,j 
            c(i,j) = a(i,j) + b(i,j)
         enddo
      enddo
      return

   end subroutine

end module

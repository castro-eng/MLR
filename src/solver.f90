!=================================================================================
!
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Garcia,        bidu@lncc.br [4]
!         Tuane Lopes,           tuane@lncc.br [5]
!         Eduardo Castro,        ecastro@lncc.br
!
!         LNCC/MCT
!         Petropolis, 19.2020
!=================================================================================
module mSolver
   ! Modulo contendo variaveis e operações referentes ao fluxo

   implicit none

   logical :: listaSolverDisponivel(3)
   logical :: simetria


   contains

   !==========================================================================================================    
   subroutine identificaSolversDisponiveis()
      ! Subrotina para identificação dos solveres disponíveis

      ! skyline
      listaSolverDisponivel(1)=.true.
      ! pardiso
#ifdef withPardiso
      listaSolverDisponivel(2)=.true.
#else
      listaSolverDisponivel(2)=.false.
#endif
      ! hypre
#ifdef withHYPRE
      listaSolverDisponivel(3)=.true.
#else
      listaSolverDisponivel(3)=.false.
#endif
   
      print*, "Solvers disponiveis:"
      print*, "                      SKYLINE"
      if(listaSolverDisponivel(2).eqv..true.) print*, "                      PARDISO"
      if(listaSolverDisponivel(3).eqv..true.) print*, "                      HYPRE"
   end subroutine

   !==========================================================================================================    
   subroutine verificarSolver(optSolver_)
      ! Verificação se o solver escolhido se encontra disponível

      character(len=7), intent(in) :: optSolver_

      if(optSolver_=='pardiso')then
         if (listaSolverDisponivel(2).eqv..false.) then
            print*, "O Solver escolhido ...,  ", optSolver_,",  não está disponível"
            stop 100
         endif
      endif

      if(optSolver_=='hypre') then
         if(listaSolverDisponivel(3).eqv..false.) then
            print*, "O Solver escolhido ...,  ", optSolver_, ",  não está disponível"
            stop 100
         endif
      endif
   end subroutine   

   !==========================================================================================================    
   subroutine solver(optSolver_ , estrutSistEq_, label_)         
      use mSolverGaussSkyline,      only: solverGaussSkyline, FACTOR, BACK
      use mSolverPardiso,           only: solverPardisoEsparso
      use mSolverHypre
      use mUtilSistemaEquacoes,     only: btod
      use mMalha,                   only: x, numnp
      use mEstruturasDadosSistEq,   only: estruturasArmazenamentoSistemaEq
      use mGlobaisEscalares,        only: passoTempo

      implicit none
      character(len=*), intent(in) :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*), intent(in) :: label_
      
      character(LEN=12) :: etapa
      integer * 4 ::  num_iterations, printSol = 2
      real*8 ::  final_res_norm, elapsedT, tol, solutionNorm
      integer*4   :: i, equacao, dof

      if (optSolver_=='skyline') then
         IF(passoTempo==1)call factor(estrutSistEq_%alhs,estrutSistEq_%idiag,estrutSistEq_%nalhs,estrutSistEq_%neq)            
         call back  (estrutSistEq_%alhs,estrutSistEq_%brhs,estrutSistEq_%idiag,estrutSistEq_%neq)
      endif
          
      if (optSolver_=='pardiso') then
         simetria=.true.
         ! if(PASSOTEMPO==1) THEN
         !    etapa='full';
         ! else
         !    etapa='back';
         ! endif
         etapa='full';
         call solverPardisoEsparso(estrutSistEq_, simetria, etapa, label_)
      endif
          
      if(optSolver_=='hypre') then
         write(*,'(2a)') ' iterativo ', optSolver_
            
         do i = 1, estrutSistEq_%neq
            rows(i) = i-1 
         end do
         call atribuirValoresVetor_HYPRE(b_HYPRE, 1, estrutSistEq_%neq, rows, estrutSistEq_%brhs)

         if (passoTempo==1) call fecharMatriz_HYPRE         (A_HYPRE, parcsr_A)
         call fecharVetor_HYPRE          (b_HYPRE, par_b   )
         call fecharVetor_HYPRE          (u_HYPRE, par_u   )

         if(.not.allocated(initialGuess)) then
            allocate(initialGuess(estrutSistEq_%neq)); initialGuess=0.0
         endif
            
         solver_id  = 1;
         precond_id = 2; tol = 1.0e-10
         write(*,*) "em hidro, precond_id = ", precond_id, ", tol = ", tol 
         call resolverSistemaAlgHYPRE (A_HYPRE, parcsr_A, b_HYPRE, par_b, u_HYPRE, par_u, &
                  solverH, solver_id, precond_id, tol,num_iterations, final_res_norm, initialGuess, &
                  estrutSistEq_%brhs, rows, estrutSistEq_%neq, myid, mpi_comm)

         call extrairValoresVetor_HYPRE(u_HYPRE, 1, estrutSistEq_%neq, rows,estrutSistEq_%BRHS)
                  initialGuess=estrutSistEq_%brhs

         ! call destruirMatriz_HYPRE(A_HYPRE)
         call destruirVetor_HYPRE (b_HYPRE)
         call destruirVetor_HYPRE (u_HYPRE)

         ! call criarMatriz_HYPRE  (A_HYPRE, Clower, Cupper, mpi_comm )
         call criarVetor_HYPRE   (b_HYPRE, Clower, Cupper, mpi_comm )
         call criarVetor_HYPRE   (u_HYPRE, Clower, Cupper, mpi_comm )

      endif

      solutionNorm = 0.0
      do i = 1, estrutSistEq_%neq
         solutionNorm = solutionNorm + estrutSistEq_%brhs(i)**2
      end do
      solutionNorm = sqrt(solutionNorm)

      call btod(estrutSistEq_%id,estrutSistEq_%u, estrutSistEq_%brhs,estrutSistEq_%ndof,numnp)

      ! 200       FORMAT(2(1PE15.8,2X))
   end subroutine solver
   
end module mSolver




!=================================================================================
!                    Modulo Gauss Skyline
!=================================================================================
module mSolverGaussSkyline

   !funcoes e subrotinas
   public :: backns, factns, back, factor
   public :: diag, load, addnsl, addlhs, addrhs
   public :: btod, kdbc, pivots, dirichletConditions, btdb, predct, colht

   public :: coldot
   public :: matadd
   integer*4:: nada

   contains
        
   !**** new **********************************************************************
   subroutine solverGaussSkyline(estrutSistEq_)
      use mEstruturasDadosSistEq 

      implicit none

      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_

      ! write(*,*) " em solverGaussSkyline, call factor"
      call factor(estrutSistEq_%alhs,estrutSistEq_%idiag,estrutSistEq_%nalhs,estrutSistEq_%neq)
      ! write(*,*) " em solverGaussSkyline, call back"
      call back  (estrutSistEq_%alhs,estrutSistEq_%brhs,estrutSistEq_%idiag,estrutSistEq_%neq)

      !write(*,'(5f10.5)') brhs(1:5)
      !write(*,'(5f10.5)') brhs(neq-4: neq)

   end subroutine solverGaussSkyline

   !**** new*********************************************************
   subroutine addnsl(alhs,clhs,eleffm,idiag,lm,nee,ldiag)
      ! program to add element left-hand-side matrix to          
      ! global left-hand-side matrix                      

      !        ldiag = .true.,  add diagonal element matrix              
      !                                                                  
      !        ldiag = .false, then                                     
      !        add full nonsymmetric element matrix                   
                                                                 
      implicit none

      integer*4:: nee
      real*8  :: alhs(*),clhs(*),eleffm(nee,*)
      integer*4:: idiag(*),lm(*)
      logical ldiag

      integer*4:: i,j,k,l,m

      if (ldiag) then
         do j=1,nee
             k = abs(lm(j))
            if (k.gt.0) then
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j)
            endif
         enddo

      else
         do j=1,nee
             k = abs(lm(j))
            if (k.gt.0) then
               do i=1,nee
                   m = abs(lm(i))
                  if (m.gt.0) then
                     if (k.gt.m) then
                        l = idiag(k) - k + m
                        alhs(l) = alhs(l) + eleffm(i,j)
                     else
                        l = idiag(m) - m + k
                        clhs(l) = clhs(l) + eleffm(i,j)
                     endif
                     if (k.eq.m) then
                        l = idiag(k)
                        alhs(l) = alhs(l) + eleffm(i,j)
                        clhs(l) = alhs(l)
                     endif
                  endif
               enddo
            endif
         enddo
      endif

      return
   end subroutine

   !**** new*********************************************************
   subroutine backns(a,c,b,idiag,neq)
      !.... program to perform forward reduction and back substitution

      implicit none

      !.... deactivate above card(s) for single-precision operation

      integer*4:: neq
      real*8  :: a(*),c(*),b(*)
      integer*4:: idiag(*)

      integer*4:: i, j, jj, jcolht, jjlast, jjnext, istart, jtemp
      real*8  :: ajj, bj

      !.... forward reduction

      jj = 0

      do j=1,neq
         jjlast = jj
         jj     = idiag(j)
         jcolht = jj - jjlast
         if (jcolht.gt.1) then
            b(j) = b(j) - coldot(c(jjlast+1),b(j-jcolht+1),jcolht-1)
         endif
      enddo

      !.... diagonal scaling
      do j=1,neq
         ajj = a(idiag(j))

         !.... warning: diagonal scaling is not performed if ajj equals zero
         if (ajj.ne.0.d0) b(j) = b(j)/ajj
      enddo
      
      !.... back substitution
      if (neq.eq.1) return
      jjnext = idiag(neq)

      do j=neq,2,-1
         jj     = jjnext
         jjnext = idiag(j-1)
         jcolht = jj - jjnext
         if (jcolht.gt.1) then
            bj = b(j)
            istart = j - jcolht + 1
            jtemp  = jjnext - istart + 1

            do i=istart,j-1
               b(i) = b(i) - a(jtemp+i)*bj
            enddo
         endif
      enddo

      return
   end subroutine

   !**** new *********************************************************
   subroutine factns(a,c,idiag,neq)
      !.... program to perform crout factorization: a = l * d * u

      !        a(i):  coefficient matrix stored in compacted column form;
      !               after factorization contains d and u
      
      !        c(i):  non-symmetric lower triangular coefficient matrix stored in
      !                compacted row form; after factorization contains l

      implicit none

      !.... deactivate above card(s) for single-precision operation

      real*8  :: a(*),c(*)
      integer*4:: idiag(*)
      integer*4:: neq

      integer*4:: i, j, ii, jj, ij, iilast, jjlast
      integer*4:: istart,  icolht, jcolht, jtemp, jm1
      integer*4:: length

      ! write(*,*) ' em factorns '

      jj = 0

      do j=1,neq

         jjlast = jj
         jj     = idiag(j)
         jcolht = jj - jjlast

         if (jcolht.gt.2) then

            !....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
            istart = j - jcolht + 2
            jm1    = j - 1
            ij     = jjlast + 2
            ii     = idiag(istart-1)

            do i=istart,jm1
               iilast = ii
               ii     = idiag(i)
               icolht = ii - iilast
               length = min0(icolht-1,i - istart + 1)
               if (length.gt.0)  then
                  a(ij) = a(ij) - coldot(a(ij-length),c(ii-length),length)
                  c(ij) = c(ij) - coldot(c(ij-length),a(ii-length),length)
               endif
               ij = ij + 1
            enddo
         endif

         if (jcolht.ge.2) then
            !....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
            !           replace a(j,j) with d(j,j).

            jtemp = j - jj
            do ij=jjlast+1,jj-1
               ii = idiag(jtemp + ij)

               !....... warning: the following calculations are skipped 
               !                 if a(ii) equals zero
               if (a(ii).ne.0.0d0) then
                  c(ij) = c(ij)/a(ii)
                  a(jj) = a(jj) - c(ij)*a(ij)
                  a(ij) = a(ij)/a(ii)
               endif
            enddo
         endif
      enddo

      return
   end subroutine

   !**** NEW ****************************************************************** 
   subroutine back(a,b,idiag,neq)
      !.... program to perform forward reduction and back substitution 

      implicit none

      real*8  :: a(*)
      integer*4:: neq
      integer*4:: idiag(neq)
      real*8  :: b(neq)

      integer*4:: i, j, jcolht, istart, jtemp
      integer*4:: jj, jjnext
      integer*4:: jjlast
      real*8  :: ajj, bj

      !.... forward reduction 

      jj = 0

      do j=1,neq
         jjlast = jj
         jj     = idiag(j)
         jcolht = jj - jjlast
         if (jcolht.gt.1) then
            b(j) = b(j) - coldot(a(jjlast+1),b(j-jcolht+1),jcolht-1) 
         end if
      enddo

      !.... diagonal scaling 

      do j=1,neq
         ajj = a(idiag(j))
         if (ajj.ne.0.d0) then
            b(j) = b(j)/ajj
         end if
      end do

      !.... back substitution 

      if (neq.eq.1) return
      jjnext = idiag(neq)

      do j=neq,2,-1
         jj     = jjnext
         jjnext = idiag(j-1)
         jcolht = jj - jjnext
         if (jcolht.gt.1) then
            bj = b(j)
            istart = j - jcolht + 1
            jtemp  = jjnext - istart + 1
  
            do i=istart,j-1
               b(i) = b(i) - a(jtemp+i)*bj
            enddo
         endif
      end do

      return
   end subroutine
   
   !**** new *************************************************************** 
   subroutine factor(a,idiag,nalhs,neq)
      !.... program to perform crout factorization: a = u(transpose) * d * u 

      !        a(i):  coefficient matrix stored in compacted column form; 
      !               after factorization contains d and u 

      implicit none

      integer*4, intent(in) :: neq, nalhs
      real*8, intent(inout)  :: a(nalhs)

      integer*4, intent(in) :: idiag(neq)

      integer*4:: i, j, jlast, icolht, jcolht, istart, jm1, jtemp
      integer*4:: ii, ij, jj, jlngth, length, iilast, jjnext
      integer*4:: jjlast
      real*8  :: ajj, bj, temp 
      ! real*8   ::  dot_product
      real*8, external   ::  ddot

      ! write(*,*) ' em factor '
      ! write(*,*) ' a(1:nalhs) =', a(1:nalhs)
      jj = 0
      i=0; j=0; jlast=0; icolht=0; jcolht=0; istart=0; jm1=0; jtemp=0;
      ii=0;ij=0;jj=0;jlngth=0;length=0;iilast=0;jjnext=0;jjlast=0;
      ajj=0.0;bj=0.0;temp=0.0;

      do j=1,neq

         jjlast = jj
         jj     = idiag(j)
         jcolht = jj - jjlast

         if (jcolht.gt.2) then
            !....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j) 

            istart = j - jcolht + 2
            jm1    = j - 1
            ij     = jjlast + 2
            ii     = idiag(istart-1)

            do i=istart,jm1

               iilast = ii
               ii     = idiag(i)
               icolht = ii - iilast
               jlngth = i - istart + 1
               length = min0(icolht-1,jlngth)
               if (length.gt.0) then
                  ! a(ij) = a(ij) - coldot(a(ii-length),a(ij-length),length)
                  ! a(ij) = a(ij) - ddot(length, a(ii-length),1, a(ij-length),1) ! external blas function
                  a(ij) = a(ij) - dot_product(a(ii-length:ii-1),a(ij-length:ij-1))
               end if
               ij = ij + 1
            enddo
         endif

         if (jcolht.ge.2) then
            !....... for column j and i.le.j-1, replace a(i,j) with u(i,j); 
            !           replace a(j,j) with d(j,j). 

            jtemp = j - jj
            do ij=jjlast+1,jj-1

               ii = idiag(jtemp + ij)
               if (a(ii).ne.0.d0) then
                  temp  = a(ij)
                  a(ij) = temp/a(ii)
                  a(jj) = a(jj) - temp*a(ij)
               endif
            enddo
         endif
      enddo

      return
   end subroutine

   ! **** new *********************************************************************
   subroutine addlhsN(alhs,eleffm,lmT, idiag, nee, nel, ldiag,lsym)
      ! .... program to add element left-hand-side matrix to
      !         global left-hand-side matrix

      !         ldiag = .true.,  add diagonal element matrix
      !         ldiag = .false., then
      !          lsym = .true., add upper triangle of full element matrix
      !          lsym = .false., add full element matrix

      implicit none

      ! .... deactivate above card(s) for single-precision operation

      integer*4:: nee, nel
      real*8  :: alhs(*),eleffm(nee,nee)
      integer*4, pointer :: lmT(:,:,:), idiag(:)
      logical :: ldiag,lsym

      integer*4:: i, j, l, k, m
      integer*4:: lm(nee)

      ! nee = size(estrutSistEqP_%lm,1)*size(estrutSistEqP_%lm,2)

      ! write(*,*) " em addlnsN, nel=", nel, ", nee = ", nee

      lm(:) =reshape(lmT(:,:,nel),(/nee/)); 
      ! write(*,*) "lm = ", lm
      ! nel = nel + 1
      if (ldiag) then
         do j=1,nee
            k = lm(j)
            if (k.gt.0) then
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j)
            endif
         enddo

      else
         do j=1,nee
            k = lm(j)
            if (k.gt.0) then
               do i=1,j
                  m = lm(i)
                  if (m.gt.0) then
                     if (k.ge.m) then
                        l = idiag(k) - k + m
                     else
                        l = idiag(m) - m + k
                     endif
                     alhs(l) = alhs(l) + eleffm(i,j)
                  endif
               enddo

               if (.not. lsym) then
                  do i = j,nee
                     m = lm(i)
                     if (m .gt. 0) then
                        if (k .ge. m) then
                           l = idiag(k) - k + m
                        else
                           l = idiag(m) - m + k
                        endif
                     endif
                  enddo
               endif
            endif
         enddo
      endif

      return
   end subroutine

   !**** new **********************************************************************
   subroutine addrhsN (brhs,elresf,lmT,nee, nel)
      !.... program to add element residual-force vector to
      !        global right-hand-side vector

      implicit none

      real*8  :: brhs(*),elresf(*)
      integer*4, pointer :: lmT(:,:,:)
      integer*4:: nee, nel
      integer*4:: k, j
      integer*4 :: lm(nee)
   
      ! write(*,*) " em addrhsN, nel=", nel, ", nee = ", nee
      lm(:) =reshape(lmT(:,:,nel),(/nee/)); 
             
      do j=1,nee
         k = lm(j)
         if (k.gt.0) then 
            brhs(k) = brhs(k) + elresf(j)
            ! write(*,*) k, brhs(k), j, elresf(j)
         end if
      enddo

      return
   end subroutine

   ! **** new *********************************************************************
   subroutine addlhs(alhs,eleffm,idiag,lm,nee,ldiag,lsym)
      ! .... program to add element left-hand-side matrix to
      !         global left-hand-side matrix

      !         ldiag = .true.,  add diagonal element matrix
      !         ldiag = .false., then
      !           lsym = .true., add upper triangle of full element matrix
      !           lsym = .false., add full element matrix

      implicit none

      ! .... deactivate above card(s) for single-precision operation

      integer*4:: nee
      real*8  :: alhs(*),eleffm(nee,nee)
      integer*4:: idiag(*),lm(*)
      logical :: ldiag,lsym

      integer*4:: i, j, l, k, m
      integer*4, save ::  nel = 0
       
      nel = nel + 1
      ! write(*,*) "nel =", nel, " lm = ", lm((nel-1)*nee+1:nel*nee)
      if (ldiag) then
         do j=1,nee
            k = lm(j)
            if (k.gt.0) then
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j)
            endif
         enddo

      else
         do j=1,nee
            k = lm(j)
            if (k.gt.0) then
               do i=1,j
                  m = lm(i)
                  if (m.gt.0) then
                     if (k.ge.m) then
                        l = idiag(k) - k + m
                     else
                        l = idiag(m) - m + k
                     endif
                     alhs(l) = alhs(l) + eleffm(i,j)
                  endif
               enddo

               if (.not. lsym) then
                  do i = j,nee
                     m = lm(i)
                     if (m .gt. 0) then
                        if (k .ge. m) then
                           l = idiag(k) - k + m
                        else
                           l = idiag(m) - m + k
                        endif
                     endif
                  enddo
               endif
            endif
         enddo
      endif

      return
   end subroutine

   !**** new **********************************************************************
   subroutine addrhs (brhs,elresf,lm,nee)
      !.... program to add element residual-force vector to
      !        global right-hand-side vector

      implicit none

      real*8  :: brhs(*),elresf(*)
      integer*4:: lm(*)
      integer*4:: nee

      integer*4:: k, j

      do j=1,nee
         k = lm(j)
         if (k.gt.0) then 
            brhs(k) = brhs(k) + elresf(j)
            ! write(*,*) k, brhs(k), j, elresf(j)
         end if
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

   ! **** new *********************************************************************
   subroutine escreverSistemaSkylineEmMTX(alhs,brhs,idiag,lm,id,ien,nee,nen,numel,numnp,neq,nomeArq)
      ! .... program to add element left-hand-side matrix to
      !         global left-hand-side matrix

      !         ldiag = .true.,  add diagonal element matrix
      !         ldiag = .false., then
      !           lsym = .true., add upper triangle of full element matrix
      !           lsym = .false., add full element matrix

      implicit none

      ! .... deactivate above card(s) for single-precision operation

      integer*4, intent(in) :: nee, nen, numel, numnp, neq
      real*8, intent(in)  :: alhs(:), brhs(:)
      integer*4:: idiag(*),lm(nee,numel)
      integer*4:: id(1,numnp),ien(nen,numel)
      character (len=*)  :: nomeArq
      logical :: ldiag,lsym

      integer*4:: i, j, l, k, m, p, eq
      integer*4::  nel, no, lAnterior, ned, node
      integer*4:: jjlast, jj, jcolht   
      integer*4:: luSist = 2000
      character(len=40), parameter :: formatoEscritaA='(2(i0,1x),e15.8)'
      character(len=40), parameter :: formatoEscritaB='(e15.8)'
      character(len=40), parameter :: formatoEscritaC='(e15.8, a)'

      ! formato mtx, para o matlab
      ! %%MatrixMarket matrix coordinate real symmetric
      ! % Generated 14-Aug-2014
      ! 12 12 41
      ! 1 1  1.333333333
      ! 3 1 -0.3333333333
      ! 8 1 -0.1666666667

      write(*,*) ' +++, em subroutine escreverSistemaSkylineEmMTX(alhs,idiag,lm,id,ien,nee,nen,numel,numnp)'
      open(file=nomeArq, unit=luSist) 

      write(luSist,'(a)')'%% matriz A de coeficientes reais simetrica positiva definida ' 
      write(luSist,'(a)')'% produzida pelo metodo classico de galerkin para o metodo de elementos finitos '
      write(luSist,'(a)')'% armazenamento em skyline '
      write(luSist,'(a,3(i0,a))') '%  matriz ',neq, 'X',neq, ' com ', idiag(neq), ' coefs armazenados' 
      write(luSist,'(3(i0,1x))')  neq, neq, idiag(neq)

      j = 1; i = 1
      write(luSist,formatoEscritaA)  j, j-(idiag(j)-i), alhs(i)
      do j=1,neq
         jj     = idiag(j)
         jjlast = jj
         jcolht = jj - jjlast

         if (jcolht.gt.1) then
            do i =  jjlast+1, jcolht-1 + jjlast+1
               write(luSist,formatoEscritaA) j,  j-(idiag(j)-i), alhs(i)
            end do
         end if
      enddo

      j = j - 1; i = i - 1
      ! write(*,*) j,  j -( idiag(j) - i), alhs(i)
      ! write(luSist,'( (a,i0,a))') '% lado direito com ',neq,' elementos' 
      ! write(luSist,'(  3(i10))' )  neq 
      i = 1
      write(luSist, formatoEscritaC ) brhs(i), "   BRHS"
      do i = 2, neq
         write(luSist, formatoEscritaB ) brhs(i)
      end do

      close(luSist) 

      return
   end subroutine

   !**** new **********************************************************************
   subroutine calcularResiduoSkyline(alhs,nalhs,vetX,vetBOriginal,neq,idiag,res)

      implicit none

      integer*4, intent(in)    :: nalhs, neq
      real*8,  intent(in)    :: alhs(nalhs), vetX(neq), vetBOriginal(neq)
      real*8,  intent(inout) :: res(neq)
      integer*4, intent(in)    :: idiag(neq)

      real*8  :: vetBcalculado(neq), norma
      integer*4:: i

      call matVecMulSkyLine(alhs,nalhs,vetX,neq,idiag,vetBcalculado)
      res=vetBOriginal-vetBcalculado
      norma = 0.0
      do i=1, neq
         norma= norma + res(i) * res(i)
      enddo
      norma = sqrt(norma) 
      write(*,'(a,1pe12.4,a)') "r=Au-b, norma do residuo da solucao, r=" , norma
 
   end subroutine

   !**** new **********************************************************************
   subroutine matVecMulSkyLine(alhsV,nalhsV,vetX,neqV,idiagV,vetBcalc)

      implicit none

      integer*4:: nalhsV, neqV
      real*8  :: alhsV(nalhsV), vetX(neqV), vetBcalc(neqV)
      integer*4:: idiagV(neqV)

      integer*4:: i,j, altura, jcolht, cont
      real*8  :: aij, xj
      integer*4:: jj, jjnext, istart,jtemp,jjlast

      vetBcalc=0.0
      altura=0
      cont=neqV
      jcolht=0
      jj = 0

      do j=1,neqV
         jjlast = jj
         jj     = idiagV(j)
         jcolht = jj - jjlast
         if (jcolht.gt.1) then
            vetBcalc(j) =  vetBcalc(j) + coldot(alhsV(jjlast+1),vetX(j-jcolht+1),jcolht-1)
         end if
      enddo

      do i=1,neqV
         aij=alhsV(idiagV(i))
         xj =vetX(i)
         vetBcalc(i)=vetBcalc(i)+aij*xj
      end do

      jjnext = idiagV(neqV)
      do j=neqV,2,-1
         jj     = jjnext
         jjnext = idiagV(j-1)
         jcolht = jj - jjnext
         if (jcolht.gt.1) then
            xj     = vetX(j)
            istart = j - jcolht + 1
            jtemp  = jjnext - istart + 1

            ! forma natural, por linhas: a cada novo i tem ui pronto
            ! ui =  ui + aij*xj, j =1, neq, i = 1, neq

            ! forma alternativa, por colunas: a cada novo i tem uma
            ! contribuicao de xi para o valor total ui 
            ! uj =  uj + aji*xi, j =1, neq, i = 1, neq,

            do i=istart,j-1
               vetBcalc(i) = vetBcalc(i) + alhsV(jtemp+i)*xj
            enddo
         endif
      end do
   end subroutine

end module

 module mSolverHypre

   integer*4 :: myid, num_procs, mpi_comm
      
      
   integer*8              :: A_HYPRE, parcsr_A, b_HYPRE, par_b, u_HYPRE, par_u, solverH
   integer*4              :: solver_id, precond_id
   integer*4              :: Clower, Cupper
   integer*4, allocatable :: rows(:)
   REAL*8,  allocatable   :: initialGuess(:)

#ifdef withHYPRE
   include 'mpif.h'
#endif

   integer*4 :: posPonteiro, contPonteiro, posColunas, posCoef
   integer*4 :: lda, nonzeros, nonzerosEst
   integer*8, parameter :: HYPRE_PARCSR=5555

   contains

   subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)

      integer*4, intent(inout) :: myid_, num_procs_
      integer*4 :: mpi_comm_
      integer*4 ::  ierr  

      write(*,'(a)') " ++++  em subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)"
#ifdef withHYPRE
      ! mpi_comm_ = 1000 ! MPI_COMM_WORLD
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid,      ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
      ! Convert the Fortran MPI communicator to a C version that can be used in hypre.
      ! Uncomment the line below if you are using MPICH
      myid_      = myid
      num_procs_ = num_procs
      mpi_comm_  = MPI_COMM_WORLD
      ! write(*,'(a,i12)') "i+++, mpi_comm_ =" ,mpi_comm_
      ! write(*,'(a,i12)') "i+++, MPI_COMM_WORLD =" ,MPI_COMM_WORLD
      ! Uncomment the line below if you are using LAM-MPI
      !call HYPRE_MPI_Comm_f2c(mpi_comm, MPI_COMM_WORLD, ierr)
#endif
      write(*,'(3(a,i0))')  " em subroutine inicializarMPI, myid= ", myid_,", &
               numprocs= ", num_procs_,", mpi_comm= ", mpi_comm_
   end subroutine inicializarMPI

   !=============================================================================
   subroutine finalizarMPI()
      integer*4 :: ierr=-1
#ifdef withHYPRE
      call MPI_Finalize(ierr)
#endif
   end subroutine finalizarMPI   

   !=============================================================================
   subroutine resolverSistemaAlgHYPRE  (A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_, &
               solver_id_, precond_id_, tol_, num_iterations_, final_res_norm_, &
               initialGuess_, brhs_, rows_, neq_, myid_, mpi_comm_)  

      implicit none
      integer*8, intent(in)         :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer*4, intent(in)         :: solver_id_, precond_id_
      double precision,  intent(in) :: tol_ ! = 1.0e-08 
      integer*4, intent(out)        :: num_iterations_
      double precision, intent(out) :: final_res_norm_
      integer*4, intent(in)         :: neq_
      integer*4, intent(in)         :: rows_(neq_)
      double precision, intent(in)  :: initialGuess_(neq_), brhs_(neq_)
      integer*4, intent(in)         :: myid_
      integer*4, intent(in)         :: mpi_comm_ 

      integer*4        :: precond
      integer*4        :: precond_id;
      integer*4        :: ierr
      integer*4, save  :: contador = 1
      integer*4        :: local_size , numMaxIterations = 10000
      integer*4        :: printLevel, solver_id
      integer*4        :: i, Flower, Fupper
      real*8           :: t1, t2, t3, t4, elapsedT 

      ! block
      character(len=15) :: nomeU 
      character(len=12) :: nomeB
      character(len=12) :: nomeA
          
      nomeU = "solucao00.out."
      nomeA = "Alhs00.out."
      nomeB = "brhs00.out."
      ! call HYPRE_IJVectorPrint( b_, nomeB, ierr)
      ! call HYPRE_IJMatrixPrint( A_, nomeA, ierr)
      ! end block
      write(*,*) "em subroutine resolverSistemaAlgHYPRE  (A_, parcsr_A_, b_, par_b_, "
      ! write(*,'(a,7i15)')  "s+++, ", A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      ! write(*,'(3(a,i0))') "s+++, em resolverSistemaAlgHYPRE, myid= ", myid,", &
               ! numprocs= ", num_procs,", mpi_comm= ", mpi_comm
      local_size = neq_
      printLevel = 0
#ifdef withHYPRE
      call HYPRE_IJVectorSetValues (u_, neq_, rows_, initialGuess_, ierr)
      ! call HYPRE_IJVectorSetValues (b_, neq_, rows_, brhs_, ierr)
      ! call HYPRE_IJVectorPrint( b_, nomeB, ierr)
      ! call HYPRE_IJVectorPrint( A_, nomeA, ierr)

      ! write(*,*) rows_(1:neq_)
      ! write(*,*) brhs_(1:neq_)
      ! write(*,*) "solver_id_=", solver_id_
      ! write(*,*) "precond =",  precond
      ! write(*,*) "precond_id_ =",  precond_id_
      ! write(*,*) "neq_ =", neq_
      ! write(*,*) "++ precond_id_ =", precond_id_,", tol =", tol_

      call timing(t1)
      if     (solver_id_==0) then
         write(*,*) " ... BoomerAMGCreate , BD"
         call solverOptionA ()

      elseif (solver_id_==1) then
         write(*,*) " ... PCG with AMG preconditioner, BD_dez";
         call solverOptionB ()

      elseif (solver_id_==2) then
         write(*,*) " ... PCG with ParaSails"
         call solverOptionC ()
      else
         if (myid_ .eq. 0) then 
            print *,'Invalid solver id specified'
            stop
         endif  
      endif
      call timing(t2)
      elapsedT = t2 - t1

      contains

      subroutine solverOptionA
         ! Create solver
         call HYPRE_BoomerAMGCreate          (solver_, ierr)

         ! Set some parameters (See Reference Manual for more parameters)
         ! print solve info + parameters 
         call HYPRE_BoomerAMGSetPrintLevel   (solver_, printLevel,      ierr)  
         ! Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (solver_, 6,      ierr) 
         ! G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (solver_, 6,      ierr)     
         ! Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (solver_, 1,      ierr)  
         ! maximum number of levels 
         call HYPRE_BoomerAMGSetMaxLevels    (solver_, 200,     ierr) 
         ! conv. tolerance
         call HYPRE_BoomerAMGSetTol          (solver_, tol_, ierr)    

         ! Now setup and solve!
         call HYPRE_BoomerAMGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr )
         call HYPRE_BoomerAMGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr )

         ! Run info - needed logging turned on 
         call HYPRE_BoomerAMGGetNumIterations(solver_, num_iterations_, ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(solver_, final_res_norm_, ierr)

         ! Destroy solver_
         call HYPRE_BoomerAMGDestroy         (solver_, ierr )

         ! PCG with AMG preconditioner
      end subroutine solverOptionA

      subroutine solverOptionB ()
         ! Create solver_
         call HYPRE_ParCSRPCGCreate        (MPI_COMM_WORLD, solver_, ierr)
         ! Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter    (solver_, numMaxIterations,   ierr)
         call HYPRE_ParCSRPCGSetTol        (solver_, tol_, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (solver_, 1,      ierr)
         ! call HYPRE_ParCSRPCGSetPrintLevel (solver_, 2,      ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (solver_, printLevel,      ierr)

         ! Now set up the AMG preconditioner and specify any parameters

         call HYPRE_BoomerAMGCreate          (precond, ierr)

         ! Set some parameters (See Reference Manual for more parameters)

         ! Relaxation Parameters:
         ! Visiting Grid:                     down   up  coarse
         ! Number of sweeps:                     1    1     1 
         ! Type 0=Jac, 3=hGS, 6=hSGS, 9=GE:      6    6     6 

         ! print less solver_ info since a preconditioner
         ! call HYPRE_BoomerAMGSetPrintLevel   (precond, printLevel,     ierr); 
         ! Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (precond, 6,     ierr) 
         ! SYMMETRIC G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (precond, 6,     ierr)     
         ! Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (precond, 1,     ierr)  
         ! conv. tolerance
         call HYPRE_BoomerAMGSetTol          (precond, 0.0d-6, ierr)     
         ! do only one iteration! 
         call HYPRE_BoomerAMGSetMaxIter      (precond, 1,     ierr)

         ! set amg as the pcg preconditioner
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)

         ! Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         ! Run info - needed logging turned on 
         call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
         call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)
         write(*,*) "num_iterations= ", num_iterations_; 

         ! Destroy precond and solver
         call HYPRE_BoomerAMGDestroy          (precond, ierr )
         call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
      end subroutine solverOptionB

      subroutine solverOptionC()
         ! Create solver
         call HYPRE_ParCSRPCGCreate          (MPI_COMM_WORLD, solver_, ierr)

         ! Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter      (solver_, 1000, ierr)
         call HYPRE_ParCSRPCGSetTol          (solver_, 1.0d-7, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm      (solver_, 1, ierr)
         call HYPRE_ParCSRPCGSetPrintLevel   (solver_, printLevel, ierr)
         ! call HYPRE_ParCSRPCGSetLogging      (solver_, 1, ierr)

         ! Now set up the Parasails preconditioner and specify any parameters
         call HYPRE_ParaSailsCreate          (MPI_COMM_WORLD, precond,ierr)
         call HYPRE_ParaSailsSetParams       (precond, 0.1d0, 1, ierr)
         call HYPRE_ParaSailsSetFilter       (precond, 0.05d0, ierr)
         call HYPRE_ParaSailsSetSym          (precond, 1)
         ! call HYPRE_ParaSailsSetLogging      (precond, 3, ierr)

         ! set parsails as the pcg preconditioner
         ! precond_id_ = 4
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)

         ! Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)

         ! Run info - needed logging turned on 

         call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
         call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)

         ! Destroy precond and solver
         call HYPRE_ParaSailsDestroy          (precond, ierr )
         call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
      end subroutine solverOptionC

#endif
   end subroutine resolverSistemaAlgHYPRE

   !=============================================================================
   subroutine escreverResultados_HYPRE (u_, num_iterations_, final_res_norm_, &
               elapsedT_, myId_, print_solution_)

      implicit none
      integer*8   :: u_
      integer*4        num_iterations_
      integer*4, intent(in) :: myId_
      integer*4           :: print_solution_ 

      double precision final_res_norm_, elapsedT_

      integer*4 ::  ierr  
      character (LEN=16)  :: SolverStatus="incomplete"

      ! write(*,*) 'em subroutine escreverResultados_HYPRE'
      SolverStatus="complete"
      if(SolverStatus=="complete" .and. myid_ .eq. 0) then
         print *
         print *, "Final Relative Residual Norm = ", final_res_norm_
         print *, "Iterations                   = ", num_iterations_
         print *, 'Elapsed real time            = ', elapsedT_
         print *
      endif

      ! Print the solution
      ! print_solution_ = 0 
      ! print_solution_ = 1 
      if ( print_solution_ .ne. 0 ) then
#ifdef withHYPRE
         call HYPRE_IJVectorPrint( u_, "ij.out.u", ierr)
#endif
      endif
   end subroutine escreverResultados_HYPRE           

   !*************************************************************************************************
   subroutine addnslHYPREB(matA, A_HYPRE_,eleffm_, lm_, idiag_, nee_, diag_, lsym_)
      ! call addnslHYPRE(estrutSistEqP_%A_HYPRE, eleffm, lmLocal, estrutSistEqP_%idiag, nee, diag, lsym)

      ! program to add element left-hand-side matrix to          
      ! global left-hand-side matrix                      
      ! ldiag = .true.,  add diagonal element matrix              
      ! ldiag = .false, then                                     
      ! add full nonsymmetric element matrix                   

      implicit none
      
      integer*8:: A_HYPRE_ 
      real*8  ::  eleffm_(:,:)
      integer*4:: lm_(:) , idiag_(:), nee_
      logical :: diag_, lsym_

      integer*4, parameter:: numCoefPorLinha=100 
      integer*4 :: i,j,k,l, eq, c, nnz, ierr
      integer*4 ::  cols(numCoefPorLinha)
      real*8  :: values(numCoefPorLinha) 
      integer*4, save :: contador = 0
      real :: matA (20,20)

      contador = contador + 1
      do j=1,nee_
         nnz = 0
         eq = abs(lm_(j))
         if(eq > 0) then 
            do i=1,nee_
               c = abs(lm_(i))
               if(c > 0) then 
                  nnz = nnz + 1
                  cols(nnz) = c
                  values(nnz) = eleffm_(j,i)
                  matA(eq,c) = matA(eq,c) + eleffm_(j,i)
               end if
            enddo
            cols = cols - 1
            eq = eq - 1
#ifdef withHYPRE
            ! write(*,'(i0, a,i3,a,10i3)') A_HYPRE_, ",  eq = ", eq,  ", ", cols(1:nnz)
            ! write(*,'(i0, a,i3,a,10e12.3)')A_HYPRE_, ",  eq = ", eq,  ", ", values(1:nnz)
            ! write(*,*) "em ne addnslHYPRE, estrutSistEq_%A_HYPRE=", estrutSistEq_%A_HYPRE, nnz, eq, cols(1:nnz), values(1:nnz) 
            call HYPRE_IJMatrixAddToValues(A_HYPRE_, 1, nnz, eq, cols, values, ierr)
#endif
         end if
      enddo

      return
   end subroutine addnslHYPREB
      

   subroutine addnslHYPRE(A_, eleffm,lmT,nee,nel)
      ! program to add element left-hand-side matrix to          
      ! global left-hand-side matrix                      
                                                                 
      ! ldiag = .true.,  add diagonal element matrix              
      ! ldiag = .false, then                                     
      ! add full nonsymmetric element matrix                   
                                                                 
      implicit none

      integer*8, intent(inout) :: A_
      integer*4:: nee, nel
      real*8  :: eleffm(nee,*)
      integer, pointer :: lmT(:,:,:)
      integer*4:: lm(nee)

      integer*4, parameter:: numCoefPorLinha=100 
      integer*4:: i,j,k,l, eq, c, nnz, ierr 
      integer*4::  cols(numCoefPorLinha)
      real*8  :: values(numCoefPorLinha) 
      integer, save :: contador = 0

#ifdef withHYPRE
      lm(:)=reshape(lmT(:,:,nel),(/nee/));
   
      contador = contador + 1
      ! write(*,*) contador, "lm = ", lm(:)
      do j=1,nee
         nnz = 0
         eq = abs(lm(j))
         if(eq > 0) then 
            do i=1,nee
               c = abs(lm(i))
               if(c > 0) then 
                  nnz = nnz + 1
                  cols(nnz) = c
                  values(nnz) = eleffm(j,i)
               end if
            enddo
            cols = cols - 1
            eq = eq - 1
            ! write(*,'(a,i5,a,10i5)') " eq = ", eq,  ", ", cols(1:nnz)
            ! call HYPRE_IJMatrixAddToValues(A_, 1, nnz, eq, cols, values, ierr)
         end if
      enddo
#endif

      return
   end subroutine addnslHYPRE

   ! *************************************************************************************************
   ! subroutine addnslHYPRE( A_HYPRE_,eleffm_, lm_, idiag_, nee_, diag_, lsym_)
      ! call addnslHYPRE(estrutSistEqP_%A_HYPRE, eleffm, lmLocal, estrutSistEqP_%idiag, nee, diag, lsym)
   
      ! program to add element left-hand-side matrix to          
      ! global left-hand-side matrix                      
      ! ldiag = .true.,  add diagonal element matrix              
      ! ldiag = .false, then                                     
      ! add full nonsymmetric element matrix                   

      ! implicit none
      
      ! integer*8:: A_HYPRE_ 
      ! real*8  ::  eleffm_(:,:)
      ! integer*4:: lm_(:) , idiag_(:), nee_
      ! logical :: diag_, lsym_

      ! integer*4, parameter:: numCoefPorLinha=100 
      ! integer*4 :: i,j,k,l, eq, c, nnz, ierr
      ! integer*4 ::  cols(numCoefPorLinha)
      ! real*8  :: values(numCoefPorLinha) 
      ! integer*4, save :: contador = 0

      ! contador = contador + 1
      ! do j=1,nee_
      !    nnz = 0
      !    eq = abs(lm_(j))
      !    if(eq > 0) then 
      !       do i=1,nee_
      !          c = abs(lm_(i))
      !          if(c > 0) then 
      !             nnz = nnz + 1
      !             cols(nnz) = c
      !             values(nnz) = eleffm_(j,i)
      !          end if
      !       enddo
      !       cols = cols - 1
      !       eq = eq - 1
! #ifdef withHYPRE
!             ! write(*,'(i0, a,i3,a,10i3)') A_HYPRE_, ",  eq = ", eq,  ", ", cols(1:nnz)
!             ! write(*,'(i0, a,i3,a,10e12.3)')A_HYPRE_, ",  eq = ", eq,  ", ", values(1:nnz)
!             ! write(*,*) "em ne addnslHYPRE, estrutSistEq_%A_HYPRE=", estrutSistEq_%A_HYPRE, nnz, eq, cols(1:nnz), values(1:nnz) 
!             call HYPRE_IJMatrixAddToValues(A_HYPRE_, 1, nnz, eq, cols, values, ierr)   
! #endif
   !       end if
   !    enddo
      
   !    return
   ! end subroutine addnslHYPRE

   ! *************************************************************************************************
   ! subroutine addnslHYPRE01(estrutSistEq_, eleffm, lm, nee)
   !    ! program to add element left-hand-side matrix to          
   !    ! global left-hand-side matrix                      
                                                                 
   !    ! ldiag = .true.,  add diagonal element matrix              
   !    ! ldiag = .false, then                                     
   !    ! add full nonsymmetric element matrix                                                                 

   !    ! implicit none

   !    ! type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_  
      
   !    ! real*8  :: eleffm(:,:)
   !    ! integer*4:: lm(:), nee

   !    ! integer*4, parameter:: numCoefPorLinha=100 
   !    ! integer*4 :: i,j,k,l, eq, c, nnz, ierr
   !    ! integer*4 ::  cols(numCoefPorLinha)
   !    ! real*8  :: values(numCoefPorLinha) 
   !    ! integer*4, save :: contador = 0

   !    ! integer*4, allocatable:: lmLocal(:)
   !    ! nee = size(estrutSistEq_%lm,1)*size(estrutSistEq_%lm,2)
   !    ! allocate (lmLocal(nee))
   !    ! lmLocal(:)=reshape(estrutSistEq_%lm(:,:,nel),(/nee/))
   
   !    ! contador = contador + 1
   !    ! write(*,*) nel, "lm = ", lm(:)
   !    do j=1,nee
   !       nnz = 0
   !       eq = abs(lm(j))
   !       if(eq > 0) then 
   !          do i=1,nee
   !             c = abs(lm(i))
   !             if(c > 0) then 
   !                nnz = nnz + 1
   !                cols(nnz) = c
   !                values(nnz) = eleffm(j,i)
   !             end if
   !          enddo
   !          cols = cols - 1
   !          eq = eq - 1
! #ifdef withHYPRE
!             ! write(*,'(i0, a,i3,a,10i3)') estrutSistEq_%A_HYPRE, ",  eq = ", eq,  ", ", cols(1:nnz)
!             ! write(*,'(i0, a,i3,a,10e12.3)') estrutSistEq_%A_HYPRE, ",  eq = ", eq,  ", ", values(1:nnz)
!             ! write(*,*) "em ne addnslHYPRE, estrutSistEq_%A_HYPRE=", estrutSistEq_%A_HYPRE, nnz, eq, cols(1:nnz), values(1:nnz) 
!             call HYPRE_IJMatrixAddToValues(estrutSistEq_%A_HYPRE, 1, nnz, eq, cols, values, ierr)
! #endif
   !       end if
   !    enddo

   !    return
   ! end subroutine addnslHYPRE01

   !**** new **********************************************************************
   subroutine addrhsHYPRE (b_, brhs_, elresf, lm, nee_)
      !.... program to add element residual-force vector to
      ! global right-hand-side vector

      implicit none

      integer*8 , intent(inout) ::  b_
      real*8, intent(in)  :: brhs_(*)
      integer*4:: nee_
      real*8  :: elresf(nee_)
      integer*4:: lm(:)

      integer*4:: k, j, nnz, ierr
      integer*4:: rows(nee_+20), i
      real*8   :: values(nee_+20)
      ! integer, save :: nel = 0

      write(*,*) " em subroutine addrhsHYPRE (b_, brhs_, elresf,lm,nee_)"
      values = 0.0
      ! nel = nel + 1
      ! write(*,'(i5,a,16i4)') nel, "lm = ", lm(:)
      nnz = 0     
      do j=1,nee_
         k = lm(j)
         if (k.gt.0) then
            nnz = nnz + 1
            values(nnz) = elresf(j) + brhs_(k)  
            rows  (nnz) = k-1
         end if
      enddo

      ! rows = rows - 1
      ! write(*,'(a,i3,a,16i6)') " nel = ", nel,  ", ", rows(1:nnz)
      ! write(*,'(a,i3,a,i9,16e15.4)') " nel = ", nel,  ", ", b_, values(1:nnz)
      ! call HYPRE_IJVectorSetValues (b_, nnz, rows, elresf, ierr)
#ifdef withHYPRE
      ! write(*,*) "rows =", rows(1:nnz)
      ! write(*,*) "values =", values(1:nnz)
      call HYPRE_IJVectorAddToValues (b_, nnz, rows, values, ierr)
#endif
      ! write(*,'(a,4i5)') "+++++  nee_ = ", nee_
      ! write(*,'(a,4i5)') "+++++",   lm(1:nee_)
      ! write(*,'(a,4i5)') "+++++",   rows(1:nnz)
      ! write(*,'(a,4f8.4)') "+++++", elresf(1:nee_)
      ! write(*,'(a,4f8.4)') "+++++", values(1:nnz)
      return
   end subroutine addrhsHYPRE

   subroutine lerValoresSistemaAlgHYPRE (alhs_ , Ap_, Ai_, rhs_, x_, rows_, neq_, &
            nonzerosT_, Clower_, Cupper_, A_, parcsr_A_, b_, par_b_, u_, par_u_, & 
            solver_, myid_, mpi_comm_) 

      implicit none
      integer*4, intent(in)  :: nonzerosT_
      integer*4, intent(in)  :: neq_
      integer*4, intent(in)  :: Clower_, Cupper_
      double precision       :: alhs_(nonZerosT_)  
      integer*4              :: Ap_(neq_), Ai_(nonZerosT_)  
      double precision       :: rhs_(neq_), x_(neq_)
      integer*4              :: rows_(neq_)  
      integer*8              :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer*4, intent(in)  :: myid_
      integer*4, intent(in)  :: mpi_comm_

      integer*4     :: nnz, i, j, k, eq, ierr, printLevel
      integer*4     :: local_size, Flower, Fupper
      integer*4     :: cols(90)
      integer*4     :: colsB(90), eqB
      double precision ::  values(90)
     
#ifdef withHYPRE
      ! write(*,*) " em subroutine lerValoresSistemaAlgCSR_HYPRE ( A_, parcsr_A_, b_, par_b_, u_, par_u_, .."
      print*, " atribuindo valores da matrix para o HYPRE_IJMatrix "
      Flower    = Clower_+1; Fupper    = Cupper_+1;
      ! HYPRE_IJMatrixRead( <filename>, MPI_COMM_WORLD, HYPRE_PARCSR, &A );
      local_size = Cupper_ - Clower_ + 1

      ! HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD, HYPRE_PARCSR, &b ); 
      print*, " ++++ B atribuindo valores de RHS para o HYPRE_IJVector ", Flower, Fupper
      ! call HYPRE_IJVectorSetValues (b_, local_size, rows_, rhs_, ierr )
      ! call HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD, par_b_, b_ ); 

      print*, " +++++ atribuindo valores de guess para o HYPRE_IJVector "
      ! call HYPRE_IJVectorSetValues (u_, local_size, rows_, x_, ierr)
      ! call HYPRE_IJVectorRead(  <filename>,  , MPI_COMM_WORLD, par_u_, u_ ); 
#endif
   end subroutine lerValoresSistemaAlgHYPRE

   !=============================================================================
   ! subroutine resolverSistemaAlgHYPRE01  (estrutSistEq_ , tol_, num_iterations_, final_res_norm_) 
   !    implicit none
   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
   !    double precision,  intent(in) :: tol_
   !    double precision,  intent(out) :: final_res_norm_ 
   !    integer*4, intent(out)        :: num_iterations_
   !    write(*,*) "em subroutine resolverSistemaAlgHYPRE02  (estrutSistEq_ , tol_, ... " 
   !    call resolverSistemaAlgHYPRE( estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A, & 
   !             estrutSistEq_%b_HYPRE, estrutSistEq_%par_b, estrutSistEq_%u_HYPRE, &
   !             estrutSistEq_%par_u, estrutSistEq_%solver, estrutSistEq_%solver_id, &
   !             estrutSistEq_%precond_id, tol_, num_iterations_, final_res_norm_, &
   !             estrutSistEq_%initialGuess, estrutSistEq_%brhs, estrutSistEq_%rows, &
   !             estrutSistEq_%neq, myid, mpi_comm)
   ! end subroutine resolverSistemaAlgHYPRE01

   !=============================================================================
   ! subroutine criarVetorBRHS_HYPRE  (estrutSistEq_, mpi_comm_) 
   !    implicit none
   !    integer*4, intent(in)  :: mpi_comm_
   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
   !    call criarVetor_HYPRE(estrutSistEq_%b_HYPRE, estrutSistEq_%Clower, estrutSistEq_%Cupper, mpi_comm_)
   ! end subroutine criarVetorBRHS_HYPRE  

   !=============================================================================
   ! subroutine criarVetorSolucao_HYPRE  (estrutSistEq_, mpi_comm_) 
   !    implicit none
   !    integer*4, intent(in)  :: mpi_comm_
   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
   !    call criarVetor_HYPRE(estrutSistEq_%u_HYPRE, estrutSistEq_%Clower, estrutSistEq_%Cupper, mpi_comm_)
   ! end subroutine criarVetorSolucao_HYPRE  

   !=============================================================================
   subroutine criarVetor_HYPRE(v_, Clower_, Cupper_, mpi_comm_)
      integer*8, intent(out) :: v_
      integer*4, intent(in)  :: Clower_, Cupper_
      integer*4, intent(in)  :: mpi_comm_
      integer*4 :: ierr 
      ! write(*,'(a,i10)') "cV+++, ", v_
#ifdef withHYPRE
      call HYPRE_IJVectorCreate        (mpi_comm_, Clower_, Cupper_, v_, ierr )
      call HYPRE_IJVectorSetObjectType (v_, HYPRE_PARCSR, ierr)
      call HYPRE_IJVectorInitialize    (v_, ierr)
#endif
      ! write(*,'(a,i10)') "cV+++, ", v_
   end subroutine criarVetor_HYPRE

   !=============================================================================
   subroutine criarMatriz_HYPRE(M_, Clower_, Cupper_,  mpi_comm_)
      integer * 8              :: M_
      integer * 4, intent(in)  :: Clower_, Cupper_
      integer * 4, intent(in)  :: mpi_comm_
      integer, parameter ::  HYPRE_PARCSR=5555
      integer :: ierr 
      ! write(*,*) " cM+++ , Clower_= ", Clower_, ", Cupper_ " , Cupper_
#ifdef withHYPRE
      A_=-9
      ! Create the matrix.
      ! Note that this is a square matrix, so we indicate the row partition
      ! size twice (since number of rows = number of cols)
      call HYPRE_IJMatrixCreate        (mpi_comm_, Clower_, Cupper_, Clower_, Cupper_, M_, ierr )
      ! Choose a parallel csr format storage (see the User's Manual)
      call HYPRE_IJMatrixSetObjectType (M_, HYPRE_PARCSR, ierr)
      ! call HYPRE_StructMatrixSetSymmetric (M_, 1)
      ! Initialize before setting coefficients
      call HYPRE_IJMatrixInitialize    (M_, ierr)
#endif
      ! write(*,'(7(a,i0))') "cM+++, ", M_ , ", Clower_= ", Clower_, ",  Cupper_=" , Cupper_
      if(A_==0) then
         write(*,*) " em subroutine criarMatriz_HYPRE, erro, M_ = ", M_
         stop
      end if
   end subroutine criarMatriz_HYPRE

   !=============================================================================
   ! subroutine criarMatriz_HYPRE01  (estrutSistEq_, mpi_comm_) 
   !    implicit none
   !    integer*4, intent(in)  :: mpi_comm_
   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
   !    call criarMatriz_HYPRE (estrutSistEq_%A_HYPRE, estrutSistEq_%Clower, estrutSistEq_%Cupper, mpi_comm_)
   ! end subroutine criarMatriz_HYPRE01  

   !=============================================================================
   ! subroutine fecharMatriz_HYPRE01  (estrutSistEq_) 
   !    implicit none
   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
   !    write(*,*) "em subroutine fecharMatriz_HYPRE01  (estrutSistEq_ , tol_, ... " 
   !    call fecharMatriz_HYPRE( estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A)
   ! end subroutine fecharMatriz_HYPRE01

   !=============================================================================
   subroutine fecharMatriz_HYPRE(A_, parcsr_A_) !
      implicit none
      integer*8, intent(in)  :: A_
      integer*8, intent(out) :: parcsr_A_
      integer*4 :: ierr
#ifdef withHYPRE
      ! Assemble after setting the coefficients
      call HYPRE_IJMatrixAssemble( A_, ierr)
      ! Get parcsr matrix object
      call HYPRE_IJMatrixGetObject( A_, parcsr_A_, ierr)
#endif
      ! write(*,'(2(a,i0))') "fM+++, A_ ", A_, ", parcsr_A_ =", parcsr_A_
   end subroutine fecharMatriz_HYPRE

   !=============================================================================
   ! subroutine fecharVetorBRHS_HYPRE  (estrutSistEq_) 
   !    implicit none
   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
   !    call fecharVetor_HYPRE(estrutSistEq_%b_HYPRE, estrutSistEq_%par_b)
   ! end subroutine fecharVetorBRHS_HYPRE  

   !=============================================================================
   ! subroutine fecharVetorSolucao_HYPRE  (estrutSistEq_) 
   !    implicit none
   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
   !    call fecharVetor_HYPRE(estrutSistEq_%u_HYPRE, estrutSistEq_%par_u)
   ! end subroutine fecharVetorSolucao_HYPRE  

   !=============================================================================
   subroutine fecharVetor_HYPRE (v_, par_v_)
      implicit none
      integer*8, intent(in)  :: v_
      integer*8, intent(out) :: par_v_
      integer*4 :: ierr
#ifdef withHYPRE
      ! Assemble after setting the coefficients
      call HYPRE_IJVectorAssemble  (v_, ierr)

      ! Get parcsr vector object
      call HYPRE_IJVectorGetObject (v_, par_v_, ierr)
#endif
      ! write(*,'(2(a,i0))') "fV+++, v_ ", v_, ", par_v =", par_v_
   end subroutine fecharVetor_HYPRE

   !=============================================================================
   subroutine destruirVetor_HYPRE   ( v_)
      integer*8 , intent(in) ::  v_
      integer*4 ::  ierr  
      ! write(*,'(a,i0)') "dV+++, v=",  v_ 
#ifdef withHYPRE
      call HYPRE_IJVectorDestroy(v_, ierr)
#endif
   end subroutine destruirVetor_HYPRE       

   !=============================================================================
   subroutine destruirMatriz_HYPRE   (M_)
      integer*8 , intent(inout) :: M_
      integer*4 ::  ierr  
      ! write(*,*) " ++, em destruirMatriz_HYPRE   (M_)"
      ! write(*,'(a,i0)') "dM+++, ", M__ 
#ifdef withHYPRE
      call HYPRE_IJMatrixDestroy(M_, ierr)
#endif
      M_=0;
   end subroutine destruirMatriz_HYPRE       

   !=============================================================================
   subroutine extrairValoresVetor_HYPRE(v_HYPRE_, Flower_, Fupper_, rows_, destino_) 
      integer*8, intent(out) :: v_HYPRE_
      integer*4, intent(in)  :: Flower_, Fupper_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(out)    :: destino_(:)
      integer*4 :: ierr 
      write(*,'(3(a,i0))') "extV+++, ", v_HYPRE_, ", Flower_= ", Flower_, ", Fupper_= ", Fupper_
#ifdef withHYPRE
      call HYPRE_IJVectorGetValues (v_HYPRE_, Fupper_- Flower_ + 1, rows_, destino_(Flower_), ierr)
#endif
   end subroutine extrairValoresVetor_HYPRE

   subroutine atribuirValoresVetor_HYPRE(v_HYPRE_, Flower_, Fupper_, rows_, destino_) 
      integer*8, intent(out) :: v_HYPRE_
      integer*4, intent(in)  :: Flower_, Fupper_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(in)    :: destino_(:)
      integer*4 :: ierr 
      write(*,'(3(a,i0))') "atrV+++, ", v_HYPRE_, ", Flower_= ", Flower_, ", Fupper_= ", Fupper_
#ifdef withHYPRE
      call HYPRE_IJVectorSetValues (v_HYPRE_, Fupper_- Flower_ + 1, rows_, destino_(Flower_), ierr)
#endif
   end subroutine atribuirValoresVetor_HYPRE

   subroutine adicionarValoresVetor_HYPRE(v_HYPRE_, Flower_, Fupper_, rows_, destino_) 
      integer*8, intent(out) :: v_HYPRE_
      integer*4, intent(in)  :: Flower_, Fupper_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(in)    :: destino_(:)
      integer*4 :: ierr 
      write(*,'(3(a,i0))') "addV+++, ", v_HYPRE_, ", Flower_= ", Flower_, ", Fupper_= ", Fupper_
#ifdef withHYPRE
      call HYPRE_IJVectorAddToValues (v_HYPRE_, Fupper_- Flower_ + 1, rows_, destino_(Flower_), ierr)
#endif
   end subroutine adicionarValoresVetor_HYPRE

   subroutine escreverResultadosHYPRE (x_, num_iterations_, final_res_norm_, elapsedT_, &
            myId_, print_solution_)

      implicit none
      integer*8   :: x_
      integer*4        num_iterations_
      integer*4, intent(in) :: myId_
      integer*4           :: print_solution_ 

      double precision final_res_norm_, elapsedT_

      integer*8 ::  ierr  
      character (LEN=16)  :: SolverStatus="incomplete"

      write(*,*) 'em subroutine escreverResultadosHYPRE'
      SolverStatus="complete"
      if(SolverStatus=="complete" .and. myid_ .eq. 0) then
         print *
         print *, "Final Relative Residual Norm = ", final_res_norm_
         print *, "Iterations                   = ", num_iterations_
         print *, 'Elapsed real time            = ', elapsedT_
         print *
      endif

      ! Print the solution
      ! print_solution_ = 0 
      if ( print_solution_ .ne. 0 ) then
#ifdef withHYPRE
         call HYPRE_IJVectorPrint( x_, "ij.out.x", ierr)
#endif
      endif

   end subroutine escreverResultadosHYPRE           

end module mSolverHypre
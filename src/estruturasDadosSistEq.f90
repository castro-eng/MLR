!> A module for storing and controlling the method's data structures.
module mEstruturasDadosSistEq
   !> Structure for storing the model data

   !> @param optSolver - A flag for wich type of solver are to be used!
   !>                    Skyline (LU decompsition), HYPRE, Pardiso
   !> @param        pt - A vector of integers
   !> @param     iparm - A vector of integers
   !> @param listaElemPorNo - Variable that correlates the node number
   !>                         with its respective elements. 
   !> @param     dparm - A vector of doubles

   !> @param  uInicial - Inicial reference for the velocity
   !> @param       neq - Number of equations in the system
   !> @param     nalhs - A vector of doubles
   !> @param      ndof - Number of degrre of freddom in the system
   !> @param    nlvect - Number of neumann aplied conditions
   !> @param        lm - 
   !> @param lmFratura - 
   !> @param   lmBorda - 
   !> @param lmFraturaNaBorda - 
   !> @param        id - 
   !> @param     idiag - A vector for the matrix manipulation
   !> @param         u - A vector to storage the atual velocity
   !> @param uTempoAnt - A vector to storage the last time step for the velocity
   !> @param  uIterAnt - A vector to storage the last iteration for the velocity
   !> @param      uRef - A vector of reference velocity
   !> @param         f - A vector for the aplied forces
   !> @param      brhs - 
   !> @param      alhs - 
   !> @param        Ap - 
   !> @param        Ai - 
   !> @param LMstencilEq - 
   !> @param numCoefPorLinha - 
   !> @param nVizinMax - 

   type estruturasArmazenamentoSistemaEq
   character(len=7)  :: optSolver 

   integer*4          :: pt(64), iparm(64)
   integer*4, pointer :: listaElemPorNo(:,:)

   real*8 :: dparm(64)
   real*8 :: uInicial

   integer*4          :: neq, nalhs                         ! skyline + CSR + HYPRE
   integer*4          :: ndof, nlvect                       ! skyline + CSR + HYPRE
   integer*4, pointer :: lm(:,:,:)=>null()                  ! skyline + CSR + HYPRE
   integer*4, pointer :: lmFratura(:,:,:)=>null()           ! skyline + CSR + HYPRE
   integer*4, pointer :: lmBorda(:,:,:)=>null()             ! skyline + CSR + HYPRE
   integer*4, pointer :: lmFraturaNaBorda(:,:,:)=>null()    ! skyline + CSR + HYPRE
   integer*4, pointer :: id(:,:) => null()                  ! skyline + CSR + HYPRE   
   integer*4, pointer :: idiag(:)=>null()                   ! Skyline + CSR

   real*8, pointer :: u(:,:),  uTempoAnt(:,:), uIterAnt(:,:)   ! skyline + CSR + HYPRE
   real*8, pointer :: uRef(:,:), f(:,:,:)                      ! skyline + CSR + HYPRE
   real*8, pointer :: brhs(:)=>null()                          ! skyline + CSR + HYPRE
   real*8, pointer :: alhs(:)=>null()                          ! skyline + CSR + HYPRE

   integer*4, pointer :: Ap(:)=>null()               ! CSR
   integer*4, pointer :: Ai(:)=>null()               ! CSR
   integer*4, pointer :: LMstencilEq(:,:) => null()  ! CSR
   integer*4          :: numCoefPorLinha, nVizinMax  ! CSR

   ! integer*4              :: solver_id, precond_id
   ! integer*4              :: Clower, Cupper
   ! integer*4, allocatable :: rows(:)
   ! integer*8              :: A_HYPRE, parcsr_A, b_HYPRE, par_b, u_HYPRE, par_u, solverH

   ! real*8,  allocatable   :: initialGuess(:)                    
   end type

   contains

   !=============================================================================
      !> Routine to form lm array
      !> @param    id - 
      !> @param conecElem -
      !> @param    lm -
      !> @param  ndof -
      !> @param   ned -
      !> @param   nen -
      !> @param numel -
      subroutine formlm (id,conecElem,lm,ndof,ned,nen,numel)
         
         implicit none
      
         integer*4:: ndof,ned,nen,numel
         integer*4:: id(ndof,*),conecElem(nen,*),lm(ned,nen,*)
         integer*4:: i,j,k,node
      
         do k=1,numel
            do j=1,nen
               node=conecElem(j,k)
               do i=1,ndof
                  lm(i,j,k) = id(i,node)
               enddo
            enddo
         enddo
      
         return
      end subroutine

      !**** new **********************************************************************
      !> Routine to compute diagonal addresses of left-hand-side matrix
      !> @param idiag - 
      !> @param   neq -
      !> @param     n -
      subroutine diag(idiag,neq,n)

         implicit none

         integer*4:: neq
         integer*4:: n
         integer*4:: idiag(neq)

         integer*4:: i

         n = 1
         idiag(1) = 1 
         if (neq.eq.1) return

         do i=2,neq
            idiag(i) = idiag(i) + idiag(i-1) + 1
         enddo
         n = idiag(neq)

         return
      end subroutine

      !**** new **********************************************************************
      !> Routine to compute column heights in global left-hand-side matrix
      !> @param idiag -
      !> @param    lm -
      !> @param   ned -
      !> @param   nen - 
      !> @param numel -
      !> @param   neq -
      subroutine colht(idiag,lm,ned,nen,numel,neq)

         implicit none

         integer*4:: ned, nen, numel
         integer*4:: neq
         integer*4:: idiag(*),lm(ned,nen,*)

         integer*4:: i, j, k
         integer*4:: m, minimo, num

         do k=1,numel
            minimo = neq
            do j=1,nen
               do i=1,ned
                  num = lm(i,j,k)
                  if (num.gt.0) minimo = min(minimo,num)  ! min= min0(min,num) 
               enddo
            enddo

            do j=1,nen
               do i=1,ned
                  num = lm(i,j,k)
                  if (num.gt.0) then
                     m = num - minimo
                     if (m.gt.idiag(num)) then
                        idiag(num) = m
                     endif
                  endif
               enddo
            enddo
         enddo
         
         return
      end subroutine

      !******************************************************************************
      !> Routine to create the sparse matrix pointers CSR
      !> @param estrutSistEq_ -
      !> @param           nsd -
      !> @param   conectsElem -
      !> @param listaDosElems - 
      !> @param   numConexoes -
      !> @param           nen -
      !> @param     nVizinMax -
      !> @param      simetria -
      subroutine criarPonteirosMatEsparsa_CSR(estrutSistEq_, nsd, conectsElem, listaDosElems, &
               numConexoes, nen, nVizinMax, simetria)

         implicit none

         type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
        
         integer*4,  intent(in) :: nsd, numConexoes, nen, nVizinMax
         logical, intent(in) :: simetria
         integer*4 :: conectsElem(nen,*), listaDosElems(nVizinMax,*)

         write(*,*) " em subroutine criarPonteirosMatEsparsa_CSR(nsd, ndof, neq, num"

         ! if(.not.associated(LMstencilEq_))
         allocate(estrutSistEq_%LMstencilEq(estrutSistEq_%neq,estrutSistEq_%numCoefPorLinha)); 
         estrutSistEq_%LMstencilEq=0
         
         print*, estrutSistEq_%numCoefPorLinha
         call montarLmStencilNodal_CSR (estrutSistEq_%LMstencilEq,listaDosElems, estrutSistEq_%id, &
                  conectsElem, estrutSistEq_%numCoefPorLinha, estrutSistEq_%ndof, numConexoes, nen, &
                  estrutSistEq_%nVizinMax, estrutSistEq_%neq, simetria)

         ! if(.not.associated(LMstencilEq_))
         ! write(*,*) "+++", associated(LMstencilEq_)
         ! allocate(Ap_(neq_+1));    Ap_=0
         ! write(*,'(a, 20i3)')"+++", size( Ap_),  Ap_ (:) ;  
         call montarPonteiroAp_CSR(estrutSistEq_%Ap, estrutSistEq_%LMstencilEq, &
                  estrutSistEq_%numCoefPorLinha, estrutSistEq_%neq, estrutSistEq_%nalhs)

         allocate(estrutSistEq_%Ai(estrutSistEq_%nalhs)); estrutSistEq_%Ai=0
         call montarPonteiroAi_CSR(estrutSistEq_%Ai, estrutSistEq_%LMstencilEq, &
                  estrutSistEq_%numCoefPorLinha, estrutSistEq_%neq, estrutSistEq_%nalhs)

         contains

         !**** new *************************************************************
         !> Routine for creating the stencil of the LM matrix
         !> @param     LMstencilEq -
         !> @param   listaDosElems -
         !> @param              id -
         !> @param     conectsElem - 
         !> @param numCoefPorLinha -
         !> @param            ndof -
         !> @param     numConexoes -
         !> @param             nen -
         !> @param       nVizinMax -
         !> @param             neq -
         !> @param        simetria -
         subroutine montarLmStencilNodal_CSR(LMstencilEq, listaDosElems, id, &
                  conectsElem, numCoefPorLinha, ndof, numConexoes, nen, nVizinMax, neq, simetria) 

            implicit none

            integer*4,  intent(in)    :: neq, nVizinMax
            integer*4,  intent(in)    :: numCoefPorLinha, numConexoes,ndof,nen
            integer*4,  intent(inout) :: LMstencilEq(neq,numCoefPorLinha)
            integer*4,  intent(in)    :: listaDosElems(nVizinMax,*), conectsElem(nen,*)
            integer*4,  intent(in)    :: id(ndof,*)
            logical    :: simetria
            ! logical, intent(in)    :: simetria

            integer*4 :: nc, cont, nel, i, j, k, l, m,  numEq, dir
            logical :: propriaEq, difZero
            integer*4 :: numCoefAux, numViz, menorQueNumEq(ndof)
            integer*4,  allocatable :: LMstencilEqAux(:)
            integer*4 :: shift
      
            print*, nVizinMax

            ! simetria=.false.
            write(*,*) " ++ subroutine montarLmStencilNodal_CSR(...., numCoefPorLinha = ", numCoefPorLinha     
            shift=10
            numCoefAux=numCoefPorLinha+shift 
            allocate(LMstencilEqAux(numCoefAux))

            numViz=nVizinMax
            numEq=0
            LMstencilEq=0
            LMstencilEqAux=0

            do nc=1, numConexoes
               do dir=1, ndof 
   
                  if(id(dir,nc).ne.0) then
                     numEq=numEq+1
                     LMstencilEqAux=0
                     cont=1
                     propriaEq=.false.

                     do k=1, numViz
                        if(listaDosElems(k,nc).ne.0) then
                           nel=listaDosElems(k,nc)

                           do i=1, nen!4!numConexoesPorElem

                              menorQueNumEq=0
                              difZero=.false.
                              do j=1, ndof
                                 if(id(j,conectsElem(i,nel)).ne.0) difZero=.true.
                                 if(id(j,conectsElem(i,nel)).ne.numEq) menorQueNumEq(j)=1
                              enddo

                              if(difZero.eqv..true.) then !condicional verdadeira se algum dos id's for diferentes de zero
                                 if(sum(menorQueNumEq)==ndof) then !condicional verdadeiro se os id's de todas as direcoes nao for o numero da equacao que estah sendo avaliada
                           
                                    !inclui no LMStencil as equacoes que contribuem para a linha da matriz exceto o id que contem a propria equacao
                                    do m=1, ndof
                                       if(simetria.eqv..true.) then
                                          if(id(m,conectsElem(i,nel))>=numEq) then
                                             LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                             cont=cont+1
                                          endif
                                       else
                                          if(id(m,conectsElem(i,nel)).ne.0) then
                                             LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                             cont=cont+1
                                          endif
                                       endif
                                    enddo

                                 else
                                    !inclui no LMStencil as equacoes onde o id contem a propria equacao 
                                    if(propriaEq.eqv..false.) then
                                       do m=1, ndof
                                          propriaEq=.true.
                                          if(simetria.eqv..true.) then
                                             if(id(m,conectsElem(i,nel))>=numEq) then
                                                LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                                cont=cont+1
                                             endif
                                          else
                                             if(id(m,conectsElem(i,nel)).ne.0) then
                                                LMstencilEqAux(cont)=id(m,conectsElem(i,nel))
                                                cont=cont+1
                                             endif
                                          endif
                                       enddo
                                    endif

                                 end if
                              end if
                           end do
                        end if
                     end do !k

                     call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)
                     do i=1, numCoefAux-1
                        if(LMstencilEqAux(i)==LMstencilEqAux(i+1)) LMstencilEqAux(i)=0
                     end do

                     call ordenarLMstencil(LMstencilEqAux(:),numCoefAux)
                     LMstencilEq(numEq,1:numCoefPorLinha)=LMstencilEqAux(shift+1:numCoefAux)
                     ! write(*,'(a,i0,a,81(1x,i0))'), "LmStencil ", numEq, " ->",LMstencilEq(numEq,1:numCoefPorLinha)
                  end if

               end do !dir
            end do !gl

            deallocate(LMstencilEqAux)
         end subroutine montarLmStencilNodal_CSR

         !----------------------------------------------------------------------
         !> Routine for sorting the stencil
         !> @param     LMstencilEq -
         !> @param numCoefPorLinha -
         subroutine ordenarLMstencil(LMstencilEq,numCoefPorLinha)

            implicit none

            integer*4,  intent(in)    :: numCoefPorLinha      
            integer*4,  intent(inout) :: LMstencilEq(numCoefPorLinha)

            integer*4 :: tmp
            integer*8 :: menorEq, n, nn 

            do n = 1, numCoefPorLinha
               menorEq=n
               do nn = n+1, numCoefPorLinha
                  if(LMstencilEq(nn)<LMstencilEq(menorEq)) menorEq = nn
               end do

               if(n == menorEq) cycle
               tmp                  = LMstencilEq(n)
               LMstencilEq(n)       = LMstencilEq(menorEq)
               LMstencilEq(menorEq) = tmp
                
            enddo
         end subroutine ordenarLMstencil

         !**** new *************************************************************
         !> Routine to assemble the pointeros Ap do CSR
         !> @param              Ap -
         !> @param     LMstencilEq -
         !> @param numCoefPorLinha -
         !> @param             neq -
         !> @param        nonzeros -
         subroutine montarPonteiroAp_CSR (Ap, LMstencilEq, numCoefPorLinha, neq, nonzeros)

            implicit none

            integer*4 :: Ap(:)
            integer*4 :: nonzeros, neq
            integer*4 :: numCoefPorLinha
            integer*4 :: LMstencilEq(neq,numCoefPorLinha)
            integer*8 :: l

            ! Montando Ap                  
            call montarListaPonteiros(Ap, LMstencilEq,neq,numCoefPorLinha)

            ! Contando os valores nao nulos
            nonzeros=0
            do l=2, neq+1
               nonzeros=nonzeros+(Ap(l)-Ap(l-1))
            end do

         end subroutine montarPonteiroAp_CSR

         !**** new *************************************************************
         !> Routine to assemble the pointeros Ai do CSR
         !> @param              Ai -
         !> @param     LMstencilEq -
         !> @param numCoefPorLinha -
         !> @param             neq -
         !> @param        nonzeros -
         subroutine montarPonteiroAi_CSR (Ai, LMstencilEq, numCoefPorLinha, neq, nonzeros)

            implicit none
            integer*4,  intent(out) :: Ai(:)
            integer*4 :: nonzeros, neq
            integer*4 :: numCoefPorLinha
            integer*4 :: LMstencilEq(neq,numCoefPorLinha)

            call montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)
         end subroutine montarPonteiroAi_CSR

         !**** new *************************************************************
         !> Routine to assemble the pointeros 
         !> @param              Ap -
         !> @param     LMstencilEq -
         !> @param numCoefPorLinha -
         !> @param             neq -
         !> @param        nonzeros -
         subroutine montarListaPonteiros(Ap, LMstencilEq, neq, numCoefPorLinha)

            implicit none
            integer*4,  intent(in) :: neq
            integer*4,  intent(out) :: Ap(neq+1)
            integer*4,  intent(in) :: numCoefPorLinha
            integer*4,  intent(in) :: LMstencilEq(neq,numCoefPorLinha)
            integer*4 :: LMstencilEqTemp(0:numCoefPorLinha)

            integer*8 :: n, l
            integer*4 :: posPonteiro, contPonteiro

            posPonteiro=0
            contPonteiro=0
            ! write(*,'(13i3)') Ap (:); 
            do l=1, neq
               ! write(*,'(9i3)') lmstencilEQ (l,:); 
               posPonteiro=posPonteiro+1
               Ap(posPonteiro)=contPonteiro+1
               LMstencilEqTemp=0
               LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(l,:)
               do n = 1, numCoefPorLinha
                  if(LMstencilEqTemp(n) == LMstencilEqTemp(n-1)) cycle
                  contPonteiro= contPonteiro+1
               enddo    

               if(l==neq) then
                  posPonteiro=posPonteiro+1
                  Ap(posPonteiro)=contPonteiro+1
               end if
            end do
         end subroutine montarListaPonteiros

         !**** new *************************************************************
         !> Routine to assemble a list of indices
         !> @param              Ai -
         !> @param     LMstencilEq -
         !> @param             neq -
         !> @param        nonzeros -
         !> @param numCoefPorLinha -
         subroutine montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)

            implicit none
            integer*4,  intent(in)    :: nonzeros, neq
            integer*4,  intent(inout) :: Ai(nonzeros)
            integer*4,  intent(in)    :: numCoefPorLinha
            integer*4,  intent(in)    :: LMstencilEq(neq,numCoefPorLinha)

            integer*4 :: LMstencilEqTemp(0:numCoefPorLinha)
            integer*4 :: i, n, posColunas

            posColunas=0
            do i=1, neq
               LMstencilEqTemp=0
               LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(i,:)
               do n = 1, numCoefPorLinha
                  if(LMstencilEqTemp(n).ne.LMstencilEqTemp(n-1).and.LMstencilEqTemp(n).ne.0 ) then
                     posColunas=posColunas+1
                     Ai(posColunas)= LMstencilEqTemp(n)
                  end if
               enddo
            end do

         end subroutine montarListaIndices

      end subroutine criarPonteirosMatEsparsa_CSR

      !=======================================================================
      !> Routine for writing the matrix system in CSR format
      !> @param     alhs -
      !> @param     brhs -
      !> @param       Ap -
      !> @param       Ai -
      !> @param nonzeros -
      !> @param      neq -
      !> @param  nomeArq -
      subroutine escreverSistemaAlgCSRemMTX(alhs, brhs, Ap, Ai,  nonzeros, neq, nomeArq)

         implicit none 
         real*8,  intent(in)   :: Alhs(:), brhs(:)
         integer*4,  intent(in)  :: Ap(:), Ai(:)
         integer*4,  intent(in)  :: neq, nonzeros
         character (len=*)  :: nomeArq

         integer*8 :: i, j, k
         integer*8 :: luSist = 1836 
         character(len=40), parameter :: formatoEscritaA='(2(i0,1x),e23.16)'
         character(len=40), parameter :: formatoEscritaB='(e23.16)'
         character(len=40), parameter :: formatoEscritaC='(e23.16,a)'
         real*8 :: t1, t2

         write(*,*) " em subroutine escreverSistemaAlgCSRemMTX(alhs, brhs, Ap, Ai, ... "
         call timing(t1)
         open(file=nomeArq, unit=luSist) 

         write(luSist,'(a)')'%% matriz A de coeficientes reais simetrica positiva definida ' 
         write(luSist,'(a)')'% produzida pelo metodo classico de galerkin para o metodo de elementos finitos '
         write(luSist,'(a)')'% armazenamento esparso CSR '
         write(luSist,'(a,3(i0,a))' )  '%  matriz ',neq, 'X',neq, ' com ', nonzeros, ' coefs diferentes de zero' 
         write(luSist,'(  3(i0,1x))' )  neq, neq, nonzeros 
         ! write(*,*) size(alhs), size(brhs), size(Ap), size(Ai)  
         k = 1
         do i = 1, neq
            ! write(luA, *) Ap(i) ,  Ap(i+1)
            do j = Ap(i) ,  Ap(i+1) - 1  
               ! write(luSist, formatoEscritaA   ) i, Ai(j), alhs(k)
               ! write(luA+luAux, '( 3(i10,2x), e20.10, 2x,i5)'     ) i, j,  Ai(j), alhs(k), k
               k = k + 1
            end do
         end do
         ! write(luSist,'( (a,i0,a))') '% lado direito com ',neq,' elementos' 
         ! write(luSist,'(  3(i10))' )  neq 
         i = 1
         write(luSist, formatoEscritaC )  brhs(i), "  BRHS"
         do i = 2, neq
            write(luSist, formatoEscritaB )  brhs(i)
         end do
         close(luSist)
         call timing(t2)
         write(*,*) " tempo de escrita =", t2 - t1

      end subroutine escreverSistemaAlgCSRemMTX

      !=======================================================================
      !> Routine for reading the matrix system in CSR format
      !> @param     alhs -
      !> @param     brhs -
      !> @param       Ap -
      !> @param       Ai -
      !> @param      neq - 
      !> @param nonzeros -
      !> @param  nomeArq -
      subroutine lerSistemaAlgMTXemCSR(alhs, brhs, Ap, Ai, neq, nonzeros,  nomeArq) 

         implicit none 
         real*8   ,intent(inout) :: alhs(:) ! matriz do sistema 
         real*8   ,intent(inout) :: brhs(:) ! lado direito do sistema
         integer*4,intent(inout) :: Ap(:) ! vetor apontadores para  inicio de uma linha
         integer*4,intent(inout) :: Ai(:) ! vetor com as colunas dos coefcientes de alhs
         integer*8,intent(inout) :: nonzeros
         integer*8,intent(inout) :: neq 
         character(len=*), intent(in) :: nomeArq

         integer   :: luSist = 1000, luSistLido = 2000
         integer*4 :: k 
         integer*8 :: i, j, ieq, ieqAnterior
         character(len=100) :: label 
         character(len=40), parameter :: formatoEscritaA='(2(i0,1x),e15.8)'
         character(len=40), parameter :: formatoEscritaB='(1(i0,1x),e15.8)'

         write(*,*) " em subroutine lerSistemaMTXemCSR(alhs, brhs, Ap, Ai, neq, nonzeros,",nomeArq,")"
         write(*,*) " FUNCIONANDO PARA MATRIZES SIMETRICAS COM ELEMENTOS FORNECIDOS POR LINHAS"
         open(file=nomeArq, unit=luSist) 

         read(luSist,'(a)' ) label; write(*,* ) label 
         read(luSist,'(a)' ) label; write(*,* ) label 
         read(luSist,'(a)' ) label; write(*,* ) label 
         read(luSist,'(a)' ) label; write(*,* ) label 
         read(luSist,* )  neq, neq, nonzeros
         write(*,'(a,3(i10))' ) "+++", neq, neq, nonzeros

         ! o arquivo que será lido estah escrito com os coeficientes de uma linha 
         ! considerando matrix simetrica
         ! as linhas estao em ordem crescente de 1 ateh neq (num de equacoes) 

         Ap(:) = 1
         k = 0; ieqAnterior = 0
         do i = 1, nonzeros
            k = k + 1 ! contador de nao zeros lidos
            read(luSist, * ) ieq,  Ai(k), alhs(k)
            if(ieq > ieqAnterior) Ap(ieq) = k  
            ieqAnterior = ieq
         end do
         k = k + 1
         Ap(ieq+1) = k

         write(*,'( (a,i0,a))') 'lado direito com ',neq,' elementos' 
         do i = 1, neq
            read(luSist, * ) brhs(i)
         end do

         write(luSistLido,'(  3(i0,1x))' )  neq, neq, nonzeros 
         k = 1
         do i = 1, neq
            do j = Ap(i) ,  Ap(i+1) - 1  
               write(luSistLido, formatoEscritaA   ) i, Ai(j), alhs(k)
               k = k + 1
            end do
         end do
         do i = 1, neq
            write(luSistLido, formatoEscritaB ) i, brhs(i)
         end do

         close(luSist); 
      end subroutine lerSistemaAlgMTXemCSR
end module mEstruturasDadosSistEq
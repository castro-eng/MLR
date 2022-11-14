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
!         Petropolis, 02.2014
!=================================================================================
module mFluxo
   ! Modulo contendo variaveis e operações referentes ao fluxo

   use mEstruturasDadosSistEq
   implicit none

   type (estruturasArmazenamentoSistemaEq) :: estrutSistEqF
   type (estruturasArmazenamentoSistemaEq) :: estrutSistEqFluxoNormal
   real*8, private :: vm_x_dir, vm_x_esq, vm_y_atras, vm_y_frente, vm_z_inf, vm_z_sup
   real*8, private :: gradPx, gradPy, gradPz, vazao(6)

   real*8, allocatable :: velocidadeBordas(:,:), velocidadeBordasFratura(:,:), velocidadeBordasFalhas(:,:,:)
   real*8, allocatable :: vm_integro(:,:), vm_reduzido(:,:,:)
   real*8, allocatable :: vm_bloco(:,:,:,:), vm_falha(:,:,:,:)

   real*8, allocatable :: vazaoAnt(:), prodAcum(:)

   ! TO DO -> vm_falha - old, variavel descartada.

   contains
   !=======================================================================
   !-------------- ROTINAS PARA O CONTROLE DO MODULO ----------------------
   !=======================================================================
   subroutine calcularVelocidade()
      ! Rotina para determinação das velocidades no domínio.
      ! A determinação da velocidade é realizada utilizando um pós-processamento
      ! para o qual em posse das pressões (p= E Nd), temos a velocidade como 
      ! função da lei de Darcy (u=K dp).

      use mMalha, only: numel, nen, numnp, nsd, numelBordas
      use mMalha, only: numelFratura, numnpFratura, nenFratura
      use mMalha, only: conecNodaisElem, conecNodaisFratura, conecNodaisBordas
      use mGlobaisArranjos,  only: tiposElemBordas
      use mGlobaisEscalares, only: npintFratura!, simulandoReferencia
      use mPotencial,        only: estrutSistEqP

      real*8, dimension(nsd,numnpFratura,(estrutSistEqP%ndof-2)) :: vm_mediaF
      
      ! //TODO:: retirar o uso da estrutura estrutSistEqF para armazenar a velocidade nos nós
      
      ! Calculo da velocidade nos nós e determinação das velocidades nos bordos
      call calcVelocidadeDominioIntegro(estrutSistEqF%u,estrutSistEqF%ndof)
      
      ! Descontinuidades - Falha
      if (numelFratura>0) then
         ! Calculo da velocidade nos nós e determinação das velocidades nos bordos
         call calcVelocidadeDominioReduzido(velocidadeBordasFalhas, estrutSistEqP%ndof,&
              (estrutSistEqP%ndof-2))

         ! Corrigir velocidades no domínio integro
         ! call corrigirVelocidadeDominioIntegro(estrutSistEqF%u,estrutSistEqF%ndof)
      end if
      
      ! Maciço Rochoso
      ! 1 - Determinação das velocidades nos pontos de integração
      ! call calcVelocElemBloco (estrutSistEqF%u, conecNodaisElem, &
      !                          numel, numnp, nsd, nen, estrutSistEqF%ndof)
      ! 2 - Determinação da velocidades nos nós
      ! call calcularMediasNos(vm_bloco, estrutSistEqF%u, conecNodaisElem, &
      !                        numel, numnp, nsd, nen, estrutSistEqF%ndof, 1)
      ! 3 - Determinação das velocidades nas bordas
      ! call calcularVelocBordas(estrutSistEqF%u, velocidadeBordas, &
      !                          conecNodaisBordas, tiposElemBordas, &
      !                          numelBordas, numnp, nsd, nenFratura)
      ! Descontinuidades - Falha
      ! if (numelFratura>0) then
      !    ! 1 - Determinação das velocidades nos pontos de integração e na bordas
      !    call calcVelocElemFalha (estrutSistEqF%u, velocidadeBordasFalhas, conecNodaisFratura, &
      !    numnp, numelFratura, nsd, nenFratura, estrutSistEqP%ndof, npintFratura, &
      !    (estrutSistEqP%ndof-2))
      !    ! 2 - Determinação da velocidades nos nós
      !    call calcularMediasNosFalha(vm_falha, vm_mediaF, conecNodaisFratura, numelFratura, &
      !    numnpFratura, nsd, nenFratura, estrutSistEqP%ndof, (estrutSistEqP%ndof-2))
      ! end if

   end subroutine
   !=======================================================================
   !-------------- ROTINA PARA VELOCIDADE NO DOMÍNIO INTEGRO --------------
   !=======================================================================
   subroutine calcVelocidadeDominioIntegro(u,ndof)
      ! Rotina para determinação da velocidade no domínio integro dado nos nós
      ! Armazena as velocidades na matriz "vm_integro"
      ! - u é o campo de pressão dados nos pontos nodais
      use mGlobaisEscalares, only: npint, nrowsh, viscosidade
      use mGlobaisArranjos,  only: tiposElemBordas, campoPermeabilidade
      use mMalha,            only: x, numel, nen, numnp, nsd
      use mMalha,            only: numelFratura, nenFratura, numelBordas
      use mMalha,            only: conecNodaisElem, conecNodaisBordas, normalElemBorda
      use mMalha,            only: local, normal
      use mFuncoesDeForma,   only: shlt, shlq, SHLQEN, SHLTEN, shlqen3D,  shlq3d, shlten3D, shlt3D
      use mFuncoesDeForma,   only: oneshg, SHGQ, shg3d
      use mPotencial,        only: estrutSistEqP

      IMPLICIT NONE

      integer, intent(in) :: ndof
      real*8, intent(out) :: u(ndof,numnp)
      
      LOGICAL :: QUAD

      integer :: I, J, K, L, NEL, noGlobal, m, numOcorrenciaNos(numnp)
      integer :: nelBordas, noBordas, tipoBorda
      integer :: dof

      real*8 :: gradP
      real*8 :: xl(nsd,nen), dl(estrutSistEqP%ndof,nen)
      real*8 :: det(nen), W(npint), di(nrowsh)
      real*8, dimension(nrowsh,nen,npint) :: SHL, SHG
      real*8  :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1)
      real*8  :: vetorParaNo(nsd), normalNoLocal(nsd,nen)

      real*8  :: erro;

      QUAD=.TRUE.

      ! GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
      w=0d0
      shl=0d0
      wCentroElem=0d0
      shlCentroElem=0d0

      if(nen==3) then                        ! Elemento(triangulo) - Espaco()
         CALL SHLTEN(SHL,NEN)
         call shlt(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==2) then        ! Elemento(quadrado) - Espaco(2D)
         CALL SHLQEN(SHL,NEN)
         call shlq(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==3) then        ! Elemento(tetraedro) - Espaco(3D)
         CALL SHLTEN3D(SHL,NEN)
         call shlt3D(shlCentroElem,wCentroElem,1,nen)
      else if(nen==8) then                   ! Elemento(cubo) - Espaco(3D)
         CALL SHLQEN3D(SHL,NEN)
         call shlq3d(shlCentroElem,wCentroElem,1,nen)
      end if

      !.... BEGIN LOOP OVER THE MESH ELEMENTS
      vm_integro(:,:)=0.d0
      numOcorrenciaNos(:)=0
      do nel=1,numel

         !....    LOCALIZE UNKNOWNS AND COORDINATES
         ! Armazena os valores das posições (x) do elemento em xl(nsd,nen)
         ! Armazena os valores das pressoes (estrutSistEqP%u) do elemento em dl(ndof,nen)
         CALL LOCAL(conecNodaisElem(1,NEL),x,xl,NEN,NSD,NSD)
         CALL LOCAL(conecNodaisElem(1,NEL),estrutSistEqP%u,DL,NEN,estrutSistEqP%ndof,estrutSistEqP%ndof)

         if (numelFratura > 0) then
            call local(conecNodaisElem(:,nel), normal, normalNoLocal, nen, nsd, nsd)

            do i=1,nsd
               xCentroElem(i) = dot_product(xl(i,:), shlCentroElem(nrowsh,:,1))
            end do
         end if

         !        1D-2D
         !....    EVALUATE GLOBAL SHAPE FUNCTION
         if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,1,nel,1)   ! Elemento(linha) - Espaco()
         if(nen==3) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)      ! Elemento(triangulo) - Espaco()
         if(nen==4.and.nsd==2) call shgq  (xl,det,shl,shg,npint,nel,quad,nen) ! Elemento(quadrado) - Espaco(2D)
         if(nen==4.and.nsd==3) call shg3d (xl,det,shl,shg,npint,nel,nen)      ! Elemento(tetraedro) - Espaco(3D)
         if(nen==8) call shg3d (xl,det,shl,shg,npint,nel,nen)           ! Elemento(cubo) - Espaco(3D)

         !....    COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
         !....    BEGIN LOOP OVER ELEMENT NODES
         do L=1,nen
            noGlobal = conecNodaisElem(L,nel)

            ! COMPUTE THE DERIVATIVES IN FUNCTION OF THE SPACE DIMENSION
            do k=1,nsd

               ! COMPUTE FIRST DERIVATIVES OF NODAL DISPLACEMENT
               gradP=0.D0
               do j=1,nen
                  ! Para o caso do elemento na interface com a descontinuidade
                  dof = 1
                  if (numelFratura > 0) then
                     if ((estrutSistEqP%lm(2,J,nel) > 0) .or. (estrutSistEqP%lm(3,J,nel) > 0)) then
                        vetorParaNo(:) = xl(:,J) - xCentroElem(:)                   
                        if (dot_product(normalNoLocal(:,J), vetorParaNo(:)) < 0d0) then
                           dof = 2
                        end if
                     end if
                  end if
                  gradP=gradP + shg(k,J,L)*DL(dof,J)
               enddo
               vm_integro(k,noGlobal)=vm_integro(k,noGlobal)-(campoPermeabilidade(k,nel)/viscosidade)*gradP
            end do
            numOcorrenciaNos(noGlobal)=numOcorrenciaNos(noGlobal)+1
         end do
      enddo ! END OF LOOP OVER MESH ELEMENTS
      
      ! COMPUTING AVERAGE velocity AT NODAL POINTS
      do j=1,numnp
         do i=1,nsd
            ! Obs.. a variável numOcorrenciaNos fica repetida pelo numero de dimensão do espaço analisado
            vm_integro(i,j)=(vm_integro(i,j)/(numOcorrenciaNos(j)))
            u(i,j)=vm_integro(i,j)
         end do
      end do

      ! Calculando velocidadeBordas que eh usada no calculo da vazao
      velocidadeBordas(:,:)=0.d0
      do i=1, numelBordas  !------------------- Loop nos elementos pertencentes a borda
         do j=1, nenFratura   !---------------- Loop nos nós dos elementos
            noBordas=conecNodaisBordas(j, i) !- Ponteiro para no do elemento da borda
            tipoBorda = tiposElemBordas(i)   !- Ponteiro para o tipo de borda
            if(tipoBorda==1) then   !---------- Armazenagem do velocidade da borda na variavel 
               velocidadeBordas(j,i)=-vm_integro(2,noBordas)
            elseif(tipoBorda==2) then 
               velocidadeBordas(j,i)=vm_integro(1,noBordas)
            elseif(tipoBorda==3) then 
               velocidadeBordas(j,i)=vm_integro(2,noBordas)
            elseif(tipoBorda==4) then 
               velocidadeBordas(j,i)=-vm_integro(1,noBordas)
            elseif(tipoBorda==5) then 
               velocidadeBordas(j,i)=-vm_integro(3,noBordas)
            elseif(tipoBorda==6) then 
               velocidadeBordas(j,i)=vm_integro(3,noBordas)
            endif
         end do          
      end do

      RETURN      
   end subroutine
   !=======================================================================
   !------------- ROTINA PARA VELOCIDADE NO DOMÍNIO REDUZIDO --------------
   !=======================================================================
   subroutine calcVelocidadeDominioReduzido(velbord,ndof,nlayer)
      ! Rotina para determinação da velocidade no domínio reduzido dado nos nós
      ! Armazena as velocidades na matriz "vm_reduzido"

      ! A metodologia para calculo da velocidade no domínio reduzido tem a 
      ! velocidade tangencial calculada pela função de interpolação do 
      ! elemento reduzido, e a velocidade normal calculada pelo salto da pressão

      use mGlobaisEscalares, only: nrowsh, npintFratura, viscosidade
      use mGlobaisArranjos,  only: c, matFratura
      use mMalha,            only: x, numnp, nsd, normal, local
      use mMalha,            only: numelFratura, nenFratura, numnpFratura
      use mMalha,            only: conecNodaisFratura, conecNodaisReduzido
      use mMalha,            only: numelBordas, conecNodaisBordas, normalElemBorda
      use mFuncoesDeForma,   only: shlt, shlq, shl1D_intPoints
      use mFuncoesDeForma,   only: shgCodimOneSurface_points
      use mPotencial,        only: estrutSistEqP
      use mMLayer,           only: LayerMaterial

      IMPLICIT NONE

      integer, intent(in) :: ndof, nlayer
      real*8, intent(out) :: velbord(nenFratura, numelBordas, nlayer)

      ! u -> estrutSistEqF%u
      ! velbord -> velocidadeBordasFalhas
      ! conec -> conecNodaisFratura
      
      !    real*8, intent(in)  :: u(ndof,numnp)
      !    integer, intent(in) :: conecNodaisFratura(nenFratura, numelFratura)
      !    real*8, intent(out) :: velbord(nenFratura, numelBordas, nlayer)

      LOGICAL :: QUAD
      integer :: I, J, K, L, NEL, SP, noGlobal, m, nelBordas, noBordas, tipoBorda
      integer :: dof, ncamadas, numOcorrenciaNos(numnpFratura)
      real*8  :: salto,gradP
      real*8 :: velLocal(nsd)

      real*8 :: shlCentroElem(nrowsh,nenFratura,1), wCentroElem(1)

      real*8, dimension(nsd) :: xCentroElem, vetorParaNo
      real*8, dimension(npintFratura) :: det, W
      real*8, dimension(nsd,nenFratura) :: xl, normalNoLocal
      real*8, dimension(ndof,nenFratura) :: dl
      real*8 :: SHL(nsd,nenFratura,npintFratura), SHG(nrowsh,nenFratura,npintFratura) 

      real*8, dimension(:), allocatable :: layKN, layDh, layKT
      ! real*8, dimension(:,:), allocatable :: layKT

      QUAD=.TRUE.

      ! GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
      shl=0d0
      if(nenFratura==2) then        ! Elemento(linha) - Espaco(2D)
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then   ! Elemento(triangulo) - Espaco(3D)
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then   ! Elemento(quadrado) - Espaco(3D)
         call shlq(shl,w,npintFratura,nenFratura)
      end if

      !.... BEGIN LOOP OVER THE MESH ELEMENTS
      vm_reduzido(:,:,:)=0d0
      numOcorrenciaNos(:)=0
      DO NEL=1,numelFratura

         !.... LOCALIZE UNKNOWNS AND COORDINATES
         ! Armazena os valores das posições (x) do elemento em xl(nsd,nenFratura)
         ! Armazena os valores das pressoes (estrutSistEqP%u) do elemento em dl(ndof,nenFratura)
         ! Armazena os vetor normal (normal) do elemento em normalNoLocal(nsd,nenFratura)
         CALL LOCAL(conecNodaisFratura(1,NEL),x,XL,nenFratura,nsd,nsd)
         CALL LOCAL(conecNodaisFratura(1,NEL),estrutSistEqP%u,DL,nenFratura,ndof,ndof)
         CALL LOCAL(conecNodaisFratura(:,NEL), normal, normalNoLocal, nenFratura, nsd, nsd)

         ! print*, "conecNodaisFratura->", conecNodaisFratura(1,NEL)
         ! print*, "nenFratura->", nenFratura
         ! print*, "ndof->", ndof
         
         ! do i=1,ndof
         !    do j=1,nenFratura
         !       WRITE(*,'("dl(",i3,",",i3,")=",e12.5)') i,j,dl(i,j)
         !    enddo
         ! enddo

         ! do i=1,nrowsh
         !    do j=1,nenFratura
         !       WRITE(*,'("shg(",i3,",",i3,",",i3,")=",e12.5)') i,j,1,shg(i,j,1)
         !    enddo
         ! enddo
         ! print*, "------"
   
         !... Properties
         m = matFratura(nel)
         ncamadas = LayerMaterial(m,1)
         allocate(layDh(ncamadas))
         allocate(layKN(ncamadas))
         allocate(layKT(ncamadas))
         do k=1, ncamadas
            layDh(k) = c(2,LayerMaterial(m,k+1))
            layKN(k) = c(3,LayerMaterial(m,k+1))/viscosidade
            layKT(k) = c(4,LayerMaterial(m,k+1))/viscosidade
            ! if (nsd.eq.3) layKT(k,2) = c(5,LayerMaterial(m,k+1))/viscosidade
         enddo

         ! print*,'ncamadas',ncamadas

         ! 1D-2D
         !.... EVALUATE GLOBAL SHAPE FUNCTION
         shg=0d0
         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)

         !.... COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
         !.... BEGIN LOOP OVER ELEMENT NODES
         DO L=1,nenFratura
            noGlobal=conecNodaisReduzido(L,NEL)
            
            !.... BEGIN LOOP OVER EACH LAYER FOR THE ELEMENT
            do k=1, ncamadas
               
               ! COMPUTE THE DERIVATIVES IN FUNCTION OF THE SPACE DIMENSION
               do SP=1,nsd
                  gradP=0d0

                  ! COMPUTE TANGENCIAL DERIVATIVES OF THE NODAL DISPLACEMENT
                  do j=1,nenFratura
                     if ((estrutSistEqP%lmFratura(2,J,nel)==0).or.(estrutSistEqP%lmFratura(3,J,nel)==0))then
                        dof=1
                     else
                        dof=2+k
                     endif
                     gradP=gradP+SHG(sp,j,L)*DL(dof,j)
                  enddo
                  ! if(k==3) print*, "grad", gradP
                  
                  ! COMPUTE NORMAL JUMP OF THE NODAL DISPLACEMENT
                  if(ncamadas.eq.1)then
                     salto=DL(2,l)-DL(1,l)
                  else
                     if(k==1)then
                        salto=DL(4,l)-DL(1,l)
                     elseif(k==ncamadas)then
                        salto=DL(2,l)-DL((k+1),l)
                     else
                        salto=DL((k+3),l)-DL((k+1),l)
                     endif
                  endif
                  ! old
                  ! if (k==1)then
                  !    if (ncamadas.eq.1) then
                  !       salto=DL(2,l)-DL(1,l)
                  !    else
                  !       salto=DL(4,l)-DL(1,l)
                  !       ! salto=(DL(4,l)+DL(3,l))/2.d0-DL(1,l)
                  !    endif
                  ! elseif (k==ncamadas) then
                  !    salto=DL(2,l)-DL((k+1),l)
                  ! else
                  !    salto=DL((k+3),l)-DL((k+1),l)
                  !    ! salto=DL((k+2),l)-DL((k+1),l)
                  ! endif

                  velLocal(sp)=-layKN(sp)*salto*normalNoLocal(SP,l)/layDh(k)
                  velLocal(sp)=velLocal(sp)-layKT(k)*gradP
                  vm_reduzido(sp,noGlobal,k)=vm_reduzido(sp,noGlobal,k)+velLocal(sp)
               enddo
               numOcorrenciaNos(noGlobal)=numOcorrenciaNos(noGlobal)+1
               
               ! Determinação das velocidades nas bordas
               do nelBordas=1, numelBordas
                  do J=1, nenFratura
                     noBordas=conecNodaisBordas(J, nelBordas)
                     if(noBordas==noGlobal) then
                        velbord(j,nelBordas,k)=dot_product(velLocal(:),normalElemBorda(:,nelBordas))
                     endif
                  enddo
               enddo
            enddo
         enddo
         DEALLOCATE(layDh)
         DEALLOCATE(layKN)
         DEALLOCATE(layKT)
      enddo

      ! COMPUTING AVERAGE velocity AT NODAL POINTS
      do j=1,numnpFratura
         do i=1,nsd
            do sp=1,nlayer
               ! Obs.. a variável numOcorrenciaNos fica repetida pelo numero de dimensão do espaço analisado
               vm_reduzido(i,j,sp)=(vm_reduzido(i,j,sp)/(numOcorrenciaNos(j)))
               ! WRITE(*,'("Vm_reduzido(",i3,",",i3,",",i3,")=",e12.5)') i, j, sp, &
                     ! vm_reduzido(i,j,sp)
            end do
         end do
      end do

   end subroutine
   !=======================================================================
   !------------- ROTINA PARA VELOCIDADE NO DOMÍNIO REDUZIDO --------------
   !=======================================================================
   subroutine corrigirVelocidadeDominioIntegro(u,ndof)
      ! Rotina para corrigir as velocidades nodais nos nós com elementos reduzidos
      ! - u é o campo de pressão dados nos pontos nodais

      use mMalha,            only: nsd, numnp
      use mMalha,            only: numelFratura, nenFratura, numnpFratura
      use mMalha,            only: conecNodaisFratura, conecNodaisReduzido
      use mGlobaisArranjos,  only: c,matFratura
      use mMLayer,           only: LayerMaterial

      integer, INTENT(IN) :: ndof
      real*8, intent(out) :: u(ndof,numnp)

      integer :: ielm, inen, isd, icamadas, ncamadas, noGlobal, noGlobalR
      integer :: ncamadasNo(numnpFratura)
      DOUBLE PRECISION :: espCamada,espDiscont

      ! Loop em todos os elementos reduzidos
      ncamadasNo(:)=0
      do ielm=1,numelFratura
         
         ! Loop em todos os nós do elemento reduzido
         ncamadas = LayerMaterial(matFratura(ielm),1)
         do inen=1,nenFratura
            
            ! Determinação do nó do elemento em relação ao domínio integro e reduzido
            noGlobal =conecNodaisFratura(inen,ielm)
            noGlobalR=conecNodaisReduzido(inen,ielm)
            
            ! Caso o número de camadas seja maior que já avaliado sobrescrever
            if (ncamadas.GT.ncamadasNo(noGlobalR)) then
               ncamadasNo(noGlobalR)=ncamadas
               do isd=1,nsd
                  ! WRITE(*,'("Vi(",i3,",",i3,")=",e12.5)') isd, noGlobal, vm_integro(isd,noGlobal)
                  vm_integro(isd,noGlobal)=0d0
                  espDiscont=0.d0
                  do icamadas=1,ncamadas
                     espCamada = c(2,LayerMaterial(matFratura(ielm),icamadas+1))
                     vm_integro(isd,noGlobal)=vm_integro(isd,noGlobal)+vm_reduzido(isd,noGlobalR,icamadas)*espCamada
                     espDiscont=espDiscont+espCamada
                     ! WRITE(*,'("Vm_reduzido(",i3,",",i3,",",i3,")=",e12.5)') isd, noGlobalR, icamadas, vm_reduzido(isd,noGlobal,icamadas)
                  enddo
                  vm_integro(isd,noGlobal)=vm_integro(isd,noGlobal)/espDiscont
                  u(isd,noGlobal)=vm_integro(isd,noGlobal)
                  ! WRITE(*,'("Vf(",i3,",",i3,")=",e12.5)') isd, noGlobal, vm_integro(isd,noGlobal)
               enddo               
            endif
         end do
      end do
      
   end subroutine


   !=======================================================================
   !----------------- ROTINAS PARA O BLOCO --------------------------------
   !=======================================================================
   subroutine calcVelocElemBloco(u, conecNodaisElem, numel, numnp, nsd, nen, ndof)
      ! Rotina para determinação da velocidade nos pontos de integração
      ! Armazena as velocidades na matriz "vm_bloco" 
      ! - u é o campo de pressão dados nos pontos nodais
      ! - conecNodaisElem é a matriz de conectividades dos elementos
      ! - numel é o número de elementos para o bloco
      ! - numnp é o número de nós para o bloco
      ! - nsd é o número de dimensão de trabalho (1D,2D,3D)
      ! - nen é o número de nós para cada elemento
      ! - ndof é o número de graus de liberdade por nó
      use mMalha,            only: x, numelFratura
      use mGlobaisEscalares, only: ntype, numat, npint, nicode, iprtin, nrowsh, viscosidade
      use mGlobaisEscalares, only: carregamento, carregamento_ref, saltoPressao
      use mGlobaisArranjos,  only: c, mat, tiposElemBordas, campoPermeabilidade
      use mMalha,            only: local, normal, normalElemBorda, numelBordas, nenFratura, conecNodaisBordas
      use mFuncoesDeForma,   only: oneshl, oneshg, shlt, shlq, SHGQ, SHLQEN, SHLTEN, shlqen3D,  shlq3d, shg3d
      use mFuncoesDeForma,   only: shlten3D, shlt3D
      use mPotencial,        only: estrutSistEqP

      IMPLICIT NONE

      real*8, intent(inout) :: u(ndof,numnp)
      integer, intent(in) :: numel, numnp, nsd, nen, ndof
      integer, intent(in) :: conecNodaisElem(nen, numel)

      LOGICAL :: QUAD

      integer :: I, J, K, L, NEL, noGlobal, m, numOcorrenciaNos(NUMNP)
      integer :: nelBordas, noBordas, tipoBorda
      integer :: dof

      real*8 :: gradP
      real*8 :: xl(nsd,nen), dl(estrutSistEqP%ndof,nen)
      real*8 :: det(nen), W(npint), di(nrowsh)
      real*8, dimension(nrowsh,nen,npint) :: SHL, SHG
      real*8  :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1), vetorParaNo(nsd), normalNoLocal(nsd,nen)

      real*8  :: erro;

      QUAD=.TRUE.

      ! GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
      w=0d0
      shl=0d0
      wCentroElem=0d0
      shlCentroElem=0d0

      if(nen==3) then                        ! Elemento(triangulo) - Espaco()
         CALL SHLTEN(SHL,NEN)
         call shlt(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==2) then        ! Elemento(quadrado) - Espaco(2D)
         CALL SHLQEN(SHL,NEN)
         call shlq(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==3) then        ! Elemento(tetraedro) - Espaco(3D)
         CALL SHLTEN3D(SHL,NEN)
         call shlt3D(shlCentroElem,wCentroElem,1,nen)
      else if(nen==8) then                   ! Elemento(cubo) - Espaco(3D)
         CALL SHLQEN3D(SHL,NEN)
         call shlq3d(shlCentroElem,wCentroElem,1,nen)
      end if       

      !.... BEGIN LOOP OVER THE MESH ELEMENTS
      do nel=1,numel

         !....    LOCALIZE UNKNOWNS AND COORDINATES
         ! Armazena em xl os valores das posições (x) do elemento
         ! Armazena em DL os valores das Pressoes (estrutSistEqP%u) do elemento
         CALL LOCAL(conecNodaisElem(1,NEL),x,xl,NEN,NSD,NSD)
         CALL LOCAL(conecNodaisElem(1,NEL),estrutSistEqP%u,DL,NEN,estrutSistEqP%ndof,estrutSistEqP%ndof)       
         m = mat(nel)

         !        1D-2D
         !....    EVALUATE GLOBAL SHAPE FUNCTION
         if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,1,nel,1)   ! Elemento(linha) - Espaco()
         if(nen==3) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)      ! Elemento(triangulo) - Espaco()
         if(nen==4.and.nsd==2) call shgq  (xl,det,shl,shg,npint,nel,quad,nen) ! Elemento(quadrado) - Espaco(2D)
         if(nen==4.and.nsd==3) call shg3d (xl,det,shl,shg,npint,nel,nen)      ! Elemento(tetraedro) - Espaco(3D)
         if(nen==8) call shg3d (xl,det,shl,shg,npint,nel,nen)           ! Elemento(cubo) - Espaco(3D)

         !....    COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
         !....    BEGIN LOOP OVER ELEMENT NODES
         do L=1,nen
            noGlobal = conecNodaisElem(L,nel)

            ! COMPUTE THE DERIVATIVES IN FUNCTION OF THE SPACE DIMENSION
            do k=1,nsd

               ! COMPUTE FIRST DERIVATIVES OF NODAL DISPLACEMENT
               gradP=0.D0
               do j=1,nen
                  ! Para o caso do elemento na interface com a descontinuidade
                  dof = 1
                  if ( (numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) then
                     if ((estrutSistEqP%lm(2,J,nel) > 0) .or. (estrutSistEqP%lm(3,J,nel) > 0)) then
                        vetorParaNo(:) = xl(:,J) - xCentroElem(:)                   
                        if (dot_product(normalNoLocal(:,J), vetorParaNo(:)) < 0d0) then
                           dof = 2
                        end if
                     end if
                  end if
                  gradP=gradP + shg(k,J,L)*DL(dof,J)
               enddo
               vm_bloco(L,nel,k,1)=-(campoPermeabilidade(k,nel)/viscosidade)*gradP
            end do
         end do
      enddo ! END OF LOOP OVER MESH ELEMENTS

   end subroutine
   !=======================================================================
   subroutine calcularMediasNos(veloc, vm, conect, numel, numnp, nsd, nen, ndof, nlayer)
      ! Rotina para determinação de médias de velocidade nos nós
      ! Armazena as velocidades na matriz "vm_bloco" 
      ! veloc -> matriz de velocidades por elemento nos pontos de integração
      ! vm    -> matriz de velocidades nos nós
      ! conect-> matriz de conectividades dos elementos
      ! numel -> quantidade de elementos
      ! numnp -> quantidade de nós
      ! nen   -> quantidade de nós por elemento
      ! ndof  -> número de graus de liberdade por nó
      ! nsd   -> dimensão do espaço (1D,2D,3D)
      ! nlayer-> quantidade de camadas
      IMPLICIT NONE

      integer, intent(in) :: numel, numnp, nsd, nen, ndof, nlayer
      integer, intent(in) :: conect(nen, numel)
      double precision, intent(in)  :: veloc(nen,numel,nsd,nlayer)
      double precision, intent(out) :: vm(nsd,numnp)

      integer :: i, j, l, nel, noglobal
      integer :: numocorrencianos(numnp)
      double precision :: velAcumulada(nsd,numnp)

      ! calculando a velocidade acumulada
      do nel=1,numel
         do L=1,nen
            noGlobal = conect(L,nel)
            do i=1,nlayer
               do j=1,nsd
                  velAcumulada(j,noGlobal) = velAcumulada(j,noGlobal) + veloc(L,nel,j,i)
               enddo
            end do
            numOcorrenciaNos(noGlobal)=numOcorrenciaNos(noGlobal)+1
         end do
      enddo

      ! COMPUTING AVERAGE velocity AT NODAL POINTS
      do j=1,numnp
         do i=1,nsd
             vm(i,j)=(velAcumulada(i,j)/numOcorrenciaNos(j))
         end do
      end do

   end subroutine
   !=======================================================================
   subroutine calcularVelocBordas(veloc, vborda, conect, tipoB, numel, numnp, nsd, nen)
      ! Rotina retirada das velocidades nas bordas
      ! veloc -> matriz de velocidades nos nós
      ! vborda-> matriz de velocidades normais a borda
      ! conect-> matriz de conectividades dos elementos
      ! tipoB -> matriz com os tipos de bordas
      ! numel -> quantidade de elementos - para a borda
      ! numnp -> quantidade de nós
      ! nen   -> quantidade de nós por elemento - para a borda
      ! nsd   -> dimensão do espaço (1D,2D,3D)
      IMPLICIT NONE

      integer, intent(in) :: numel, numnp, nsd, nen
      integer, intent(in) :: conect(nen, numel), tipoB(numel)
      double precision, intent(in)  :: veloc(nsd,numnp)
      double precision, intent(out) :: vborda(nen,numel)

      integer :: i, j, noBordas, tipoBorda

      ! Calculando velocidadeBordas que eh usada no calculo da vazao
      vborda=0.d0
      do i=1, numel        ! Loop nos elementos pertencentes a borda
         do j=1, nen       ! Loop nos nos dos elementos <- porque nao e relativo somente a borda
            noBordas=conect(j, i)   ! Ponteiro para no do elemento da borda
            tipoBorda = tipoB(i)    ! Ponteiro para o tipo de borda
            if(tipoBorda==1) then   ! Armazenagem do velocidade da borda na variavel 
               vborda(j,i)=-veloc(2,noBordas)
            elseif(tipoBorda==2) then 
               vborda(j,i)=veloc(1,noBordas)
            elseif(tipoBorda==3) then 
               vborda(j,i)=veloc(2,noBordas)
            elseif(tipoBorda==4) then 
               vborda(j,i)=-veloc(1,noBordas)
            elseif(tipoBorda==5) then 
               vborda(j,i)=-veloc(3,noBordas)
            elseif(tipoBorda==6) then 
               vborda(j,i)=veloc(3,noBordas)
            endif
         end do          
      end do

   end subroutine
   !=======================================================================
   subroutine calcularVelocidadeBloco(u, conecNodaisElem, numel, numnp, nsd, nen, ndof)

      use mMalha,            only: x, numelFratura
      use mGlobaisEscalares, only: ntype, numat, npint, nicode, iprtin, nrowsh, viscosidade
      use mGlobaisEscalares, only: carregamento, carregamento_ref, saltoPressao
      use mGlobaisArranjos,  only: c, mat, tiposElemBordas, campoPermeabilidade
      use mMalha,            only: local, normal, normalElemBorda, numelBordas, nenFratura, conecNodaisBordas
      use mFuncoesDeForma,   only: oneshl, oneshg, shlt, shlq, SHGQ, SHLQEN, SHLTEN, shlqen3D,  shlq3d, shg3d
      use mFuncoesDeForma,   only: shlten3D, shlt3D
      use mPotencial,        only: estrutSistEqP

      IMPLICIT NONE

      real*8, intent(inout) :: u(ndof,numnp)
      integer, intent(in) :: numel, numnp, nsd, nen, ndof
      integer, intent(in) :: conecNodaisElem(nen, numel)

      LOGICAL :: QUAD

      REAL*8 :: GRADPX, GRADPY, GRADPZ, GRADP0, velAcumulada(nsd,numnp)
      real*8 :: permeabilidadeX, permeabilidadeY, permeabilidadeZ
      real*8 :: condutividadeX, condutividadeY, condutividadeZ,  velLocal(nsd)

      INTEGER   :: I, J, L, NEL, noGlobal, m, numOcorrenciaNos(NUMNP), nelBordas, noBordas, tipoBorda

      real*8 :: xl(nsd,nen), dl(estrutSistEqP%ndof,nen)
      real*8 :: det(nen), W(npint), di(nrowsh)
      real*8, dimension(nrowsh,nen,npint) :: SHL, SHG
      real*8  :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1), vetorParaNo(nsd), normalNoLocal(nsd,nen)
 
      integer :: dof
      real*8  :: erro;

      QUAD=.TRUE.

      velAcumulada=0.d0
      numOcorrenciaNos=0

      ! GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
      w=0d0
      shl=0d0
      wCentroElem=0d0
      shlCentroElem=0d0

      if(nen==3) then                        ! Elemento(triangulo) - Espaco()
         CALL SHLTEN(SHL,NEN)
         ! call shlt(shl,w,npint,nen)
         call shlt(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==2) then        ! Elemento(quadrado) - Espaco(2D)
         CALL SHLQEN(SHL,NEN)
         ! call shlq(shl,w,npint,nen)
         call shlq(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==3) then        ! Elemento(tetraedro) - Espaco(3D)
         CALL SHLTEN3D(SHL,NEN)
         ! call shlq(shl,w,npint,nen)
         call shlt3D(shlCentroElem,wCentroElem,1,nen)
      else if(nen==8) then                   ! Elemento(cubo) - Espaco(3D)
         CALL SHLQEN3D(SHL,NEN)
         ! call shlq3d(shl,w,npint,nen)
         call shlq3d(shlCentroElem,wCentroElem,1,nen)
      end if       

      !.... BEGIN LOOP OVER THE MESH ELEMENTS
      DO NEL=1,NUMEL

         !....    LOCALIZE UNKNOWNS AND COORDINATES
         ! Armazena em xl os valores das posições (x) do elemento
         ! Armazena em DL os valores das Pressoes (estrutSistEqP%u) do elemento
         CALL LOCAL(conecNodaisElem(1,NEL),x,xl,NEN,NSD,NSD)
         CALL LOCAL(conecNodaisElem(1,NEL),estrutSistEqP%u,DL,NEN,estrutSistEqP%ndof,estrutSistEqP%ndof)       

         ! Para os elementos com fronteira nos elementos de fraturas, 
         if (numelFratura > 0) then
            call local(conecNodaisElem(:,nel), normal, normalNoLocal, nen, nsd, nsd)
            do i=1,nsd
               xCentroElem(i) = dot_product(xl(i,:), shlCentroElem(nrowsh,:,1))
            end do
         end if

         m = mat(nel)
         permeabilidadeX = campoPermeabilidade(1,nel) !c(3,m)
         permeabilidadeY = campoPermeabilidade(2,nel) !c(4,m)
         permeabilidadeZ = campoPermeabilidade(3,nel) !c(5,m)
         condutividadeX = permeabilidadeX/viscosidade
         condutividadeY = permeabilidadeY/viscosidade
         condutividadeZ = permeabilidadeZ/viscosidade

         ! 1D-2D
         !....    EVALUATE GLOBAL SHAPE FUNCTION
         if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,1,nel,1)      ! Elemento(linha) - Espaco()
         if(nen==3) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)         ! Elemento(triangulo) - Espaco()
         if(nen==4.and.nsd==2) call shgq  (xl,det,shl,shg,npint,nel,quad,nen) ! Elemento(quadrado) - Espaco(2D)
         if(nen==4.and.nsd==3) call shg3d (xl,det,shl,shg,npint,nel,nen)      ! Elemento(tetraedro) - Espaco(3D)
         if(nen==8) call shg3d (xl,det,shl,shg,npint,nel,nen)              ! Elemento(cubo) - Espaco(3D)

         !....    COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
         !....    BEGIN LOOP OVER ELEMENT NODES
         do L=1,nen
            noGlobal = conecNodaisElem(L,nel)         

            ! COMPUTE FIRST DERIVATIVES OF NODAL DISPLACEMENT
            GRADPX=0.D0
            GRADPY=0.D0
            if(nsd==3)GRADPZ=0.D0

            DO J=1,NEN
               dof = 1
               if ( (numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) then
                  if ((estrutSistEqP%lm(2,J,nel) > 0) .or. (estrutSistEqP%lm(3,J,nel) > 0)) then
                     vetorParaNo(:) = xl(:,J) - xCentroElem(:)                   
                     if (dot_product(normalNoLocal(:,J), vetorParaNo(:)) < 0d0) then
                        dof = 2
                     end if
                  end if
               end if

               GRADPX = GRADPX + shg(1,J,L)*DL(dof,J)
               GRADPY = GRADPY + shg(2,J,L)*DL(dof,J)
               if(nsd==3) GRADPZ = GRADPZ + shg(3,J,L)*DL(dof,J)
            ENDDO

            velLocal(1)=-condutividadeX * GRADPX
            velLocal(2)=-condutividadeY * GRADPY
            if(nsd==3) velLocal(3)=-condutividadeZ * GRADPZ
              
            velAcumulada(1,noGlobal) = velAcumulada(1,noGlobal) + velLocal(1)
            velAcumulada(2,noGlobal) = velAcumulada(2,noGlobal) + velLocal(2)
            if(nsd==3) velAcumulada(3,noGlobal) = velAcumulada(3,noGlobal) + velLocal(3)

            numOcorrenciaNos(noGlobal)=numOcorrenciaNos(noGlobal)+1 ! usado para cálculo da média

         end do
      ENDDO ! END OF LOOP OVER MESH ELEMENTS

      ! COMPUTING AVERAGE velocity AT NODAL POINTS
      DO J=1,NUMNP
         DO I=1,nsd
             u(I,J)=(velAcumulada(I,J)/numOcorrenciaNos(J))
         END DO
      END DO

      ! Calculando velocidadeBordas que eh usada no calculo da vazao
      velocidadeBordas=0.d0
      do i=1, numelBordas                    ! Loop nos elementos pertencentes a borda
         do j=1, nenFratura                  ! Loop nos nos dos elementos <- porque nao e relativo somente a borda
            noBordas=conecNodaisBordas(j, i) ! Ponteiro para no do elemento da borda
            tipoBorda = tiposElemBordas(i)   ! Ponteiro para o tipo de borda
            if(tipoBorda==1) then            ! Armazenagem do velocidade da borda na variavel 
               velocidadeBordas(j,i)=-u(2,noBordas)
            elseif(tipoBorda==2) then 
               velocidadeBordas(j,i)=u(1,noBordas)
            elseif(tipoBorda==3) then 
               velocidadeBordas(j,i)=u(2,noBordas)
            elseif(tipoBorda==4) then 
               velocidadeBordas(j,i)=-u(1,noBordas)
            elseif(tipoBorda==5) then 
               velocidadeBordas(j,i)=-u(3,noBordas)
            elseif(tipoBorda==6) then 
               velocidadeBordas(j,i)=u(3,noBordas)
            endif
         end do          
      end do

   end subroutine

   !=======================================================================
   !---------- ROTINAS PARA AS DESCONTINUIDADES/FALHAS --------------------
   !=======================================================================
   subroutine calcVelocElemFalha(u, velbord, conec, numnp, numelF, nsdF, nenF, ndofF, npintF, nlayer)
      ! Rotina para determinação da velocidade nos pontos de integração
      ! Armazena as velocidades na matriz "vm_falha" 
      ! - u é o campo de pressão dados nos pontos nodais
      ! - conec é a matriz de conectividades dos elementos
      ! - numel é o número de elementos para o bloco
      ! - numnp é o número de nós para o bloco (valor global)
      ! - numelF é o número de elementos para a falha
      ! - nsd é o número de dimensão de trabalho (1D,2D,3D)
      ! - nen é o número de nós para cada elemento
      ! - ndof é o número de graus de liberdade por nó
      ! - nlayer é o número de camadas

      use mGlobaisEscalares, only: nrowsh, viscosidade, saltoPressao
      use mGlobaisArranjos,  only: c, matFratura
      use mMalha,            only: local, normal, nenFratura, numelBordas, conecNodaisBordas
      use mMalha,            only: x, normalElemBorda
      use mFuncoesDeForma,   only: shlt, shlq
      use mFuncoesDeForma,   only: shl1D_intPoints, shgCodimOneSurface_points
      use mPotencial,        only: estrutSistEqP
      use mMLayer,           only: LayerMaterial

      IMPLICIT NONE

      real*8, intent(in)  :: u(ndofF,numnp)
      integer, intent(in) :: numnp, numelF, nsdF, nenF, ndofF, npintF, nlayer
      integer, intent(in) :: conec(nenF, numelF)
      real*8, intent(out) :: velbord(nenF, numelBordas, nlayer)

      LOGICAL :: QUAD
      integer :: I, J, K, L, NEL, SP, noGlobal, m, nelBordas, noBordas, tipoBorda
      integer :: dof, ncamadas
      real*8  :: salto,gradP

      real*8  :: shlCentroElem(nrowsh,nenF,1), wCentroElem(1)

      real*8, dimension(nsdF)      :: xCentroElem, vetorParaNo
      real*8, dimension(npintF)    :: det, W
      real*8, dimension(nsdF,nenF) :: xl
      real*8, dimension(ndofF,nenF) :: dl
      real*8, dimension(nsdF,nenF) :: normalNoLocal
      real*8, dimension(nsdF,nenF,npintF)   :: SHL
      real*8, dimension(nrowsh,nenF,npintF) :: SHG

      ! real*8, dimension(nlayer)      :: layKN, layDh
      real*8, dimension(:), allocatable :: layKN, layDh

      QUAD=.TRUE.

      ! GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
      shl = 0.0d0
      if(nenF==2) then        ! Elemento(linha) - Espaco(2D)
         call shl1D_intPoints(shl, w, npintF, nenF)
      else if(nenF==3) then   ! Elemento(triangulo) - Espaco(3D)
         call shlt(shl,w,npintF,nenF)
      else if(nenF==4) then   ! Elemento(quadrado) - Espaco(3D)
         call shlq(shl,w,npintF,nenF)
      end if

      !.... BEGIN LOOP OVER THE MESH ELEMENTS
      ! para armazenar as propridades na falha
      DO NEL=1,numelF

         !.... LOCALIZE UNKNOWNS AND COORDINATES
         ! Armazena em xl os valores das posições (x) do elemento
         ! Armazena em DL os valores das Pressoes (estrutSistEqP%u) do elemento
         ! Armazena em normalNoLocal os vetor normal (normal) do elemento
         CALL LOCAL(conec(1,NEL),x,XL,nenF,nsdF,nsdF)
         CALL LOCAL(conec(1,NEL),estrutSistEqP%u,DL,nenF,ndofF,ndofF)
         CALL LOCAL(conec(:,NEL), normal, normalNoLocal, nenF, nsdF, nsdF)

         ! do i=1,ndofF
         !    do j=1,nenF
         !       WRITE(*,'("dl(",i3,",",i3,")=",e12.5)') i,j,dl(i,j)
         !    enddo
         ! enddo

         !... Properties
         m = matFratura(nel)

         ncamadas = LayerMaterial(m,1)
         allocate(layDh(ncamadas))
         allocate(layKN(ncamadas))
         do k=1, ncamadas
            layDh(k) = c(2,LayerMaterial(m,k+1))
            layKN(k) = c(3,LayerMaterial(m,k+1))/viscosidade
         enddo

         ! 1D-2D
         !.... EVALUATE GLOBAL SHAPE FUNCTION
         shg  = 0.0
         call shgCodimOneSurface_points(xl, det, shl, shg, nenF, npintF, nsdF)

         !.... COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
         !.... BEGIN LOOP OVER ELEMENT NODES
         DO L=1,nenF
            noGlobal=conec(L,NEL)

            ! COMPUTE THE DERIVATIVES IN FUNCTION OF THE SPACE DIMENSION
            do k=1, ncamadas

               !.... BEGIN LOOP OVER EACH LAYER FOR THE ELEMENT
               do SP=1,nsdF
                  ! COMPUTE FIRST DERIVATIVES OF NODAL DISPLACEMENT
                  gradP=0.D0

                  ! if ( saltoPressao .EQV. .true. ) then
                  !    dof=2+k
                  ! else     
                  !    dof=1
                  ! end if

                  ! Determinação do gradiente de pressão
                  do j=1,nenF
                     if ((estrutSistEqP%lmFratura(2,J,nel)==0).or.(estrutSistEqP%lmFratura(3,J,nel)==0))then
                        dof=1
                     else
                        dof=2+k
                     endif
                     gradP=gradP+SHG(SP,J,L)*DL(dof,J)
                  enddo
                  ! print*, "grad", gradP

                  ! //TODO:: retirar if (saltoPressao)
                  ! //TODO:: armazenar velocidade tan e normal de forma correta
                  ! Determinação da velocidade na camada
                  if (saltoPressao .EQV. .true.) then
                     if (k==1)then
                        if (ncamadas.eq.1) then
                           salto=DL(2,l)-DL(1,l)
                        else
                           salto=DL(4,l)-DL(1,l)
                        endif
                     elseif (k==ncamadas) then
                       salto=DL(2,l)-DL((k+1),l)
                     else
                       salto=DL((k+3),l)-DL((k+1),l)
                     endif
                     ! print*, "salto", salto
                     vm_falha(L,nel,SP,k)=- layKN(k) * (gradP + salto*normalNoLocal(SP,l)/layDh(k))
                  else
                     vm_falha(L,nel,SP,k)=- layKN(k) * gradP
                  endif
               end do
               ! write(*,*) "salto->", salto, "normalNoLocal->", normalNoLocal(SP,l), "layDh->", layDh(k)
               ! vm_falha(nenFratura,numelFratura,nsd,ncamadas)

               ! Determinação das velocidades nas bordas
               do nelBordas=1, numelBordas
                  do J=1, nenFratura
                     noBordas=conecNodaisBordas(J, nelBordas)
                        if(noBordas==noGlobal) then
                           velbord(j,nelBordas,k)=dot_product(vm_falha(L,nel,:,k),normalElemBorda(:,nelBordas))
                        endif
                  enddo
               enddo

            enddo
         enddo
         DEALLOCATE(layDh)
         DEALLOCATE(layKN)
      enddo

   end subroutine
   !=======================================================================
   ! estrutura limitada para uma falha apenas
   subroutine calcularMediasNosFalha(veloc, vm, conect, numel, numnp, nsd, nen, ndof, nlayer)
      ! Rotina para determinação de médias de velocidade nos nós
      ! Armazena as velocidades na matriz "vm_bloco" 
      ! veloc -> matriz de velocidades por elemento nos pontos de integração
      ! vm    -> matriz de velocidades nos nós
      ! conect-> matriz de conectividades dos elementos
      ! numel -> quantidade de elementos
      ! numnp -> quantidade de nós
      ! nen   -> quantidade de nós por elemento
      ! ndof  -> número de graus de liberdade por nó
      ! nsd   -> dimensão do espaço (1D,2D,3D)
      ! nlayer-> quantidade de camadas
      use mMLayer,            only: LayerMaterial
      use mGlobaisArranjos,   only: matFratura

      IMPLICIT NONE

      integer, intent(in) :: numel, numnp, nsd, nen, ndof, nlayer
      integer, intent(in) :: conect(nen, numel)
      double precision, intent(in)  :: veloc(nen,numel,nsd,nlayer)
      double precision, intent(out) :: vm(nsd,numnp,nlayer)

      integer :: i, j, l, nel, noLocal, ncamadas
      integer :: numocorrencianos(numnp)
      double precision :: velAcumulada(nsd,numnp,nlayer)

      ! //TODO:: verificar estrutura 
      ! calculando a velocidade acumulada
      do nel=1,numel
         ncamadas = LayerMaterial(matFratura(nel),1)
         do L=1,nen
            noLocal  = nel+L-1
            do i=1,ncamadas
               do j=1,nsd
                  velAcumulada(j,noLocal,i) = velAcumulada(j,noLocal,i) + veloc(L,nel,j,i)
               enddo
            end do
            numOcorrenciaNos(noLocal)=numOcorrenciaNos(noLocal)+1
         end do
      enddo

      ! COMPUTING AVERAGE velocity AT NODAL POINTS
      do i=1,numnp
         do j=1,nsd
            do l=1,nlayer
               vm(j,i,l)=(velAcumulada(j,i,l)/numOcorrenciaNos(i))
            end do
         end do
      end do
   end subroutine

   !=======================================================================
   subroutine calcularVelocidadeFalha(u, conecNodaisElem, numel, numnp, nsd, nenFratura, ndof)
      use mMalha,            only: x
      use mGlobaisEscalares, only: ntype, numat, npint, npintFratura, nicode, iprtin, nrowsh, zero, viscosidade
      use mGlobaisEscalares, only: carregamento, carregamento_ref, saltoPressao
      use mGlobaisArranjos,  only: c, matFratura, tiposElemBordas
      use mMalha,            only: local, normal, numelBordas, conecNodaisBordas, normalElemBorda
      use mFuncoesDeForma,   only: oneshg, oneshl, shl1D_intPoints, shgCodimOneSurface_points
      use mFuncoesDeForma,   only: shlq, shlt,SHGQ, SHLQEN, SHLTEN, shl1D_elemNodes, shg1D_points
      use mPotencial,        only: estrutSistEqP
      use mMLayer,           only: LayerMaterial

      IMPLICIT NONE

      real*8,  intent(inout) :: u(nsd,numnp)
      integer, intent(in) :: numel, numnp, nsd, nenFratura, ndof
      integer, intent(in) :: conecNodaisElem(nenFratura, numel)

      LOGICAL QUAD

      REAL*8 :: GRADP_X, GRADP_Y, GRADP_Z, velAcumulada(nsd,numnp)
      real*8 :: porosidade_ref, permeabilidade_ref, coef_poisson, mod_young, rigidezNormalInicial
      real*8 :: beta, tensao_total_media, tensao_total_media_ref, div_u, porosidade
      real*8 :: abertura, abertura_ref
      real*8 :: permeabilidade, condutividade

      INTEGER :: I, J, L, NEL, noGlobal, m, numOcorrenciaNos(NUMNP), nelBordas, noBordas

      real*8 :: xl(nsd,nenFratura), dl(estrutSistEqP%ndof,nenFratura)
      real*8 :: det(npint), detp(npint), W(npint)
      real*8, dimension(nsd,nenFratura,nenFratura)    :: SHL
      real*8, dimension(nrowsh,nenFratura,nenFratura) :: SHG

      real*8  :: diff(nsd), length, sinTheta, cosTheta
      real*8  :: normalNoLocal(nsd,nenFratura), deslocamentoCritico
      real*8  :: velLocal(nsd), salto
      integer :: dof

      real*8, dimension(:), allocatable :: layKN, layDh, gradP_FX, gradP_FY, gradP_FZ
      integer*4:: ncamadas, k, posI1, posI2, posI3, posJ1, posJ2, posJ3

      QUAD=.TRUE.

      velAcumulada=0.d0
      numOcorrenciaNos=0

      ! GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
      shl = 0.0d0
      if(nenFratura==2) then           ! Elemento(linha) - Espaco(2D)
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then      ! Elemento(triangulo) - Espaco(3D)
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then      ! Elemento(quadrado) - Espaco(3D)
         call shlq(shl,w,npintFratura,nenFratura)
      end if

      !.... BEGIN LOOP OVER THE MESH ELEMENTS
      ! para armazenar as propridades na falha
      ! ncamadas = estrutSistEqP%ndof-2
      ! allocate(layDh(ncamadas))
      ! allocate(layKN(ncamadas))
      ! allocate(gradP_FX(ncamadas))
      ! allocate(gradP_FY(ncamadas))
      ! allocate(gradP_FZ(ncamadas))
      DO NEL=1,NUMEL

         !.... LOCALIZE UNKNOWNS AND COORDINATES
         ! Armazena em xl os valores das posições (x) do elemento
         ! Armazena em DL os valores das Pressoes (estrutSistEqP%u) do elemento
         ! Armazena em normalNoLocal os vetor normal (normal) do elemento
         CALL LOCAL(conecNodaisElem(1,NEL),x,XL,nenFratura,NSD,NSD)
         CALL LOCAL(conecNodaisElem(1,NEL),estrutSistEqP%u,DL,nenFratura,estrutSistEqP%ndof,estrutSistEqP%ndof)
         CALL LOCAL(conecNodaisElem(:,nel), normal, normalNoLocal, nenFratura, nsd, nsd)

         !....... form stiffness matrix
         !... Properties
         m = matFratura(nel)
         ncamadas = LayerMaterial(m,1)
         allocate(layDh(ncamadas))
         allocate(layKN(ncamadas))
         allocate(gradP_FX(ncamadas))
         allocate(gradP_FY(ncamadas))
         allocate(gradP_FZ(ncamadas))
         do k=1, ncamadas
            layDh(k) = c(2,LayerMaterial(m,k+1))
            layKN(k) = c(3,LayerMaterial(m,k+1))/viscosidade
         enddo

         ! Propriedades da camada 1 e 3
         ! m = matFratura(nel)
         ! do k=1, ncamadas
         !    layDh(k) = c(2,m)
         !    layKN(k) = c(3,m)/viscosidade
         ! enddo

         ! 1D-2D
         !.... EVALUATE GLOBAL SHAPE FUNCTION
         shg  = 0.0
         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)

         !.... COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
         !.... BEGIN LOOP OVER ELEMENT NODES
         DO L=1,nenFratura
            noGlobal=conecNodaisElem(L,NEL)
            ! condutividade = permeabilidade/viscosidade ! sera trocado por layKN

            !.... BEGIN LOOP OVER EACH LAYER FOR THE ELEMENT
            do k=1, ncamadas
               ! COMPUTE FIRST DERIVATIVES OF NODAL DISPLACEMENT
               gradP_FX(k)=0.D0
               gradP_FY(k)=0.D0
               gradP_FZ(k)=0.D0

               if ( saltoPressao .EQV. .true. ) then
                  dof=2+k
               else     
                  dof=1
               end if

               DO J=1,nenFratura
                  if ((estrutSistEqP%lmFratura(2,J,nel)==0).or.(estrutSistEqP%lmFratura(3,J,nel)==0))then
                     dof=1
                  else
                     dof=2+k
                  endif
          
                  gradP_FX(k) = gradP_FX(k) + SHG(1,J,L)*DL(dof,J)
                  gradP_FY(k) = gradP_FY(k) + SHG(2,J,L)*DL(dof,J)
                  if(nsd==3) gradP_FZ(k) = gradP_FZ(k) + SHG(3,J,L)*DL(dof,J)
               ENDDO
            
               if (saltoPressao .EQV. .true.) then
                  if (k==1)then
                    salto=DL(4,l)-DL(1,l)
                  elseif (k==ncamadas) then
                    salto=DL(2,l)-DL((k+1),l)
                  else
                    salto=DL((k+3),l)-DL((k+1),l)
                  endif
               
                  velLocal(1)=- layKN(k) * (gradP_FX(k) + salto*normalNoLocal(1,l)/layDh(k))
                  velLocal(2)=- layKN(k) * (gradP_FY(k) + salto*normalNoLocal(2,l)/layDh(k))
                  if(nsd==3) velLocal(3)=- layKN(k) * (gradP_FZ(k) + salto*normalNoLocal(3,l)/layDh(k))
               else
                  velLocal(1)=- layKN(k) * gradP_FX(k) 
                  velLocal(2)=- layKN(k) * gradP_FY(k)
                  if(nsd==3) velLocal(3)=- layKN(k) * gradP_FZ(k)
               endif
               
               
               do nelBordas=1, numelBordas
                  do J=1, nenFratura
                     noBordas=conecNodaisBordas(J, nelBordas)
                        if(noBordas==noGlobal) then
                           velocidadeBordasFalhas(j,nelBordas,k)=dot_product(velLocal(:),normalElemBorda(:,nelBordas))
                           print*,'velo(',j,',',nelBordas,',',k,')=',velocidadeBordasFalhas(j,nelBordas,k)
                        endif
                  enddo
               enddo

               velAcumulada(1,noGlobal)= velAcumulada(1,noGlobal) + velLocal(1)
               velAcumulada(2,noGlobal)= velAcumulada(2,noGlobal) + velLocal(2)
               if(nsd==3)velAcumulada(3,noGlobal)= velAcumulada(3,noGlobal) + velLocal(3)

               numOcorrenciaNos(noGlobal)=numOcorrenciaNos(noGlobal)+1 ! usado para cálculo da média

            enddo ! Para as camadas
         ENDDO ! END OF LOOP OVER ELEMENT NODES
         DEALLOCATE(layDh)
         DEALLOCATE(layKN)
         DEALLOCATE(gradP_FX)
         DEALLOCATE(gradP_FY)
         DEALLOCATE(gradP_FZ)
      ENDDO ! END OF LOOP OVER MESH ELEMENTS

      ! deallocate(layDh)
      ! deallocate(layKN)
      ! deallocate(gradP_FX)
      ! deallocate(gradP_FY)
      ! deallocate(gradP_FZ)

      ! COMPUTING AVERAGE STRESSES AT NODAL POINTS
      DO NEL=1,NUMEL
         DO L=1,nenFratura
             noGlobal=conecNodaisElem(L,NEL)
             if(numOcorrenciaNos(noGlobal).ne.0) then
                u(:,noGlobal)=(velAcumulada(:,noGlobal)/numOcorrenciaNos(noGlobal))
                numOcorrenciaNos(noGlobal)=0
             endif
         END DO
      END DO

   end subroutine

   !=======================================================================
   !---------- ROTINAS PARA AS DESCONTINUIDADES/FRATURAS ------------------
   !=======================================================================
   subroutine calcularVelocidadeFratura(u, conecNodaisElem, numel, numnp, nsd, nenFratura, ndof)
      use mMalha,            only: x
      use mGlobaisEscalares, only: ntype, numat, npint, npintFratura, nicode, iprtin, nrowsh, zero, viscosidade
      use mGlobaisEscalares, only: carregamento, carregamento_ref, saltoPressao
      use mGlobaisArranjos,  only: c, matFratura, tiposElemBordas
      use mMalha,            only: local, normal, numelBordas, conecNodaisBordas, normalElemBorda
      use mFuncoesDeForma,   only: oneshg, oneshl, shl1D_intPoints, shgCodimOneSurface_points
      use mFuncoesDeForma,   only: shlq, shlt,SHGQ, SHLQEN, SHLTEN, shl1D_elemNodes, shg1D_points
      use mPotencial,        only: estrutSistEqP

      IMPLICIT NONE

      real*8,  intent(inout) :: u(nsd,numnp)
      integer, intent(in) :: numel, numnp, nsd, nenFratura, ndof
      integer, intent(in) :: conecNodaisElem(nenFratura, numel)

      LOGICAL QUAD

      REAL*8 :: GRADP_X, GRADP_Y, GRADP_Z, velAcumulada(nsd,numnp)
      real*8 :: porosidade_ref, permeabilidade_ref, coef_poisson, mod_young, rigidezNormalInicial
      real*8 :: beta, tensao_total_media, tensao_total_media_ref, div_u, porosidade
      real*8 :: abertura, abertura_ref
      real*8 :: permeabilidade, condutividade

      INTEGER :: I, J, L, NEL, noGlobal, m, numOcorrenciaNos(NUMNP), nelBordas, noBordas

      real*8 :: xl(nsd,nenFratura), dl(estrutSistEqP%ndof,nenFratura)
      real*8 :: det(npint), detp(npint), W(npint)
      real*8, dimension(nsd,nenFratura,nenFratura)    :: SHL
      real*8, dimension(nrowsh,nenFratura,nenFratura) :: SHG

      real*8  :: diff(nsd), length, sinTheta, cosTheta
      real*8  :: normalNoLocal(nsd,nenFratura), deslocamentoCritico
      real*8  :: velLocal(nsd), salto
      integer :: dof

      QUAD=.TRUE.

      velAcumulada=0.d0
      numOcorrenciaNos=0

      ! GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
      shl = 0.0d0
      if(nenFratura==2) then
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then
         call shlq(shl,w,npintFratura,nenFratura)
      end if

      !.... BEGIN LOOP OVER THE MESH ELEMENTS
      DO NEL=1,NUMEL

         !.... LOCALIZE UNKNOWNS AND COORDINATES
         CALL LOCAL(conecNodaisElem(1,NEL),x,XL,nenFratura,NSD,NSD)
         CALL LOCAL(conecNodaisElem(1,NEL),estrutSistEqP%u,DL,nenFratura,estrutSistEqP%ndof,estrutSistEqP%ndof)
         CALL LOCAL(conecNodaisElem(:,nel), normal, normalNoLocal, nenFratura, nsd, nsd)

         !....... form stiffness matrix
         !... Properties
         m = matFratura(nel)
         abertura = c(2,m)
         permeabilidade = c(3,m)

         ! 1D-2D
         !.... EVALUATE GLOBAL SHAPE FUNCTION
         shg  = 0.0
         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)

         !.... COMPUTE FIRST DERIVATIVES OF NODAL PRESSURE
         !.... BEGIN LOOP OVER ELEMENT NODES
         DO L=1,nenFratura
            noGlobal=conecNodaisElem(L,NEL)
            condutividade = permeabilidade/viscosidade

            ! COMPUTE FIRST DERIVATIVES OF NODAL DISPLACEMENT
            GRADP_X=0.D0
            GRADP_Y=0.D0      
            GRADP_Z=0.D0  
            
            if ( saltoPressao .EQV. .true. ) then
               dof=3
            else     
               dof=1
            end if

            DO J=1,nenFratura
               if ((estrutSistEqP%lmFratura(2,J,nel)==0).or.(estrutSistEqP%lmFratura(3,J,nel)==0))then
                  dof=1
               else
                  dof=3
               endif
          
               GRADP_X = GRADP_X + SHG(1,J,L)*DL(dof,J)
               GRADP_Y = GRADP_Y + SHG(2,J,L)*DL(dof,J)
               if(nsd==3) GRADP_Z = GRADP_Z + SHG(3,J,L)*DL(dof,J)
            ENDDO
            
            if (saltoPressao .EQV. .true.) then 
               salto=DL(2,l)-DL(1,l)
               
               velLocal(1)=- condutividade * (GRADP_X + salto*normalNoLocal(1,l)/abertura)
               velLocal(2)=- condutividade * (GRADP_Y + salto*normalNoLocal(2,l)/abertura)
               if(nsd==3) velLocal(3)=- condutividade * (GRADP_Z + salto*normalNoLocal(3,l)/abertura)
            else
               velLocal(1)=- condutividade * GRADP_X 
               velLocal(2)=- condutividade * GRADP_Y
               if(nsd==3) velLocal(3)=- condutividade * GRADP_Z
            endif
               
            do nelBordas=1, numelBordas
               do J=1, nenFratura
                  noBordas=conecNodaisBordas(J, nelBordas)
                  if(noBordas==noGlobal) then
                     velocidadeBordasFratura(j,nelBordas)=dot_product(velLocal(:),normalElemBorda(:,nelBordas))
                  endif
               enddo
            enddo

            velAcumulada(1,noGlobal)= velAcumulada(1,noGlobal) + velLocal(1)
            velAcumulada(2,noGlobal)= velAcumulada(2,noGlobal) + velLocal(2)
            if(nsd==3)velAcumulada(3,noGlobal)= velAcumulada(3,noGlobal) + velLocal(3)

            numOcorrenciaNos(noGlobal)=numOcorrenciaNos(noGlobal)+1 ! usado para cálculo da média

         ENDDO ! END OF LOOP OVER ELEMENT NODES
      ENDDO ! END OF LOOP OVER MESH ELEMENTS

      ! COMPUTING AVERAGE STRESSES AT NODAL POINTS
      DO NEL=1,NUMEL
         DO L=1,nenFratura
             noGlobal=conecNodaisElem(L,NEL)
             if(numOcorrenciaNos(noGlobal).ne.0) then
                u(:,noGlobal)=(velAcumulada(:,noGlobal)/numOcorrenciaNos(noGlobal))
                numOcorrenciaNos(noGlobal)=0
             endif
         END DO
      END DO

   end subroutine calcularVelocidadeFratura

   !=======================================================================
   !-------------- ROTINAS PARA AS INTERSECCOES ---------------------------
   !=======================================================================
   ! subroutine calcularInterseccoesBordaDescont(flagGerarConectividades)
   !    use mGlobaisArranjos,  only: tiposElemBordas
   !    use mMalha,            only: nenFratura, numelBordas, conecNodaisFratura, numelFraturasNaBorda
   !    use mMalha,            only: conectividadesInterseccao, conecNodaisBordas, numelFratura, conecNodaisFraturasNaBorda
   !    use mGlobaisArranjos,  only: matFratura, matFraturasNaBorda,  tipoBordaElemFratNaBorda
   !    use mLeituraescrita,   only: iconectIntersecao
         
   !    logical :: flagGerarConectividades
          
   !    integer:: nelBordas,  i, n, nel
   !    integer:: nosArestaElemBorda(2), nosArestaElemFratura(2)
   !    logical:: interceptou
   !    integer :: arestaElemBorda, arestaElemFratura, nelFratura, tipoBorda
   !    logical, allocatable:: elemFraturaInterceptou(:)   
   !    integer :: iconectNodaisFratNaBorda, iTipoBordaElemFratNaBorda

      
   !    flagGerarConectividades = .false.
      
   !    !! Para este caso 3D, com tetraedros, quero que a primeira coluna seja 3; 
   !    !! Não sei se nenFratura é o mais correto.
   !    allocate(conectividadesInterseccao(nenFratura,2,numelBordas))
   !    conectividadesInterseccao=0
           
   !    open(unit=iconectIntersecao,  file='conectividadesIntersecao.dat', status='old', err=100) 
   !    print*, "Encontrado arquivo conectividadesIntersecao.dat"
   !    do nelBordas=1, numelBordas
   !       do i=1,nenFratura
   !          read(iconectIntersecao, '(1i10,2i10)') n,conectividadesInterseccao(i,:,nelBordas)
   !          !  write(*,'(2i10,2i10)') nelBordas, n, conectividadesInterseccao(i,:,nelBordas)
   !       end do   
   !    end do 
   !    close(iconectIntersecao)
            
   !    open(unit=iconectNodaisFratNaBorda,  file='conectividadesFraturasNasBordas.inc', status='old', err=100)
   !    print*, "arquivo conectividadesFraturasNasBordas.inc encontrado"
   !    ! gerar as conectividades através da rotina leituraGeracaoConectividadesDS()
   !    flagGerarConectividades = .true.
   !    close(iconectNodaisFratNaBorda)
      
   ! end subroutine

   !=======================================================================
   subroutine calcularInterseccoesDesconBorda(flag)
      use mGlobaisArranjos,  only: tiposElemBordas, aberturasBorda
      use mGlobaisArranjos,  only: c, matFratura, matDescontNaBorda, tipoBordaElemDescontNaBorda
      use mMalha,            only: local, x, nenFratura, nsd, numelBordas, conecNodaisFratura
      use mMalha,            only: conectividadesInterseccaoFLuxCons, conecNodaisBordas, numelFratura, normal
      use mMalha,            only: numelDescontNaBorda, conecNodaisDesconNaBorda
      use mLeituraescrita,   only: iconectIntersecao

      logical :: flag      ! flag para indicar a criação dos arquivos de conectividades das bordas
      logical:: interceptou
      logical, allocatable:: elemFraturaInterceptou(:)

      integer:: nelBordas,  i, j, n, nel
      integer:: nosArestaElemBorda(2), nosArestaElemFratura(2)
      integer :: arestaElemBorda, arestaElemFratura, nelFratura, tipoBorda
      integer :: iconectNodaisFratNaBorda, itipoBordaElemDescontNaBorda

      ! Inicialização de variáveis
      flag=.false.      
      allocate(conectividadesInterseccaoFLuxCons(nenFratura,2,numelBordas))
      conectividadesInterseccaoFLuxCons=0

      ! Leitura do arquivo com as conectividades na interseção
      open(unit=iconectIntersecao,  file='conectividadesIntersecao.dat', status='old', err=100) 
      print*, "Encontrado arquivo conectividadesIntersecao.dat"
      do nelBordas=1, numelBordas
         ! TODO: corrigir variavel (conectividadesInterseccao)
         ! read(iconectIntersecao, '(1i10,2i10)') n,conectividadesInterseccao(:,nelBordas)
         do i=1,nenFratura
             read(iconectIntersecao, '(1i10,2i10)') n,conectividadesInterseccaoFLuxCons(i,:,nelBordas) ! original
          end do   
      end do 
      close(iconectIntersecao)

      ! Leitura do arquivo com as conectividades na interseção
      open(unit=iconectNodaisFratNaBorda,  file='conectividadesFraturasNasBordas.inc', status='old', err=100)
      print*, "arquivo conectividadesFraturasNasBordas.inc encontrado"
      flag = .true. ! gerar as conectividades através da rotina leituraGeracaoConectividadesDS()
      close(iconectNodaisFratNaBorda)
      
      return

      100  print*, "Gerar condições de fraturas na bordo"
      allocate(elemFraturaInterceptou(numelFratura)); elemFraturaInterceptou=.false.
      numelDescontNaBorda=0
      tipoBorda=0      
      
      ! por ser a primeira vez, não sei qual é o tamanho necessário  
      allocate(conecNodaisDesconNaBorda(nenFratura-1, numelFratura))
      allocate(matDescontNaBorda(numelFratura))
      allocate(tipoBordaElemDescontNaBorda(numelFratura))
     
      do nelBordas=1, numelBordas
         if ( tipoBorda .NE. tiposElemBordas(nelBordas) ) elemFraturaInterceptou=.false.
         tipoBorda = tiposElemBordas(nelBordas)
         
         
         do arestaElemBorda=1, nenFratura
            interceptou=.false.
            nosArestaElemBorda(1)= conecNodaisBordas(arestaElemBorda,nelBordas) 
            if(arestaElemBorda==nenFratura) then            
               nosArestaElemBorda(2)= conecNodaisBordas(1,nelBordas) 
            else
               nosArestaElemBorda(2)= conecNodaisBordas(arestaElemBorda+1,nelBordas) 
            endif
         
            do nelFratura=1, numelFratura
               do arestaElemFratura=1, nenFratura
                  nosArestaElemFratura(1)= conecNodaisFratura(arestaElemFratura,nelFratura)
                  if(arestaElemFratura==nenFratura) then            
                     nosArestaElemFratura(2)= conecNodaisFratura(1,nelFratura)
                  else
                     nosArestaElemFratura(2)= conecNodaisFratura(arestaElemFratura+1,nelFratura)
                  endif
                  
                  if((nosArestaElemBorda(1)==nosArestaElemFratura(1)).or.(nosArestaElemBorda(1)==nosArestaElemFratura(2))) then
                     if((nosArestaElemBorda(2)==nosArestaElemFratura(1)).or.(nosArestaElemBorda(2)==nosArestaElemFratura(2))) then
                        ! Corrigir variavel conectividadesInterseccao
                        ! conectividadesInterseccaoFLuxCons(arestaElemBorda,:,nelBordas)=nosArestaElemBorda(:)
                        interceptou=.true.
                        if ( elemFraturaInterceptou(nelFratura) .EQV. .false. ) then          
                           elemFraturaInterceptou(nelFratura)=.true.
                           numelDescontNaBorda = numelDescontNaBorda + 1
                           conecNodaisDesconNaBorda(1, numelDescontNaBorda) = nosArestaElemFratura(1) 
                           conecNodaisDesconNaBorda(2, numelDescontNaBorda) = nosArestaElemFratura(2) 
                           matDescontNaBorda(numelDescontNaBorda) = matFratura(nelFratura)
                           tipoBordaElemDescontNaBorda(numelDescontNaBorda) = tipoBorda
!                         write(1000, *)  nelFratInterceptou, mat, nosArestaElemFratura(1), nosArestaElemFratura(2)                           
!                         print*,"aqui", nelBordas, nosArestaElemBorda, nelFratura, nosArestaElemFratura
                        end if ! elemFraturaInterceptou   
                     endif
                  endif
               
                  if(interceptou) exit
               end do !arestaElemFratura
               
               if(interceptou) exit               
            enddo !nelFratura
            
            if(interceptou) exit         
         end do !arestaElemBorda
      
      end do !nelBordas

      open(unit=iconectIntersecao,  file='conectividadesIntersecao.dat', status='new', err=200)
      do nelBordas=1, numelBordas
         do i=1,nenFratura
            ! Corrigir variavel conectividadesInterseccao
            !  write(iconectIntersecao,'(1i10,2i10)') nelBordas, conectividadesInterseccaoFLuxCons(i,:,nelBordas)
          end do   
      end do        
      close(iconectIntersecao)
             
      iconectNodaisFratNaBorda = 50
      open(unit=iconectNodaisFratNaBorda,  file='conectividadesFraturasNasBordas.inc', status='new', err=300)
      write (iconectNodaisFratNaBorda, "('*numelDescontNaBorda{')") 
      write (iconectNodaisFratNaBorda, "(I10)") numelDescontNaBorda
      write (iconectNodaisFratNaBorda, "('}')") 
      write (iconectNodaisFratNaBorda, "('*conectividades_nodais_fraturas_borda{')")    
      do nel=1, numelDescontNaBorda
         write (iconectNodaisFratNaBorda, "(4I10)") nel, matDescontNaBorda(nel), &
            conecNodaisDesconNaBorda(1, nel), conecNodaisDesconNaBorda(2, nel)
      end do   
      write (iconectNodaisFratNaBorda, "('         0')") 
      write (iconectNodaisFratNaBorda, "('}')") 
      write (iconectNodaisFratNaBorda, "(' ')")        
      close(iconectNodaisFratNaBorda)                                      
      
      
      itipoBordaElemDescontNaBorda = 60
      open(unit=itipoBordaElemDescontNaBorda,  file='tipoBordaElemDescontNaBorda.dat', status='new', err=400)
      do nel=1, numelDescontNaBorda
         write (itipoBordaElemDescontNaBorda, "(2I10)") nel, tipoBordaElemDescontNaBorda(nel)
      end do         
      close(itipoBordaElemDescontNaBorda)  
      return
      
      200 write(*,*) ' arquivo: conectividadesIntersecao.dat NAO pode ser criado'
      stop 
      
      300 write(*,*) ' arquivo: conectividadesFraturasNasBordas.inc NAO pode ser criado'
      stop
      
      400 write(*,*) ' arquivo: tipoBordaElemDescontNaBorda.dat NAO pode ser criado'       
      stop 
      

      deallocate(elemFraturaInterceptou)

      
   end subroutine calcularInterseccoesDesconBorda
   !=======================================================================
   subroutine calcularInterseccoesFraturaBorda0(estrutSistEqP)
      use mGlobaisEscalares, only: npintFratura, nrowsh
      use mGlobaisArranjos,  only: tiposElemBordas, aberturasBorda
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shgq, shlq, shgCodimOneSurface_points, shl1D_intPoints
      use mMalha,            only: local, x, nenFratura, nsd, numelBordas, conecNodaisFratura
      use mMalha,            only: conectividadesInterseccao, conecNodaisBordas, numelFratura, normal
      use mGlobaisArranjos,  only: c, matFratura
         
      type (estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqP
          
      real*8 :: xl(nsd,nenFratura)
      real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)
      real*8 :: det(npintFratura), w(npintFratura)
      integer:: nelBordas,  i, j
      integer:: noFratura, noDeInterseccao(2), cont
      logical:: elementoDeFratura

      w=0.0
      shl=0.0
      
      allocate(conectividadesInterseccao(2,numelBordas)); conectividadesInterseccao=0
     
      do nelBordas=1, numelBordas

         elementoDeFratura=.false.
         noDeInterseccao=0
         cont=1

         do i=1, numelFratura
            do j=1, nenFratura
               noFratura=conecNodaisFratura(j,i)
               ! print*, nelBordas, noFratura
               if((noFratura==conecNodaisBordas(1,nelBordas)).or.(noFratura==conecNodaisBordas(2,nelBordas)).or.&
                        (noFratura==conecNodaisBordas(3,nelBordas)).or.(noFratura==conecNodaisBordas(4,nelBordas))) then
                  elementoDeFratura=.true.                 
                  ! print*, nelBordas
                  
                  if(noFratura==noDeInterseccao(1).or.noFratura==noDeInterseccao(2))then
                     ! não faz nada
                  else
                     noDeInterseccao(cont)=noFratura
                     cont=cont+1 
                     ! print*, "cont" cont
                  endif
               endif
            enddo
         enddo

         if(elementoDeFratura.eqv..true.) then
            ! print*, nelBordas, noDeInterseccao
            conectividadesInterseccao(:,nelBordas)=noDeInterseccao(:)
         endif

      end do
      
   end subroutine calcularInterseccoesFraturaBorda0

   !=======================================================================
   !---------------- ROTINAS PARA CALCULO DA VAZAO ------------------------
   !=======================================================================
   subroutine calcularVazao(estrutSistEqP, Lx, Ly, Lz)
      use mGlobaisEscalares, only: npintFratura, nrowsh, tempo, deltaT, pi, viscosidade
      use mGlobaisArranjos,  only: tiposElemBordas, aberturasBorda, campoPermeabilidade, vazao
      use mGlobaisArranjos,  only: c, matFratura, tiposElemBordas
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shgq, shlq, shgCodimOneSurface_points
      use mfuncoesDeForma,   only: shl1D_intPoints, shlqen, shlten
      use mMalha,            only: local, x, nenFratura, nsd, numelBordas, conecNodaisBordas
      use mMalha,            only: numelFratura, normal, conectividadesInterseccao
      use mMLayer,           only: LayerMaterial

      type (estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqP
          
      real*8 :: xl(nsd,nenFratura)
      real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura), shlParaArea(nsd,nenFratura,npintFratura)
      real*8 :: det(npintFratura), w(npintFratura)
      real*8  :: wDet, di(nrowsh)
      integer*4:: nee, nelBordas, tipoBorda, l, i, indexI, sd_i, ndofDesloc      
      real*8  :: xCentroElem(nsd), shlCentroElem(nsd,nenFratura,1), wCentroElem(1), vetorParaNo(nsd), normalNoLocal(nsd,nenFratura)
      real*8  :: area(6), vn_ptoInt
      real*8  :: pressaoDir, pressaoEsq, pressaoInf, pressaoSup, pressaoFrente, pressaoAtras
      real*8  :: diffP, Lx, Ly, Lz, distanciaEntrePontos(3), distancia, areaFratura
      real*8  :: tempoAdimensionalAtual, tempoAdimensionalAnterior, vazaoAdimensionalAtual, vazaoAdimensionalAnterior
      real*8 :: vazaoAnterior(14), mltil, Kestimado, Ppoco, Pinicial, permMatriz, beta

      real*8, dimension(:), allocatable :: layDh
      integer*4:: ncamadas, k, m
      real*8 :: abertura

      w=0.0
      shl=0.0
      wCentroElem=0d0
      vazaoAnterior(:)=vazao(:)

      ! Determinação das funções de forma.
      if(nenFratura==2) then
         ! Elemento(Linha) - Espaco(2D)
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then
         ! Elemento(triangulo) - Espaco(3D)
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then
         ! Elemento(quadrado) - Espaco(3D)
         call shlq(shl,w,npintFratura,nenFratura)
         call shlqen(shl,nenFratura)
      end if
      
      ! Loop sobre os elementos das bordas. 
      area =0.d0
      vazao=0.d0
      do nelBordas=1,numelBordas

         ! Armazena em xl os valores das posições (x) do elemento         
         call local(conecNodaisBordas(:,nelBordas), x, xl, nenFratura, nsd, nsd)         
         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)

         !... Identificacao do tipo de borda
         ! 1: Inferior
         ! 2: Direita
         ! 3: Superior
         ! 4: Esquerda
         tipoBorda = tiposElemBordas(nelBordas)
         if(tipoBorda==1) then
            pressaoFrente=estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==2) then
            pressaoDir=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==3) then
            pressaoAtras= estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==4) then
            pressaoEsq=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==5) then
            pressaoInf=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         else 
            pressaoSup=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         endif
         area(tipoBorda) = area(tipoBorda) + dot_product(w(:),det(:))
         ! if(tipoBorda==2) then
         !    print*, "det=", det
         !    print*, "w=", w
         !    print*, "dot=", dot_product(w(:),det(:))
         ! endif

         !.... loop on integration points         
         ! Computo da vazao para o macico rochoso
         do l=1,npintFratura
            !velocidade normal avaliada no ponto de integracao! !!!ou será no nó?
            vn_ptoInt = dot_product(velocidadeBordas(:,nelBordas), shg(nrowsh,:,l))
            vazao(tipoBorda) = vazao(tipoBorda) + vn_ptoInt*w(l)*det(l)
         end do

         ! Computo da vazao para a descontinuidade
         if(numelFratura>0) then
            ! Espaco(2D)
            if(nsd==2) then
               ! Loop nos nós do elemento
               ! do i=1, nenfratura
               !    ! Testa se o nó do elemento de bordo contém descontinuidade
               !    if()then
               !       ! Determinação das propriedades das camadas.
               !       ncamadas = estrutSistEqP%ndof-2
               !       ! ncamadas = LayerMaterial(matFratura(nel),1)
               !       do k=1, ncamadas
               !          m = matFratura(nelBordas)
               !          abertura=c(2,LayerMaterial(m,k+1))
               !          area(tipoborda)=area(tipoborda)+abertura/2.d0
               !          vazao(tipoborda)=vazao(tipoborda)+velocidadeBordasFalhas(i,nelbordas,k)*abertura/2.d0
               !       end do
               !    end if
               ! end do

               ! OLD
               do i=1, nenfratura
                  ! write(*,*) 'abert->',aberturasborda(i,nelbordas)
                  ! print*,'aberturasborda',aberturasborda(i,nelbordas)
                  ! write(*,*) 'Veloc->',velocidadeBordasFalhas(i,nelbordas,1)
                  ! write(*,*) 'Veloc->',velocidadeBordasFalhas(i,nelbordas,3)
                  if(aberturasborda(i,nelbordas)>1.d-10) then
                     area(tipoborda)=area(tipoborda)+aberturasborda(i,nelbordas)/2.d0
                     vazao(tipoborda)=vazao(tipoborda)+velocidadeBordasFalhas(i,nelbordas,1)*aberturasborda(i,nelbordas)/2.d0
                     ! print*,'VAZAO->',vazao(tipoborda)
                  end if
               end do
            ! Espaco(3D)
            else                 
               if((conectividadesInterseccao(1,nelBordas).ne.0)) then   ! Caso com intersecao de descontinuidade
                  distanciaEntrePontos(:)=(x(:,conectividadesInterseccao(1,nelBordas)) - x(:,conectividadesInterseccao(2,nelBordas)))
                  distancia=sqrt(dot_product(distanciaEntrePontos(:),distanciaEntrePontos(:)))
                  
                  do i=1, nenFratura
                     if(aberturasBorda(i,nelBordas)>1.d-10) exit
                  enddo
                  areaFratura=(aberturasBorda(i,nelBordas)/2.d0)*distancia
                  area(tipoBorda)=area(tipoBorda)+areaFratura
                  vazao(tipoBorda)=vazao(tipoBorda)+velocidadeBordasFratura(i,nelBordas)*areaFratura
               endif
            endif
         end if
            
      end do

      !TODO: Acho que a vazão esta calculada errada, deveria ser feito a regularização da forma:
      ! vazao(1) = vazao(1)*(area(2)/area(1))
      ! vazao(2) = vazao(2)*(area(1)/area(2))
      ! vazao(3) = vazao(3)*(area(4)/area(3))
      ! vazao(4) = vazao(4)*(area(3)/area(4))

      ! write(657,'(5es20.8)') tempo, abs(vazao(1:4))
      write(*,'(a, 6e15.5)') "Areas", area
      ! write(*,'(a, 6e15.5)') "Pressões", pressaoFrente, pressaoDir, pressaoAtras, pressaoEsq, pressaoInf, pressaoSup
      
      write(*,'(a,1es20.8)') "Vazao na borda da frente   = ", abs(vazao(1))
      write(*,'(a,1es20.8)') "Vazao na borda da direita  = ", abs(vazao(2))
      write(*,'(a,1es20.8)') "Vazao na borda de trás     = ", abs(vazao(3))
      write(*,'(a,1es20.8)') "Vazao na borda da esquerda = ", abs(vazao(4))
      if(nsd==3) then
         write(*,'(a,1es20.8)') "Vazao na borda inferior    = ",abs(vazao(5))
         write(*,'(a,1es20.8)') "Vazao na borda superior    = ",abs(vazao(6))
      endif
   
      Lx=(maxval(x(1,:)) - minval(x(1,:)))
      Ly=(maxval(x(2,:)) - minval(x(2,:)))
      if(nsd==3) Lz=(maxval(x(3,:)) - minval(x(3,:)))
      ! print*, "Lx=", Lx
      ! print*, "Ly=", Ly
      ! if(nsd==3) print*, "Lz=", Lz

      !!!!!!!!!!! Estimativa da Permeabiliba por Teste de Poço !!!!!!!!!!!!!!!!!!!!!!

      permMatriz=campoPermeabilidade(1,1)
      Pinicial=estrutSistEqP%uInicial
      Ppoco=0.d0
      beta=c(6,1)
      
      tempoAdimensionalAnterior=(permMatriz*(tempo-deltaT))/(viscosidade*beta*Ly*Ly)
      tempoAdimensionalAtual=(permMatriz*tempo)/(viscosidade*beta*Ly*Ly)
      
      vazaoAdimensionalAnterior=(Ly*vazaoAnterior(1)*viscosidade)/(permMatriz*Lx*(Pinicial-Ppoco))
      vazaoAdimensionalAtual=(Ly*vazao(1)*viscosidade)/(permMatriz*Lx*(Pinicial-Ppoco))

      mltil=(vazao(1)-vazaoAnterior(1))/(1/sqrt(tempo)-1/sqrt(tempo-deltaT))
      Kestimado=(pi*viscosidade*mltil**2)/(beta*Lx**2*(Pinicial-Ppoco)**2)

      ! write(658,'(5es20.8)') tempoAdimensionalAtual, vazaoAdimensionalAtual
      ! write(659,'(5es20.8)') tempoAdimensionalAtual, (log(vazaoAdimensionalAtual)-log(vazaoAdimensionalAnterior))/(log(tempoAdimensionalAtual)-log(tempoAdimensionalAnterior)), Kestimado

      !! GradP = (p2-p1)/L   
      gradPx=(pressaoEsq-pressaoDir)/Lx
      gradPy=(pressaoFrente-pressaoAtras)/Ly
      write(*,'(a,1es20.8)') "gradPx:", gradPx
      write(*,'(a,1es20.8)') "gradPy:", gradPy
      if(nsd==3) then
         gradPz=(pressaoInf-pressaoSup)/Lz
         write(*,'(a,1es20.8)') "gradPz:", gradPz
      endif

      vm_x_dir    = vazao(2)/area(2)
      vm_x_esq    = vazao(4)/area(4)
      vm_y_atras  = vazao(3)/area(3)
      vm_y_frente = vazao(1)/area(1)
      if(nsd==3) then
         vm_z_sup = vazao(6)/area(6)
         vm_z_inf = vazao(5)/area(5)    
      endif
      !TODO: Acho que a vazão esta calculada errada, deveria ser feito a regularização da forma:
      ! vm_x_dir    = vazao(2)/area(1)
      ! vm_x_esq    = vazao(4)/area(3)
      ! vm_y_atras  = vazao(3)/area(4)
      ! vm_y_frente = vazao(1)/area(2)

      write(*,'(a, 6e15.5)') "Vel. Medias", vm_x_dir, vm_x_esq, vm_y_atras, vm_y_frente, vm_z_inf, vm_z_sup

   end subroutine calcularVazao









































   ! --------------------------------------------------------------------------------
   subroutine calcVazao(estrutSistEqP, Lx, Ly, Lz)
      ! Rotina para calcular a vazão nos bordos e os gradientes de pressão em relação aos eixos

      use mGlobaisEscalares, only: npintFratura, nrowsh, tempo, deltaT, pi, viscosidade
      use mGlobaisArranjos,  only: tiposElemBordas, aberturasBorda, campoPermeabilidade, vazao
      use mGlobaisArranjos,  only: c, matFratura, tiposElemBordas
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shgq, shlq, shgCodimOneSurface_points
      use mfuncoesDeForma,   only: shl1D_intPoints, shlqen, shlten
      use mMalha,            only: local, x, nenFratura, nsd, numelBordas, conecNodaisBordas
      use mMalha,            only: numelFratura, normal, conectividadesInterseccao
      use mMLayer,           only: LayerMaterial

      type (estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqP
          
      real*8 :: xl(nsd,nenFratura)
      real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura), shlParaArea(nsd,nenFratura,npintFratura)
      real*8 :: det(npintFratura), w(npintFratura)
      real*8  :: wDet, di(nrowsh)
      integer*4:: nee, nelBordas, tipoBorda, l, i, indexI, sd_i, ndofDesloc      
      real*8  :: xCentroElem(nsd), shlCentroElem(nsd,nenFratura,1), wCentroElem(1), vetorParaNo(nsd), normalNoLocal(nsd,nenFratura)
      real*8  :: area(6), vn_ptoInt
      real*8  :: pressaoDir, pressaoEsq, pressaoInf, pressaoSup, pressaoFrente, pressaoAtras
      real*8  :: diffP, Lx, Ly, Lz, distanciaEntrePontos(3), distancia, areaFratura
      real*8  :: tempoAdimensionalAtual, tempoAdimensionalAnterior, vazaoAdimensionalAtual, vazaoAdimensionalAnterior
      real*8 :: vazaoAnterior(14), mltil, Kestimado, Ppoco, Pinicial, permMatriz, beta

      real*8, dimension(:), allocatable :: layDh
      integer*4:: ncamadas, k

      w=0.0
      shl=0.0
      wCentroElem=0d0
      ! shlCentroElem=0d0
      vazaoAnterior(:)=vazao(:)

      if(nenFratura==2) then           ! Elemento(Linha) - Espaco(2D)
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
         ! call shl1D_intPoints(shlCentroElem,wCentroElem,1,nenFratura)
      else if(nenFratura==3) then      ! Elemento(triangulo) - Espaco(3D)
         ! CALL SHLTEN(SHL,nenFratura)
         call shlt(shl,w,npintFratura,nenFratura)
         ! call shlt(shlCentroElem,wCentroElem,1,nenFratura)
      else if(nenFratura==4) then      ! Elemento(quadrado) - Espaco(3D)
         call shlq(shl,w,npintFratura,nenFratura)
         CALL SHLQEN(SHL,nenFratura)
         ! call shlq(shlCentroElem,wCentroElem,1,nenFratura)
      end if

      area=0.d0
      vazao=0.d0

      ! Propriedades da camada 1 e 3
      ! para armazenar as propridades na falha
      ncamadas = estrutSistEqP%ndof-2
      allocate(layDh(ncamadas))
      do k=1, ncamadas
         layDh(k) = c(2,2)
      enddo

      do nelBordas=1,numelBordas

         ! Armazena em xl os valores das posições (x) do elemento         
         call local(conecNodaisBordas(:,nelBordas), x, xl, nenFratura, nsd, nsd)         
         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)

         !... Identificacao do tipo de borda
         ! 1: Inferior
         ! 2: Direita
         ! 3: Superior
         ! 4: Esquerda
         tipoBorda = tiposElemBordas(nelBordas)
         if(tipoBorda==1) then
            pressaoFrente=estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==2) then
            pressaoDir=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==3) then
            pressaoAtras= estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==4) then
            pressaoEsq=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==5) then
            pressaoInf=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         else 
            pressaoSup=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         endif
         area(tipoBorda) = area(tipoBorda) + dot_product(w(:),det(:))

         !.... loop on integration points         
         ! Computo da vazao para o macico rochoso
         do l=1,npintFratura
            !velocidade normal avaliada no ponto de integracao! !!!ou será no nó?
           vn_ptoInt = dot_product(velocidadeBordas(:,nelBordas), shg(nrowsh,:,l))
           vazao(tipoBorda) = vazao(tipoBorda) + vn_ptoInt*w(l)*det(l)
         end do

         ! Computo da vazao para a descontinuidade
         if(nsd==2) then      ! Espaco(2D)
            do i=1, nenfratura
               if(aberturasborda(i,nelbordas)>1.d-10) then
                  do k=1, ncamadas
                     area(tipoborda)=area(tipoborda)+layDh(k)/2.d0
                     vazao(tipoborda)=vazao(tipoborda)+velocidadeBordasFalhas(i,nelbordas,k)*layDh(k)/2.d0
                  end do
                  ! area(tipoborda)=area(tipoborda)+aberturasborda(i,nelbordas)/2.d0
                  ! vazao(tipoborda)=vazao(tipoborda)+velocidadebordasfratura(i,nelbordas)*aberturasborda(i,nelbordas)/2.d0
               end if
            end do
         else                 ! Espaco(3D)
            if((conectividadesInterseccao(1,nelBordas).ne.0)) then   ! Caso com intersecao de descontinuidade
               distanciaEntrePontos(:)=(x(:,conectividadesInterseccao(1,nelBordas)) - x(:,conectividadesInterseccao(2,nelBordas)))
               distancia=sqrt(dot_product(distanciaEntrePontos(:),distanciaEntrePontos(:)))
               
               do i=1, nenFratura
                  if(aberturasBorda(i,nelBordas)>1.d-10) exit
               enddo
               areaFratura=(aberturasBorda(i,nelBordas)/2.d0)*distancia
               area(tipoBorda)=area(tipoBorda)+areaFratura
               vazao(tipoBorda)=vazao(tipoBorda)+velocidadeBordasFratura(i,nelBordas)*areaFratura
            endif
         endif
         
      end do

      deallocate(layDh)

      ! write(657,'(5es20.8)') tempo, abs(vazao(1:4))
      write(*,'(a, 6e15.5)') "Areas", area
      
      write(*,'(a,1es20.8)') "Vazao na borda da frente   = ", abs(vazao(1))
      write(*,'(a,1es20.8)') "Vazao na borda da direita  = ", abs(vazao(2))
      write(*,'(a,1es20.8)') "Vazao na borda de trás     = ", abs(vazao(3))
      write(*,'(a,1es20.8)') "Vazao na borda da esquerda = ", abs(vazao(4))
      if(nsd==3) then
         write(*,'(a,1es20.8)') "Vazao na borda inferior    = ",abs(vazao(5))
         write(*,'(a,1es20.8)') "Vazao na borda superior    = ",abs(vazao(6))
      endif
   
      Lx=(maxval(x(1,:)) - minval(x(1,:)))
      Ly=(maxval(x(2,:)) - minval(x(2,:)))
      if(nsd==3) Lz=(maxval(x(3,:)) - minval(x(3,:)))

      !!!!!!!!!!! Estimativa da Permeabiliba por Teste de Poço !!!!!!!!!!!!!!!!!!!!!!
      permMatriz=campoPermeabilidade(1,1)
      Pinicial=estrutSistEqP%uInicial
      Ppoco=0.d0
      beta=c(6,1)
      
      tempoAdimensionalAnterior=(permMatriz*(tempo-deltaT))/(viscosidade*beta*Ly*Ly)
      tempoAdimensionalAtual=(permMatriz*tempo)/(viscosidade*beta*Ly*Ly)
      
      vazaoAdimensionalAnterior=(Ly*vazaoAnterior(1)*viscosidade)/(permMatriz*Lx*(Pinicial-Ppoco))
      vazaoAdimensionalAtual=(Ly*vazao(1)*viscosidade)/(permMatriz*Lx*(Pinicial-Ppoco))

      mltil=(vazao(1)-vazaoAnterior(1))/(1/sqrt(tempo)-1/sqrt(tempo-deltaT))
      Kestimado=(pi*viscosidade*mltil**2)/(beta*Lx**2*(Pinicial-Ppoco)**2)


      gradPx=(pressaoEsq-pressaoDir)/Lx
      gradPy=(pressaoFrente-pressaoAtras)/Ly
      write(*,'(a,1es20.8)') "gradPx:", gradPx
      write(*,'(a,1es20.8)') "gradPy:", gradPy
      if(nsd==3) then
         gradPz=(pressaoInf-pressaoSup)/Lz
         write(*,'(a,1es20.8)') "gradPz:", gradPz
      endif

      vm_x_dir    = vazao(2)/area(2)
      vm_x_esq    = vazao(4)/area(4)
      vm_y_atras  = vazao(3)/area(3)
      vm_y_frente = vazao(1)/area(1)
      if(nsd==3) then
         vm_z_sup = vazao(6)/area(6)
         vm_z_inf = vazao(5)/area(5)    
      endif

   end subroutine calcVazao

   !=======================================================================
   !------- ROTINAS PARA CALCULO DAS PROPRIEDADES EQUIVALENTES ------------
   !=======================================================================    
   subroutine calcularPropsEfetivas()
      use mGlobaisEscalares, only: viscosidade
      use mGlobaisEscalares, only: areaCelula, areaFluido, porosidadeCelula_ref
      use mMalha,            only: nsd

      real*8 :: perm_x_efetiva_esquerda, perm_x_efetiva_direita
      real*8 :: perm_y_efetiva_frente,   perm_y_efetiva_atras
      real*8 :: perm_z_efetiva_inferior, perm_z_efetiva_superior
      double precision, dimension(6) :: out
      integer :: i

      perm_x_efetiva_esquerda = viscosidade*vm_x_esq/gradPx
      perm_x_efetiva_direita  = viscosidade*vm_x_dir/gradPx
          
      perm_y_efetiva_frente = viscosidade*vm_y_frente/gradPy
      perm_y_efetiva_atras  = viscosidade*vm_y_atras/gradPy
          
      if(nsd==3) then
         perm_z_efetiva_inferior = viscosidade*vm_z_inf/gradPz
         perm_z_efetiva_superior = viscosidade*vm_z_sup/gradPz
      endif
  
      ! Propriedades efetivas
      write(*,*) "Kx(esq)   =", abs(perm_x_efetiva_esquerda)
      write(*,*) "Kx(dir)   =", abs(perm_x_efetiva_direita)
      write(*,*) "Ky(frente)=", abs(perm_y_efetiva_frente)
      write(*,*) "Ky(atras) =", abs(perm_y_efetiva_atras)
      if(nsd==3) then
         write(*,*) "Kz(inf)   =", abs(perm_z_efetiva_inferior)
         write(*,*) "Kz(sup)   =", abs(perm_z_efetiva_superior)
      endif

      out(1)=abs(perm_x_efetiva_esquerda)
      out(2)=abs(perm_x_efetiva_direita)
      out(3)=abs(perm_y_efetiva_frente)
      out(4)=abs(perm_y_efetiva_atras)
      if(nsd==3) then
         out(5)=abs(perm_z_efetiva_inferior)
         out(6)=abs(perm_z_efetiva_superior)
      else
         out(5)=0.
         out(6)=0.
      endif
      open(unit=756, action='write', position='append', file= './out/res2.inc')
      write(756,'(6es20.8)')  (out(i),i=1,6)
      close(756)


   end subroutine calcularPropsEfetivas


   !=======================================================================
   !---------- ROTINAS PARA CALCULO DE FLUXO CONSISTENTE ------------------
   !=======================================================================

   subroutine montarSistEqAlgVelocidadeNormal()
      use mMalha,            only: x, numel, nen, numnp, nsd
      use mMalha,            only: numelFratura, nenFratura, numnpFratura
      use mMalha,            only: conecNodaisElem, conecNodaisFratura, nenFraturasNaBorda
      use mGlobaisEscalares, only: npintFraturasNaBorda, nrowshFratNaBorda
      use mPotencial,        only: estrutSistEqP

      call calcCoefSistAlgFluxoDirichletBRHS( estrutSistEqFluxoNormal, estrutSistEqP, x, conecNodaisElem, numnp, numel, nen, nsd )
      call calcCoefSistAlgFluxoDirichletALHS( estrutSistEqFluxoNormal)

      if (numelFratura>0) then 
         nenFraturasNaBorda = 2
         npintFraturasNaBorda = 2
         nrowshFratNaBorda = 2
         call calcCoefSistAlgFluxoDirichlet_DesconBRHS( estrutSistEqFluxoNormal, estrutSistEqP, &
            x,  conecNodaisFratura, numnpFratura, numelFratura, nenFratura, nsd )
         if (nsd == 2 ) call calcCoefSistAlgFluxoDirichlet_DesconALHS_2D(estrutSistEqFluxoNormal)
         if (nsd == 3 ) call calcCoefSistAlgFluxoDirichlet_DesconALHS   (estrutSistEqFluxoNormal)
      endif   
   end subroutine


   !=======================================================================
   subroutine calcCoefSistAlgFluxoDirichletALHS( estrutSistEqFD_ )
      use mGlobaisEscalares, only: transiente, passoTempo, npintFratura, viscosidade, nrowsh
      use mGlobaisArranjos,  only: tiposElemBordas, bordaPrescrita
      use mfuncoesDeForma,   only: shlt, shlq, shl1D_intPoints, shgCodimOneSurface_points
      use mMalha,            only: x, nsd, local, nenFratura, numelBordas, conecNodaisBordas

      use mUtilSistemaEquacoes,  only: kdbc,KDBC2

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre
      
      implicit none

      type (estruturasArmazenamentoSistemaEq), intent(inout) :: estrutSistEqFD_
      
      real*8 :: xl(nsd,nenFratura)
      real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)
      real*8 :: det(npintFratura), w(npintFratura)
      
      integer*4 :: nee

      real*8   :: eleffm(estrutSistEqFD_%ndof*nenFratura,estrutSistEqFD_%ndof*nenFratura)
      
      integer*4:: nel, l, i, j

      real*8  :: wDet
      logical :: diag,quad,lsym
      real*8  :: di(nrowsh), dj(nrowsh)
      integer*4 :: indexI, indexJ
      integer :: noGlobal

      nee  = nenFratura*estrutSistEqFD_%ndof
      diag = .false.

      w=0.0
      shl=0.0

      if(nenFratura==2) then
         call shl1D_intPoints( shl, w, npintFratura, nenFratura )
      else if(nenFratura==3) then
         call shlt( shl, w, npintFratura, nenFratura )
      else if(nenFratura==4) then
         call shlq( shl, w, npintFratura, nenFratura )
      end if          

      do nel=1, numelBordas
         if( bordaPrescrita(tiposElemBordas(nel)) ) then   
            ! clear stiffness matrix and force array
            eleffm=0.0

            ! LOCALIZE COORDINATes and Dirichlet b.c.
            call local( conecNodaisBordas(:,nel), x, xl, nenFratura, nsd, nsd)
            quad = .true. 
            call shgCodimOneSurface_points( xl, det, shl, shg, nenFratura, npintFratura, nsd)
         
            !.... loop on integration points
            do l=1,npintFratura
               wDet = w(l)*det(l)
               do j=1,nenFratura
                  dj(:) = shg(:,j,l)*wDet
                  indexJ = (j-1)*estrutSistEqFD_%ndof+1
                  noGlobal=conecNodaisBordas(j,nel)
                 
                  do i=1,nenFratura
                     di(:) = shg(:,i,l)
                     indexI = (i-1)*estrutSistEqFD_%ndof+1
                     if(transiente) then
                        eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ di(nrowsh)*dj(nrowsh)
                     else
                        eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ di(nrowsh)*dj(nrowsh)
                     endif            
                  end do ! i                      
               enddo ! j                   
            end do ! ponto integracao
                     
            !.... assemble element stifness matrix and force array into global
            !       left-hand-side matrix and right-hand side vector
            lsym=.true.         

            if (passoTempo==1) then
               if (estrutSistEqFD_%optSolver=='skyline')   then
                  call addlhsN( estrutSistEqFD_%alhs, eleffm, estrutSistEqFD_%lmBorda, estrutSistEqFD_%idiag, &
                           nee, nel, diag, lsym)
               endif
               if (estrutSistEqFD_%optSolver=='pardiso') then
                  call addlhsCSR( estrutSistEqFD_,eleffm, estrutSistEqFD_%lmBorda, nee, nel)
               endif
               if (estrutSistEqFD_%optSolver=='hypre')   then
                  ! call addnslHYPRE( estrutSistEqFD_%A_HYPRE, eleffm, estrutSistEqFD_%lmBorda, nee, nel)
               endif
            endif
         endif
      end do

   end subroutine calcCoefSistAlgFluxoDirichletALHS

   !=======================================================================
   subroutine calcCoefSistAlgFluxoDirichletBRHS( estrutSistEqFD_, estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
      ! numel -> numero de elementos no dominio.
      ! conecNodaisElem -> conectividades do dominio.      
      use mGlobaisEscalares, only: nrowsh, npint, transiente, viscosidade, deltaT, saltoPressao
      use mGlobaisArranjos,  only: mat, c, campoPermeabilidade
      use mfuncoesDeForma,   only: oneshg, shlt, shlq3d, shg3d, shgq, shlq, shlt3D
      use mMalha,            only: local, normal, numelFratura

      use mUtilSistemaEquacoes,  only: kdbc2

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre
      
      implicit none

      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqFD_, estrutSistEqP_
      real*8,    intent(in)    :: x(nsd, numnp)
      integer*4, intent(in)    :: numnp, numel, nen, nsd
      integer*4, intent(in)    :: conecNodaisElem(nen,numel)
      
      real*8 :: xl(nsd,nen), pl(estrutSistEqP_%ndof,nen), fl(estrutSistEqFD_%ndof,nen)
      real*8 :: plAnt(estrutSistEqP_%ndof,nen)
      real*8 :: shg(nrowsh,nen,npint), shl(nrowsh,nen,npint)
      real*8 :: det(npint), w(npint), vetorParaNo(nsd), normalNoLocal(nsd,nen)
      real*8 :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1)
      
      integer*4 :: nee
      real*8    :: elresf(estrutSistEqFD_%ndof*nen)
      
      integer*4 :: nel, m, l, j, indexJ, noGlobal, i
      real*8  :: permeabilidadeX, permeabilidadeY, permeabilidadeZ
      real*8  :: condutividadeX, condutividadeY, condutividadeZ, beta
      real*8  :: wDet
      real*8  :: pAtual, pAnterior
      logical :: diag,zerodl,quad
      
      real*8  :: dj(nrowsh)
      
      real*8  :: GRADPX, GRADPY, GRADPZ
      integer :: dof
 
      nee  = nen*estrutSistEqFD_%ndof
      diag = .false.
      
      w=0.0
      shl=0.0
      
      ! dof = estrutSistEqP_%ndof
      dof = estrutSistEqFD_%ndof

      ! if(nen==3) then
      !    call shlt(shl,w,npint,nen)
      ! else if(nen==4.and.nsd==2) then
      !    call shlq(shl,w,npint,nen)
      ! else if(nen==4.and.nsd==3) then
      !    call shlt3D(shl,w,npint,nen)
      ! else if(nen==8) then
      !    call shlq3d(shl,w,npint,nen)
      ! end if

      wCentroElem=0d0
      shlCentroElem=0d0
        
      if(nen==3) then
         call shlt(shl,w,npint,nen)
         call shlt(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==2) then
         call shlq(shl,w,npint,nen)
         call shlq(shlCentroElem,wCentroElem,1,nen)
      else if(nen==4.and.nsd==3) then
         call shlt3D(shl,w,npint,nen)
         call shlt3D(shlCentroElem,wCentroElem,1,nen)
      else if(nen==8) then
         call shlq3d(shl,w,npint,nen)
         call shlq3D(shlCentroElem,wCentroElem,1,nen)
      end if

      do nel=1,numel
         ! clear force array
         elresf=0.0

         ! LOCALIZE COORDINATes and Dirichlet b.c.
         call local(conecNodaisElem(1,nel), x, xl, nen, nsd, nsd)
         call local(conecNodaisElem(1,nel), estrutSistEqFD_%u, fl, nen, estrutSistEqFD_%ndof, estrutSistEqFD_%ndof)
         call local(conecNodaisElem(1,nel), estrutSistEqP_%u,  PL, nen, estrutSistEqP_%ndof,  estrutSistEqP_%ndof)

         ! print*, "nel", nel
         ! print*, "PL-", pl
         ! print*, "FL-", fl
         ! print*, "-------------------------------"
         
         if (numelFratura > 0) then
            call local(conecNodaisElem(:,nel), normal, normalNoLocal, nen, nsd, nsd)

            do i=1,nsd
               xCentroElem(i) = dot_product(xl(i,:), shlCentroElem(nrowsh,:,1))
            end do
         end if
         
         if(transiente.eqv..true.) then
            call local(conecNodaisElem(1,nel),estrutSistEqP_%uTempoAnt,plAnt,nen,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
         end if          

         quad = .true.
         if (nen.eq.4) then 
            if (conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
         end if             

         if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,1,nel,1)
         if(nen==3) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
         if(nen==4.and.nsd==2) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
         if(nen==4.and.nsd==3) call shg3d (xl,det,shl,shg,npint,nel,nen)
         if(nen==8) call shg3d (xl,det,shl,shg,npint,nel,nen)
         
         !....... form stiffness matrix         
         !... Properties         
         m = mat(nel)
         permeabilidadeX = campoPermeabilidade(1,nel) !c(3,m)
         permeabilidadeY = campoPermeabilidade(2,nel) !c(4,m)
         if (nsd==3) permeabilidadeZ = campoPermeabilidade(3,nel) !c(5,m)
         beta=c(6,m) 
         ! print*,"permX=",permeabilidadeX,"permY=",permeabilidadeY
         ! print*, "visco=",viscosidade

         condutividadeX = permeabilidadeX/viscosidade
         condutividadeY = permeabilidadeY/viscosidade
         condutividadeZ = permeabilidadeZ/viscosidade

         !.... loop on integration points
         do l=1,npint
            wDet = w(l)*det(l)
            
            if(transiente.eqv..true.) then
               pAnterior = 0.d0
               pAtual    = 0.d0
               DO J=1,NEN
                  pAtual    = pAtual    + SHG(nrowsh,J,L)*pl(1,J)
                  pAnterior = pAnterior + SHG(nrowsh,J,L)*plAnt(1,J)
               ENDDO
               ! if ( nel .LE. 100 ) print*, "pAtual, pAnterior", nel,l, pAtual, pAnterior, pAtual - pAnterior
            endif
            
            GRADPX=0.0
            GRADPY=0.0
            if(nsd==3)GRADPZ=0.0

            do j=1,nen
               !! NOVO
               dof = 1
               if ( (numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) then
                  if ((estrutSistEqP_%lm(2,J,nel) > 0) .or. (estrutSistEqP_%lm(3,J,nel) > 0)) then
                     vetorParaNo(:) = xl(:,J) - xCentroElem(:)                   
                     if (dot_product(normalNoLocal(:,J), vetorParaNo(:)) < 0d0) then
                        dof = 2
                     end if
                  end if
               end if   
            
               GRADPX = GRADPX + shg(1,J,L)*PL(dof,J)
               GRADPY = GRADPY + shg(2,J,L)*PL(dof,J)
               if(nsd==3) GRADPZ = GRADPZ + shg(3,J,L)*PL(dof,J)
            end do
            
            do j=1,nen
               noGlobal=conecNodaisElem(j,nel)
               
               dj(:) = shg(:,j,l)*wDet
               ! indexJ = (j-1)*estrutSistEqP_%ndof + 1
               indexJ = (j-1)*estrutSistEqFD_%ndof + 1                 

               if(estrutSistEqFD_%id(1,noGlobal).ne.0) then
                  if(transiente.eqv..true.) then
                     ! print*, "pAtual, pAnterior", nel,l,j, noGlobal, pAtual, pAnterior, pAtual - pAnterior                     
                     elresf(indexJ)=elresf(indexJ) - beta * ((pAtual-pAnterior)/deltaT) * dj(nrowsh)
                  endif      
                  elresf(indexJ)=elresf(indexJ) - (dj(1)*condutividadeX*GRADPX + dj(2)*condutividadeY*GRADPY)
                  if(nsd==3) elresf(indexJ)=elresf(indexJ) - dj(3)*condutividadeZ*GRADPZ
               endif
            end do
         end do

         ! print*,"nel", nel, "estrutSistEqFD_%lm", estrutSistEqFD_%lm(:,:,nel)
         call addrhsN (estrutSistEqFD_%brhs, elresf, estrutSistEqFD_%lm, nee, nel)
      end do


      return
   end subroutine calcCoefSistAlgFluxoDirichletBRHS

   !=======================================================================
   subroutine  calcCoefSistAlgFluxoDirichlet_DesconBRHS(estrutSistEqFD_, estrutSistEqP_, x, conecNodaisElem, &
      numnpFratura, numelFratura, nenFratura, nsd )
      
      use mGlobaisEscalares, only: nrowsh, npintFratura, transiente, viscosidade, deltaT, saltoPressao
      use mGlobaisArranjos,  only: matFratura, c
      use mfuncoesDeForma,   only: shlt,shlq, shgCodimOneSurface_points, shl1D_intPoints
      use mMalha,            only: local

      use mMalha,            only: conecNodaisReduzido

      use mUtilSistemaEquacoes,  only: kdbc2

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre
      use mMLayer,           only: LayerMaterial
      
      implicit none

      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqFD_, estrutSistEqP_
      real*8,    intent(in)    :: x(nsd, numnpFratura)
      integer*4, intent(in)    :: numnpFratura, numelFratura, nenFratura, nsd
      integer*4, intent(in)    :: conecNodaisElem(nenFratura,numelFratura)
      
      real*8 :: xl(nsd,nenFratura), fl(estrutSistEqFD_%ndof,nenFratura)
      real*8 :: pl(estrutSistEqP_%ndof,nenFratura), plAnt(estrutSistEqP_%ndof,nenFratura)
      real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)
      real*8 :: det(npintFratura), w(npintFratura)
      
      integer*4 :: nee
      real*8    :: elresf(estrutSistEqFD_%ndof*nenFratura)
      
      integer*4 :: nel, m, l, j, indexJ, noGlobal
      real*8  :: permeabilidade, condutividade, abertura, beta
      real*8  :: wDet
      logical :: zerodl,quad,diag
      real*8  :: pAtual, pAnterior, dj(nrowsh)
      
      ! real*8 :: GRADPX, GRADPY, GRADPZ
      
      integer :: dof
      real*8  :: gradP
      integer*4:: k, sp, ncamadas
      real*8, dimension(:), allocatable :: layKN, layKT, layDh, layBeta


      nee  = nenFratura*estrutSistEqFD_%ndof
      diag = .false.
      
      w=0.0
      shl=0.0
      
      dof = estrutSistEqFD_%ndof
      
      if(nenFratura==2) then
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then
         call shlq(shl,w,npintFratura,nenFratura)
      end if
      
      do nel=1,numelFratura
         elresf=0.0
         
         ! LOCALIZE COORDINATes and Dirichlet b.c.
         call local(conecNodaisElem(1,nel), x, xl, nenFratura, nsd, nsd)
         call local(conecNodaisElem(1,nel), estrutSistEqFD_%u, fl, nenFratura, estrutSistEqFD_%ndof, estrutSistEqFD_%ndof)
         call local(conecNodaisElem(1,nel), estrutSistEqP_%u,  PL, nenFratura, estrutSistEqP_%ndof,  estrutSistEqP_%ndof)
         
         if(transiente.eqv..true.) then
            call local(conecNodaisElem(1,nel),estrutSistEqP_%uTempoAnt,plAnt,nenFratura,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
         end if          

         quad = .true.
         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)            
         
         !....... form stiffness matrix
         !... Properties
         m = matFratura(nel)

         ! Propriedades das camadas
         ncamadas = LayerMaterial(m,1)
         allocate(layDh(ncamadas))
         allocate(layKN(ncamadas))
         allocate(layKT(ncamadas))
         allocate(layBeta(ncamadas))
         do k=1, ncamadas
            layDh(k) = c(2,LayerMaterial(m,k+1))
            layKN(k) = c(3,LayerMaterial(m,k+1))/viscosidade
            layKT(k) = c(4,LayerMaterial(m,k+1))/viscosidade
            layBeta(k) = c(6,LayerMaterial(m,k+1))
         enddo
         
         !.... loop on integration points
         do l=1,npintFratura
            wDet = w(l)*det(l)
            
            if(transiente.eqv..true.) then
               pAnterior = 0.d0
               pAtual    = 0.d0
               DO j=1,nenFratura
                  pAtual    = pAtual    + SHG(nrowsh,J,L)*PL(1,J)
                  pAnterior = pAnterior + SHG(nrowsh,J,L)*plAnt(1,J)
               ENDDO
            endif
            
            !.... BEGIN LOOP OVER EACH LAYER FOR THE ELEMENT
            do k=1,ncamadas
               ! COMPUTE THE DERIVATIVES IN FUNCTION OF THE SPACE DIMENSION
               do SP=1,nsd
                  
                  gradP=0d0
                  ! COMPUTE TANGENCIAL DERIVATIVES OF THE NODAL DISPLACEMENT
                  do j=1,nenFratura
                     if ((estrutSistEqP_%lmFratura(2,J,nel)==0).or.(estrutSistEqP_%lmFratura(3,J,nel)==0))then
                        dof=1
                     else
                        dof=2+k
                     endif
                     gradP=gradP+SHG(sp,j,L)*PL(dof,J)
                     ! print*,"pl(",dof,",",J,")=",PL(dof,J)
                  enddo

                  do j=1,nenFratura
                     noGlobal=conecNodaisElem(j,nel)
                     
                     dj(:) = shg(:,j,l)*wDet
                     indexJ = (j-1)*estrutSistEqFD_%ndof + k    
      
                     if(estrutSistEqFD_%id(1,noGlobal).ne.0) then
                        if(transiente.eqv..true.) then
                           elresf(indexJ)=elresf(indexJ) - layDh(k)*layBeta(k)*((pAtual-pAnterior)/deltaT)*dj(nrowsh)
                        endif      
                        elresf(indexJ)=elresf(indexJ) - layDh(k)*layKN(k)*gradP*dj(sp)
                     endif
                  end do
               enddo
            enddo
         end do

         ! computation of Dirichlet b.c. contribution
         call ztest(fl,nee,zerodl)
         call addrhsN (estrutSistEqFD_%brhs, elresf, estrutSistEqFD_%lmFratura, nee, nel)

         deallocate(layDh)
         deallocate(layKN)
         deallocate(layKT)
         deallocate(layBeta)
      end do

   end subroutine calcCoefSistAlgFluxoDirichlet_DesconBRHS   

   ! !=======================================================================
   ! subroutine calcCoefSistAlgFluxoDirichlet_DesconBRHS(estrutSistEqFD_, estrutSistEqP_, x, conecNodaisElem, &
   !    numnpFratura, numelFratura, nenFratura, nsd )
      
   !    use mGlobaisEscalares, only: nrowsh, npintFratura, transiente, viscosidade, deltaT, saltoPressao
   !    use mGlobaisArranjos,  only: matFratura, c
   !    use mfuncoesDeForma,   only: shlt,shlq, shgCodimOneSurface_points, shl1D_intPoints
   !    use mMalha,            only: local

   !    use mUtilSistemaEquacoes,  only: kdbc2

   !    use mSolverGaussSkyline,   only: addrhsN, addlhsN
   !    use mSolverPardiso,        only: addlhsCSR
   !    use mSolverHypre
   !    use mMLayer,           only: LayerMaterial
      
   !    implicit none

   !    type (estruturasArmazenamentoSistemaEq) :: estrutSistEqFD_, estrutSistEqP_
   !    real*8,    intent(in)    :: x(nsd, numnpFratura)
   !    integer*4, intent(in)    :: numnpFratura, numelFratura, nenFratura, nsd
   !    integer*4, intent(in)    :: conecNodaisElem(nenFratura,numelFratura)
      
   !    real*8 :: xl(nsd,nenFratura), fl(estrutSistEqFD_%ndof,nenFratura)
   !    real*8 :: pl(estrutSistEqP_%ndof,nenFratura), plAnt(estrutSistEqP_%ndof,nenFratura)
   !    real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)
   !    real*8 :: det(npintFratura), w(npintFratura)
      
   !    integer*4 :: nee
   !    real*8    :: elresf(estrutSistEqFD_%ndof*nenFratura)
      
   !    integer*4 :: nel, m, l, j, k, indexJ, noGlobal
   !    real*8  :: permeabilidade, condutividade, abertura, beta
   !    real*8  :: wDet
   !    logical :: diag,zerodl,quad
   !    real*8  :: pAtual, pAnterior, dj(nrowsh)
      
   !    real*8 :: GRADPX, GRADPY, GRADPZ
      
   !    integer :: dof
   !    integer*4:: ncamadas, posI1, posI2, posI3, posJ1, posJ2, posJ3
   !    real*8, dimension(:), allocatable :: layKN, layKT, layDh, layBeta


   !    nee  = nenFratura*estrutSistEqFD_%ndof
   !    diag = .false.
      
   !    w=0.0
   !    shl=0.0
      
   !    dof = estrutSistEqFD_%ndof
      
   !    if(nenFratura==2) then
   !       call shl1D_intPoints(shl, w, npintFratura, nenFratura)
   !    else if(nenFratura==3) then
   !       call shlt(shl,w,npintFratura,nenFratura)
   !    else if(nenFratura==4) then
   !       call shlq(shl,w,npintFratura,nenFratura)
   !    end if

   !    print*, 'Graus de Liberdade', estrutSistEqP_%ndof
      
   !    do nel=1,numelFratura
   !       elresf=0.0
         
   !       ! LOCALIZE COORDINATes and Dirichlet b.c.
   !       call local(conecNodaisElem(1,nel), x, xl, nenFratura, nsd, nsd)
   !       call local(conecNodaisElem(1,nel), estrutSistEqFD_%u, fl, nenFratura, estrutSistEqFD_%ndof, estrutSistEqFD_%ndof)
   !       call local(conecNodaisElem(1,nel), estrutSistEqP_%u,  PL, nenFratura, estrutSistEqP_%ndof,  estrutSistEqP_%ndof)
         
   !       if(transiente.eqv..true.) then
   !          call local(conecNodaisElem(1,nel),estrutSistEqP_%uTempoAnt,plAnt,nenFratura,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
   !       end if          

   !       quad = .true.
   !       call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)            
         
   !       !....... form stiffness matrix
   !       !... Properties
   !       m = matFratura(nel)
   !       abertura = c(2,m)
   !       permeabilidade = c(3,m)
   !       condutividade = permeabilidade/viscosidade
   !       beta=c(6,m)
   !       ! Propriedades das camadas
   !       ncamadas = LayerMaterial(m,1)
   !       allocate(layDh(ncamadas))
   !       allocate(layKN(ncamadas))
   !       allocate(layKT(ncamadas))
   !       allocate(layBeta(ncamadas))
   !       do k=1, ncamadas
   !          layDh(k) = c(2,LayerMaterial(m,k+1))
   !          layKN(k) = c(3,LayerMaterial(m,k+1))/viscosidade
   !          layKT(k) = c(4,LayerMaterial(m,k+1))/viscosidade
   !          layBeta(k) = c(6,LayerMaterial(m,k+1))
   !       enddo
         
   !       !.... loop on integration points
   !       do l=1,npintFratura
   !          wDet = w(l)*det(l)
            
   !          if(transiente.eqv..true.) then
   !             pAnterior = 0.d0
   !             pAtual    = 0.d0
   !             DO j=1,nenFratura
   !                pAtual    = pAtual    + SHG(nrowsh,J,L)*PL(1,J)
   !                pAnterior = pAnterior + SHG(nrowsh,J,L)*plAnt(1,J)
   !             ENDDO
   !          endif
            
   !          GRADPX=0.0
   !          GRADPY=0.0
   !          if(nsd==3)GRADPZ=0.0
            
   !          do j=1,nenFratura
   !             ! do k=1,ncamadas
   !             ! enddo
   !             !! NOVO
   !             if (  saltoPressao .EQV. .true. ) then
   !                if ((estrutSistEqP_%lmFratura(2,J,nel)==0).or.(estrutSistEqP_%lmFratura(3,J,nel)==0))then
   !                   dof=1
   !                else
   !                   dof = ncamadas+2
   !                   ! dof=3
   !                endif
   !             end if                
   !             ! print*, 'Salto Pressão', saltoPressao
   !             ! print*, 'Graus de liberdade', dof

   !             GRADPX = GRADPX + shg(1,J,L)*PL(dof,J)
   !             GRADPY = GRADPY + shg(2,J,L)*PL(dof,J)
   !             if(nsd==3) GRADPZ = GRADPZ + shg(3,J,L)*PL(dof,J)
   !          end do 
             
   !          do j=1,nenFratura
   !             noGlobal=conecNodaisElem(j,nel)
               
   !             dj(:) = shg(:,j,l)*wDet
   !             indexJ = (j-1)*estrutSistEqFD_%ndof + 1    

   !             if(estrutSistEqFD_%id(1,noGlobal).ne.0) then
   !                if(transiente.eqv..true.) then
   !                   elresf(indexJ)=elresf(indexJ) - abertura * beta * ((pAtual-pAnterior)/deltaT) * dj(nrowsh)
   !                endif      
   !                elresf(indexJ)=elresf(indexJ) - abertura*(dj(1)*condutividade*GRADPX + dj(2)*condutividade*GRADPY)
   !                if(nsd==3) elresf(indexJ)=elresf(indexJ) - abertura*dj(3)*condutividade*GRADPZ
   !             endif
   !          end do
   !       end do

   !       ! computation of Dirichlet b.c. contribution
   !       call ztest(fl,nee,zerodl)
   !       call addrhsN (estrutSistEqFD_%brhs, elresf, estrutSistEqFD_%lmFratura, nee, nel)

   !       deallocate(layDh)
   !       deallocate(layKN)
   !       deallocate(layKT)
   !       deallocate(layBeta)
   !    end do

   ! end subroutine calcCoefSistAlgFluxoDirichlet_DesconBRHS   

   !=======================================================================
   subroutine calcCoefSistAlgFluxoDirichlet_DesconALHS( estrutSistEqFD_)
      
      use mGlobaisEscalares, only: passoTempo, nrowshFratNaBorda, npintFraturasNaBorda
      use mGlobaisArranjos,  only: c, matDescontNaBorda, tipoBordaElemDescontNaBorda, bordaPrescrita
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shgq, shlq, shg1D_points, shl1D_intPoints
      use mfuncoesDeForma,   only: shgCodimOneSurface_points, shgSurface_points
      use mMalha,            only: local, x, nsd
      use mMalha,            only: numelDescontNaBorda, conecNodaisDesconNaBorda, nenFraturasNaBorda
      use mMalha,            only: nenFraturasNaBorda         

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre
      
      implicit none

      type (estruturasArmazenamentoSistemaEq), intent(inout) :: estrutSistEqFD_
      
      real*8 :: xl(nsd,nenFraturasNaBorda)
      real*8 :: shg(nrowshFratNaBorda,nenFraturasNaBorda,npintFraturasNaBorda), shl(nsd-1,nenFraturasNaBorda,npintFraturasNaBorda)
      real*8 :: det(npintFraturasNaBorda), w(npintFraturasNaBorda)
      
      integer*4 :: nee
      real*8   :: eleffm(estrutSistEqFD_%ndof*nenFraturasNaBorda,estrutSistEqFD_%ndof*nenFraturasNaBorda)                
      integer*4:: nel, m, l, i, j

      real*8  :: wDet,  abertura
      logical :: diag,quad,lsym
      real*8  :: di(nrowshFratNaBorda), dj(nrowshFratNaBorda)
      integer*4 :: indexI, indexJ
      integer :: noGlobal

      nee  = nenFraturasNaBorda*estrutSistEqFD_%ndof
      diag = .false.

      w=0.0
      shl=0.0
        
      if(nenFraturasNaBorda==2) then
         call shl1D_intPoints(shl, w, npintFraturasNaBorda, nenFraturasNaBorda)
      else if(nenFraturasNaBorda==3) then
         call shlt(shl,w,npintFraturasNaBorda,nenFraturasNaBorda)
      else if(nenFraturasNaBorda==4) then
         call shlq(shl,w,npintFraturasNaBorda,nenFraturasNaBorda)
      end if                  

      do nel=1, numelDescontNaBorda
      
         if ( bordaPrescrita(tipoBordaElemDescontNaBorda(nel)) ) then
            eleffm=0.0
         
            ! LOCALIZE COORDINATes and Dirichlet b.c.
            call local(conecNodaisDesconNaBorda(:,nel), x, xl, nenFraturasNaBorda, nsd, nsd)            
         
            quad = .true. 
            call shg1D_Points(xl, det, shl, shg, nenFraturasNaBorda, npintFraturasNaBorda, nsd)            
         
            m = matDescontNaBorda(nel)
            abertura = c(2,m)
         
            !.... loop on integration points
            do l=1,npintFraturasNaBorda
               wDet = w(l)*det(l)
            
               do j=1,nenFraturasNaBorda
                  dj(:) = shg(:,j,l)*wDet
                  indexJ = (j-1)*estrutSistEqFD_%ndof+1
               
                  noGlobal=conecNodaisDesconNaBorda(j,nel)
                  do i=1,nenFraturasNaBorda
                     di(:) = shg(:,i,l)
                     indexI = (i-1)*estrutSistEqFD_%ndof+1
                     ! if(transiente) then
                        ! eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ abertura*di(nenFraturasNaBorda)*dj(nenFraturasNaBorda)
                     ! else
                     eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ abertura*di(nrowshFratNaBorda)*dj(nrowshFratNaBorda)
                        ! eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ abertura*di(nrowshFratNaBorda)*dj(nrowshFratNaBorda)
                        ! eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ abertura*di(nenFraturasNaBorda)*dj(nenFraturasNaBorda)
                     ! endif            
                  end do ! i                      
               enddo ! j                   
            end do ! ponto integracao
            
            !.... assemble element stifness matrix and force array into global
            !       left-hand-side matrix and right-hand side vector
            lsym=.true.
         
            if (passoTempo==1) then
               if (estrutSistEqFD_%optSolver=='skyline')   then
                  call addlhsN    (estrutSistEqFD_%alhs, eleffm, estrutSistEqFD_%lmFraturaNaBorda,&
                     estrutSistEqFD_%idiag, nee, nel, diag, lsym)
               endif
               if (estrutSistEqFD_%optSolver=='pardiso') then
                  call addlhsCSR  (estrutSistEqFD_,eleffm, estrutSistEqFD_%lmFraturaNaBorda, nee, nel)
               endif
               if (estrutSistEqFD_%optSolver=='hypre')   then
                  ! call addnslHYPRE(estrutSistEqFD_%A_HYPRE, eleffm, estrutSistEqFD_%lmFraturaNaBorda, nee, nel)
               endif
            endif    
         endif
      end do

   end subroutine calcCoefSistAlgFluxoDirichlet_DesconALHS

   !=======================================================================
   subroutine calcCoefSistAlgFluxoDirichlet_DesconALHS_2D(estrutSistEqFD_)

      use mGlobaisEscalares, only: passoTempo
      use mGlobaisArranjos,  only: c,  matFratura, tiposElemBordas, bordaPrescrita
      use mMalha,            only: local
      use mMalha,            only: numelFratura, numelBordas, conecNodaisFratura
      use mMalha,            only: conecNodaisBordas, nenFratura      

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre
      use mMLayer,           only: LayerMaterial

      implicit none

      type (estruturasArmazenamentoSistemaEq), intent(inout) :: estrutSistEqFD_
      integer*4:: noGlobalF, noGlobalB, nelF, nelB, m, j, k
      real*8  :: abertura

      real*8 :: eleffm(1,1)
      integer*4, pointer :: lm(:,:,:)
      logical :: interceptou

      allocate(lm(1,1,1))
          
      do nelF=1, numelFratura
         m = matFratura(nelF)
         abertura=0.
         do j=1,LayerMaterial(m,1)
            abertura = abertura+c(2,LayerMaterial(m,j+1))
         enddo
         ! print*, "elm", nelF
         ! print*, "mat", m
         ! abertura = c(2,m)
         do j=1, nenFratura
            interceptou = .false.
            noGlobalF = conecNodaisFratura(j, nelF)
            do nelB =1, numelBordas
               if ( (bordaPrescrita(tiposElemBordas(nelB))) .and. (interceptou .eqv. .false.) ) then
                  do k=1, nenFratura
                     noGlobalB = conecNodaisBordas(k, nelB)
                     if ( noGlobalF == noGlobalB ) then                           
                        interceptou = .true.
                        lm(1,1,1) = estrutSistEqFD_%id(1,noGlobalF)
                        eleffm(1,1) = abertura
                        if (passoTempo==1) then
                           if (estrutSistEqFD_%optSolver=='skyline')   then
                             call addlhsN (estrutSistEqFD_%alhs, eleffm, lm, estrutSistEqFD_%idiag, 1, 1, .false., .true.)
                           endif
                           if (estrutSistEqFD_%optSolver=='pardiso') then
                             call addlhsCSR   (estrutSistEqFD_, eleffm, lm, 1, 1)
                           endif
                           if (estrutSistEqFD_%optSolver=='hypre')   then
                              ! call addnslHYPRE(estrutSistEqFD_%A_HYPRE, eleffm, lm, 1, 1)
                           endif            
                        endif    
                     end if
                  end do   
               end if                    
            end do
         end do
      end do

      ! Original.
      ! do nelF=1, numelFratura
      !    m = matFratura(nelF)
      !    abertura=0.
      !    do j=1,LayerMaterial(m,1)
      !       abertura = abertura+c(2,LayerMaterial(m,j+1))
      !    enddo
      !    ! print*, "elm", nelF
      !    ! print*, "mat", m
      !    ! abertura = c(2,m)
      !    do j=1, nenFratura
      !       interceptou = .false.
      !       noGlobalF = conecNodaisFratura(j, nelF)
      !       do nelB =1, numelBordas
      !          if ( (bordaPrescrita(tiposElemBordas(nelB))) .and. (interceptou .eqv. .false.) ) then
      !             do k=1, nenFratura
      !                noGlobalB = conecNodaisBordas(k, nelB)
      !                if ( noGlobalF == noGlobalB ) then                           
      !                   interceptou = .true.
      !                   lm(1,1,1) = estrutSistEqFD_%id(1,noGlobalF)
      !                   eleffm(1,1) = abertura
      !                   if (passoTempo==1) then
      !                      if (estrutSistEqFD_%optSolver=='skyline')   then
      !                        call addlhsN (estrutSistEqFD_%alhs, eleffm, lm, estrutSistEqFD_%idiag, 1, 1, .false., .true.)
      !                      endif
      !                      if (estrutSistEqFD_%optSolver=='pardiso') then
      !                        call addlhsCSR   (estrutSistEqFD_, eleffm, lm, 1, 1)
      !                      endif
      !                      if (estrutSistEqFD_%optSolver=='hypre')   then
      !                         ! call addnslHYPRE(estrutSistEqFD_%A_HYPRE, eleffm, lm, 1, 1)
      !                      endif            
      !                   endif    
      !                end if
      !             end do   
      !          end if                    
      !       end do
      !    end do
      ! end do

   end subroutine calcCoefSistAlgFluxoDirichlet_DesconALHS_2D


   !=======================================================================
   subroutine calcularVazaoNew(estrutSistEqP, estrutSistEqFD, Lx, Ly, Lz, passoTempo)
      use mGlobaisEscalares, only: npintFratura, nrowsh, tempo, deltaT, pi, viscosidade, transiente, fazerWellIndex
      use mGlobaisArranjos,  only: tiposElemBordas, aberturasBorda, vazao
      use mfuncoesDeForma,   only: shlt, shlq, shgCodimOneSurface_points, shl1D_intPoints, shlqen
      use mMalha,            only: local, x, nenFratura, nsd, numelBordas, conecNodaisBordas, numelFratura
      use mMalha,            only: conectividadesInterseccaoFLuxCons, numnp

      type (estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqP, estrutSistEqFD
      integer*4 :: passoTempo
      
      real*8    :: xl(nsd,nenFratura)
      real*8    :: vl(estrutSistEqFD%ndof,nenFratura)
      real*8    :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)
      real*8    :: det(npintFratura), w(npintFratura)
      integer*4 :: nelBordas, tipoBorda, l, i, teste
      real*8    :: vn_ptoInt
      real*8    :: pressaoDir, pressaoEsq, pressaoInf, pressaoSup, pressaoFrente, pressaoAtras, pressaoPoco
      real*8    :: Lx, Ly, Lz, distanciaEntrePontos(3), distancia, areaFratura
      real*8    :: prodInst, vel_media
      real*8    :: area(14), vazaofratura(14)
      real*8    :: vm_x_dir_Fratura, vm_x_esq_Fratura, vm_y_atras_Fratura, vm_y_frente_Fratura, vm_z_sup_Fratura, vm_z_inf_Fratura

      w=0.0
      shl=0.0

      fazerWellIndex = .false.

      if(nenFratura==2) then
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then
         call shlq(shl,w,npintFratura,nenFratura)
         CALL shlqen(SHL,nenFratura)
      end if

      area=0.d0
      vazao=0.d0
      vazaoFratura=0.d0  
      
      do nelBordas=1,numelBordas
         ! if ( (nelBordas /= 25) .AND. (nelBordas /= 48) .AND. (nelBordas /= 73) .AND. (nelBordas /=96) ) then
         ! print*, nelBordas, tiposElemBordas(nelBordas), estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         
         call local(conecNodaisBordas(:,nelBordas), x,               xl,nenFratura, nsd, nsd) 
         CALL local(conecNodaisBordas(:,nelBordas), estrutSistEqFD%u,vl,nenFratura, estrutSistEqFD%ndof, estrutSistEqFD%ndof)          
         
         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)
        
         !... Identificacao do tipo de borda
         ! 1: Frente     (Y)
         ! 2: Direita    (X)
         ! 3: Atrás      (Y)
         ! 4: Esquerda   (X)
         ! 5: Inferior   (Z)
         ! 6: Superior   (Z)
         tipoBorda = tiposElemBordas(nelBordas)
         if(tipoBorda==1) then
            pressaoFrente=estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==2) then
            pressaoDir=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==3) then
            pressaoAtras= estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==4) then
            pressaoEsq=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==5) then
            pressaoInf=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         elseif(tipoBorda==6) then
            pressaoSup=   estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         endif
         
         if(fazerWellIndex) then
            if(tipoBorda==7) pressaoPoco=estrutSistEqP%u(1, conecNodaisBordas(1,nelBordas))
         endif
        
         area(tipoBorda) = area(tipoBorda) + dot_product(w(:),det(:))

         !.... loop on integration points (da borda!!)
         do l=1,npintFratura
            !velocidade normal avaliada no ponto de integracao
            vn_ptoInt = dot_product(vl(1,:), shg(nrowsh,:,l))  
            vazao(tipoBorda) = vazao(tipoBorda) + vn_ptoInt*w(l)*det(l)
         end do

         ! print*,"BordT=",tipoBorda," vel=",vl(1,:)

         !! Patricia - 13/03/2019
         !! Para o caso 3D, este cálculo da vazao na área da fratura faz muuuuito pouca diferença na vazao total
         !! Para o 2D, não testei
         ! print*,'nsd->',nsd
         if (numelFratura>0) then 
            if(nsd==2) then
               do i=1, nenFratura                  
               ! XXX Repensar este caso de forma semelhante ao nsd==3, pq neste caso não tenho velocidadebordasfratura calculado               
                  if(aberturasborda(i,nelbordas)>1.d-10) then
                     area(tipoborda)=area(tipoborda)+aberturasborda(i,nelbordas)/2.d0
                     ! print*, velocidadebordasfratura(i,nelbordas)
                     vazao(tipoborda)=vazao(tipoborda)+velocidadebordasfratura(i,nelbordas)*aberturasborda(i,nelbordas)/2.d0               
                     vazaoFratura(tipoborda)=vazaoFratura(tipoborda)+velocidadebordasfratura(i,nelbordas)*aberturasborda(i,nelbordas)/2.d0               
                  end if
               end do        
            
            else
               do i=1,nenFratura
                  if((conectividadesInterseccaoFLuxCons(i,1,nelBordas).ne.0)) then
                     ! if ( tipoBorda == 5) write(1000,*) "        ", nelBordas, vazao(tipoBorda)                  
                     distanciaEntrePontos(:)=(x(:,conectividadesInterseccaoFLuxCons(i,1,nelBordas)) &
                           - x(:,conectividadesInterseccaoFLuxCons(i,2,nelBordas)))
                     distancia=sqrt(dot_product(distanciaEntrePontos(:),distanciaEntrePontos(:)))

                     areaFratura=(aberturasBorda(i,nelBordas)/2.d0)*distancia
                     
                     if (nenFratura==3) then 
                        if ( i .eq. 1 ) vel_media = (vl(1,1)+vl(1,2))/2.d0
                        if ( i .eq. 2 ) vel_media = (vl(1,2)+vl(1,3))/2.d0                         
                        if ( i .eq. 3 ) vel_media = (vl(1,3)+vl(1,1))/2.d0
                     else if (nenFratura==4) then
                        if ( i .eq. 1 ) vel_media = (vl(1,1)+vl(1,2))/2.d0
                        if ( i .eq. 2 ) vel_media = (vl(1,2)+vl(1,3))/2.d0
                        if ( i .eq. 3 ) vel_media = (vl(1,3)+vl(1,4))/2.d0
                        if ( i .eq. 4 ) vel_media = (vl(1,4)+vl(1,1))/2.d0
                     endif    
                     area(tipoBorda)         =area(tipoBorda)        +areaFratura
                     vazao(tipoBorda)        =vazao(tipoBorda)       +vel_media*areaFratura
                     vazaoFratura(tipoborda) =vazaoFratura(tipoborda)+vel_media*areaFratura
                     ! vazao(tipoBorda)        =vazao(tipoBorda)       +vl(1,i)*areaFratura
                     ! vazaoFratura(tipoborda) =vazaoFratura(tipoborda)+vel_media*areaFratura
                     ! if ( tipoBorda == 5) write(1000,*) "maisFrat", nelBordas, vazao(tipoBorda)
                  endif
               end do   
            endif
         end if
      end do
      
      if ( transiente .EQV. .true.) then 
            do i=1,nsd*2        
               prodInst = (vazaoAnt(i) + vazao(i))/2 * deltaT
               ! if (i==7) print*, "vazaoAnt(i), vazao(i), deltaT", vazaoAnt(i), vazao(i), deltaT
               vazaoAnt(i)  = vazao(i)
               prodAcum(i) = prodAcum(i) + prodInst
               write(500+i,'(1i10, 5es20.8)') passoTempo, tempo, prodInst, prodAcum(i), vazao(i), vazao(i)/area(i)
               print*, 500+i, passoTempo, tempo, prodInst, prodAcum(i), vazao(i), vazaoAnt(i), vazao(i)/area(i)
            end do
      end if
      
      write(*,'(a, 14e20.10)') "Areas", area
      
      write(*,'(a,1es20.8)') "Vazao na borda da frente   = ", vazao(1)
      write(*,'(a,1es20.8)') "Vazao na borda da direita  = ", vazao(2)
      write(*,'(a,1es20.8)') "Vazao na borda de trás     = ", vazao(3)
      write(*,'(a,1es20.8)') "Vazao na borda da esquerda = ", vazao(4)
      if(nsd==3) then
         write(*,'(a,1es20.8)') "Vazao na borda inferior    = ",vazao(5)
         write(*,'(a,1es20.8)') "Vazao na borda superior    = ",vazao(6)
      endif
      if(fazerWellIndex) write(*,'(a,1es20.8)') "Vazao no poco    = ",vazao(7)

      write(*,'(a,1es20.8)') "Vazao fratura na borda da frente   = ", abs(vazaoFratura(1))
      write(*,'(a,1es20.8)') "Vazao fratura na borda da direita  = ", abs(vazaoFratura(2))
      write(*,'(a,1es20.8)') "Vazao fratura na borda de trás     = ", abs(vazaoFratura(3))
      write(*,'(a,1es20.8)') "Vazao fratura na borda da esquerda = ", abs(vazaoFratura(4))
      if(nsd==3) then
         write(*,'(a,1es20.8)') "Vazao fratura na borda inferior    = ",abs(vazaoFratura(5))
         write(*,'(a,1es20.8)') "Vazao fratura na borda superior    = ",abs(vazaoFratura(6))
      endif      
    
      Lx=(maxval(x(1,:)) - minval(x(1,:)))
      Ly=(maxval(x(2,:)) - minval(x(2,:)))
      if(nsd==3) Lz=(maxval(x(3,:)) - minval(x(3,:)))
      
      !! GradP = (p2-p1)/L   
      gradPx=(pressaoEsq-pressaoDir)/Lx
      gradPy=(pressaoFrente-pressaoAtras)/Ly
      print*, pressaoFrente, pressaoAtras, Ly
      write(*,'(a,1es20.8)') "gradPx:", gradPx
      write(*,'(a,1es20.8)') "gradPy:", gradPy
      
      if(nsd==3) then
         gradPz=(pressaoInf-pressaoSup)/Lz
         print*, "pressaoInf", pressaoInf
         print*, "pressaoSup", pressaoSup
         write(*,'(a,1es20.8)') "gradPz:", gradPz
      endif

      vm_x_dir    = vazao(2)/area(2)
      vm_x_esq    = vazao(4)/area(4)      
      vm_y_atras  = vazao(3)/area(3)
      vm_y_frente = vazao(1)/area(1)
      if(nsd==3) then
         vm_z_sup = vazao(6)/area(6)
         vm_z_inf = vazao(5)/area(5)    
      endif
      
      vm_x_dir_Fratura    = vazaoFratura(2)/area(2)
      vm_x_esq_Fratura    = vazaoFratura(4)/area(4)      
      vm_y_atras_Fratura  = vazaoFratura(3)/area(3)
      vm_y_frente_Fratura = vazaoFratura(1)/area(1)
      if(nsd==3) then
         vm_z_sup_Fratura = vazaoFratura(6)/area(6)
         vm_z_inf_Fratura = vazaoFratura(5)/area(5)    
      endif

      write(*,'(a)') "Velocidade media Total (vazaoTotalBordo/areaTotalBordo)"
      write(*,'(a,es18.10)') "Velocidade media na borda da frente = ",   vm_y_frente
      write(*,'(a,es18.10)') "Velocidade media na borda da direita  = ", vm_x_dir
      write(*,'(a,es18.10)') "Velocidade media na borda de trás = ",     vm_y_atras
      write(*,'(a,es18.10)') "Velocidade media na borda da esquerda = ", vm_x_esq
      if(nsd==3) then      
         write(*,'(a,es18.10)') "Velocidade media na borda inferior = ",  vm_z_inf
         write(*,'(a,es18.10)') "Velocidade media na borda superior  = ", vm_z_sup
      endif      
      
      write(*,'(a)') "Velocidade media na fratura (vazaoNasFraturasDoBordo/areaTotalBordo)"
      write(*,'(a,es18.10)') "Velocidade media na borda da frente = ",   vm_y_frente_Fratura
      write(*,'(a,es18.10)') "Velocidade media na borda da direita  = ", vm_x_dir_Fratura
      write(*,'(a,es18.10)') "Velocidade media na borda de trás = ",     vm_y_atras_Fratura
      write(*,'(a,es18.10)') "Velocidade media na borda da esquerda = ", vm_x_esq_Fratura
      if(nsd==3) then      
         write(*,'(a,es18.10)') "Velocidade media na borda inferior = ",  vm_z_inf_Fratura
         write(*,'(a,es18.10)') "Velocidade media na borda superior  = ", vm_z_sup_Fratura
      endif      
      write(*,*)
   end subroutine calcularVazaoNew

   !=======================================================================


end module
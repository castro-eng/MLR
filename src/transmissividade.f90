!=================================================================================
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
  
!>        The function of this module is the control of the whole
!>        structure to determine the transmissibility.

!>         Module implemented by: 
!>         Eduardo Castro,        ecastro@lncc.br
!>
!>         LNCC/MCT
!>         Petropolis, 03.2020
!=================================================================================
!> @var trans_nelm - Number of element that compose the surface
!>                     to evaluate the transmissibility 
!> @var trans_lm - An array to storage the conectivity of the
!>                   elements that compose the surface for 
!>                   transmissibility.
!> @var trans_limite - variable to define the plane 
!> @var errDistanciaPlano - 
!> @var Material_Ignore - A list of material to be ignore for
!>                          the evaluate of the transmissibility.
!>  Eduardo Castro (ecastro@lncc.br)
module mTransmissividade

   implicit none

   !> Structure to storage a list of integer number.
   !> @param nelm - Number of itens in this list.
   !> @param list - The list of integer.
   type int_list
      integer :: nelm
      integer, POINTER :: list(:)=>null()
   end type

   ! Lista de variáveis para o modulo
   integer :: trans_nelm
   integer*4, pointer :: trans_lm(:,:)=>null()
   DOUBLE PRECISION :: trans_limite
   DOUBLE PRECISION, PARAMETER :: errDistanciaPlano=0.1d-5
   TYPE (int_list) :: Material_Ignore

   contains

   ! --------------------------------------------------------------------------------
   !> Routine to evaluate the transmissibility
   !> @param T - equivalent transmissibility between cells.
   !> @param vaz - mass flow rate through the surface 
   subroutine transmissivity(T,vaz)

      use mMalha, only: nsd
      use mGlobaisEscalares, only: viscosidade
      implicit none

      DOUBLE PRECISION, INTENT(OUT) :: T
      DOUBLE PRECISION, INTENT(IN) :: vaz(2)
      DOUBLE PRECISION :: P1, P2, dl, q
      DOUBLE PRECISION, DIMENSION(nsd) ::  x1, x2, xl, vnormal

      integer :: i
      double precision, dimension(8) :: out
     
      ! 1- Determinação da vazão
      CALL vazaoTrans(q,xl,vnormal)
      ! WRITE(*,'("normal(",e12.5,",",e12.5,")")') vnormal(1), vnormal(2)

      ! 2- Determinação da pressão média e centro de massa sobre ambas regiões
      CALL presaoTrans(P1, P2, x1, x2, xl, vnormal)

      ! 3- Determinação da distância entre as pressões médias
      dl=distancia(X1,x2,nsd)

      ! 4- Calculo da Transmissividade
      write(*,'(a)') "--------Transmissividade--------"
      write(*,'(a,1es20.8)') "                        pressaoA ->", P1
      write(*,'(a,1es20.8)') "                        pressaoB ->", P2
      T = (q*viscosidade)/(P1-P2)
      write(*,'(a,1es20.8)') "                           vazão ->", q
      write(*,'(a,1es20.8)') "                Transmissividade ->", T
      T = (vaz(1)*viscosidade)/(P1-P2)
      write(*,'(a,1es20.8)') "           vazão consistente (D) ->", vaz(1)
      write(*,'(a,1es20.8)') "Transmissividade consistente (D) ->", T
      T = (vaz(2)*viscosidade)/(P1-P2)
      write(*,'(a,1es20.8)') "           vazão consistente (E) ->", vaz(2)
      write(*,'(a,1es20.8)') "Transmissividade consistente (E) ->", T

      ! T = (vaz(1)*viscosidade)/(1e5)
      ! write(*,'(a,1es20.8)') "           vazão consistente (D) ->", vaz(1)
      ! write(*,'(a,1es20.8)') "  Permeabilidade consistente (D) ->", T
      ! T = (vaz(2)*viscosidade)/(1e5)
      ! write(*,'(a,1es20.8)') "           vazão consistente (E) ->", vaz(2)
      ! write(*,'(a,1es20.8)') "  Permeabilidade consistente (E) ->", T
      
      out(1)=vaz(1)
      out(2)=vaz(2)
      out(3)=q
      out(4)=(vaz(1)*viscosidade)/(P1-P2)
      out(5)=(vaz(2)*viscosidade)/(P1-P2)
      out(6)=(q*viscosidade)/(P1-P2)
      out(7)=P1
      out(8)=P2
      open(unit=756, action='write', position='append', file= './out/res.inc')
      write(756,'(8es20.8)')  (out(i),i=1,8)
      close(756)



      if (Material_Ignore%nelm.gt.0) deallocate(Material_Ignore%list)

   end subroutine

   ! --------------------------------------------------------------------------------
   !> Routine for determining flow within a rocky massive
   !> @param q - flow over established surface.
   !> @param xl - center of mass of the flow
   !> @param vnormal - normal vector of the established surface
   subroutine vazaoTrans(q, xl, vnormal)
      !     nsd -> tamanho do espaço
      !   numnp -> quantidade de nós
      !       x -> posição espacial dos nós
      ! estrutSistEqF -> estrutura com os valores de velocidade

      use mMalha,    only: nsd, numnp, x
      use mFluxo,    only: estrutSistEqF
      implicit none
   
      DOUBLE PRECISION, INTENT(OUT) :: q, xl(nsd), vnormal(nsd)
      INTEGER :: i
      DOUBLE PRECISION :: xtemp(nsd)

      ! 1- Calculo da vazão
      call vazaoLinha(q, vnormal, x, trans_lm, estrutSistEqF%u, trans_nelm, nsd, numnp)

      ! 2- Calculo do centro de massa
      xl(:)=0d0
      do i=1,trans_nelm
         call centroMassa(xtemp,x(:,trans_lm(i,:)), nsd, 2)
         xl=xl+xtemp
      end do      
      do i=1,nsd
         xl(i)=xl(i)/DBLE(trans_nelm)
      enddo

   end subroutine
   
   ! --------------------------------------------------------------------------------
   !> Routine for determining flow within a rocky massive
   !> @param       q - array with the flow rate value.
   !> @param vnormal - normal vector of the established surface.
   !> @param   coord - nodal coordinate array of the model.
   !> @param      LM - array of connectivities referring
   !>                  to the flow calculation surface.
   !> @param     vel - array with the surface velocities for
   !>                  flow calculation.
   !> @param    nelm - number of elements on the surface for
   !>                  flow calculation.
   !> @param     nsd - space dimensions.
   !> @param   numnp - number of nodal points.
   subroutine vazaoLinha(q, vnormal, coord, LM, vel, nelm, nsd, numnp)

      implicit none

      INTEGER, INTENT(IN):: nelm, nsd, numnp
      INTEGER, INTENT(IN):: LM(nelm,2)
      DOUBLE PRECISION, INTENT(IN) :: coord(nsd,numnp), vel(nsd,numnp)
      DOUBLE PRECISION, INTENT(OUT) :: q, vnormal(nsd)
      INTEGER :: i
      DOUBLE PRECISION :: qint(nsd)

      ! Calculo da vazão sobre uma linha
      q=0d0
      do i=1, nelm
         call vecNormal(vnormal, coord(:,LM(i,1)), coord(:,LM(i,2)) )
         qint = distancia(coord(:,LM(i,1)),coord(:,LM(i,2)),nsd) * (vel(:,LM(i,1)) + vel(:,LM(i,2)))/2
         q = q + dot_product(qint,vnormal)
         ! print*, "vazão do elemento(",i,")=", dot_product(qint,vnormal)
      end do

      RETURN
   end subroutine
   
   ! --------------------------------------------------------------------------------
   !> Routine for determining the center of mass
   !> @param      xc - center of mass.
   !> @param   coord - nodal coordinate array.
   !> @param     nsd - space dimensions.
   !> @param     nen - number of nodal points for each element.
   subroutine centroMassa(xc, coord, nsd, nen)

      implicit none

      INTEGER, INTENT(IN):: nsd, nen
      DOUBLE PRECISION, INTENT(IN) :: coord(nsd,nen)
      DOUBLE PRECISION, INTENT(OUT) :: xc(nsd)

      INTEGER :: i,j
      DOUBLE PRECISION :: val

      xc(:)=0d0
      do i=1,nsd
         val=0d0
         do j=1,nen
            val=val+coord(i,j)
         end do
         xc(i)=xc(i)+val/dble(nen)
      end do

   end subroutine

   ! --------------------------------------------------------------------------------
   !> Routine for determining the average pressures and center
   !> of mass around a surface.
   !> @param   PresA - average pressure in domain A.
   !> @param   PresB - average pressure in domain B.
   !> @param      xA - center of mass of domain A.
   !> @param      xB - center of mass of domain B.
   !> @param      xI - point belonging to the transmissivity plane.
   !> @param vnormal - normal vector defining the transmissivity plane.
   subroutine presaoTrans(PresA, PresB, xA, xB, xI, vnormal)
      use mMalha,    only: nsd, numelFratura

      implicit none

      DOUBLE PRECISION, DIMENSION(nsd), INTENT( IN) :: xI, vnormal
      DOUBLE PRECISION, DIMENSION(nsd), INTENT(OUT) :: xA, xB
      DOUBLE PRECISION, INTENT(OUT) :: PresA, PresB

      INTEGER :: i
      DOUBLE PRECISION :: CMDom(nsd,2), CMRed(nsd,2)
      DOUBLE PRECISION :: PresDom(2), PresRed(2)
      DOUBLE PRECISION :: AreaDom(2), AreaRed(2)

      ! Determinação da pressão, área e centro de massa relativos aos elementos não reduzidos
      call AcumDominio(PresDom,AreaDom,CMDom,xI,vnormal)
      
      ! ------------------------------------------------------------------
      ! Determinação da pressão, área e centro de massa relativos aos elementos reduzidos
      ! if (numelFratura > 0) then
      !    call AcumDominioReduzido(PresRed,AreaRed,CMRed,xI,vnormal)
      ! end if
      
      ! WRITE(*,'(2("PresDom",e12.5,5x))') (PresDom(i),i=1,2)
      ! WRITE(*,'(2("PresRed",e12.5,5x))') (PresRed(i),i=1,2)
      ! WRITE(*,'(2("AreaDom",e12.5,5x))') (AreaDom(i),i=1,2)
      ! WRITE(*,'(2("AreaRed",e12.5,5x))') (AreaRed(i),i=1,2)

      ! Determinação do centro de massa 
      do i=1,nsd
         xA(i)=CMDom(i,1)
         xB(i)=CMDom(i,2)
         ! if(AreaRed(1).ne.0d0) xA(i)=(CMDom(i,1)*AreaDom(1)+CMRed(i,1)*AreaRed(1))/(AreaDom(1)+AreaRed(1))
         ! if(AreaRed(2).ne.0d0) xB(i)=(CMDom(i,2)*AreaDom(2)+CMRed(i,2)*AreaRed(2))/(AreaDom(2)+AreaRed(2))
      enddo

      ! Determinação da pressão 
      PresA=PresDom(1)/AreaDom(1)
      PresB=PresDom(2)/AreaDom(2)
      ! if(AreaRed(1).ne.0d0) PresA=(PresDom(1)+PresRed(1))/(AreaDom(1)+AreaRed(1))
      ! if(AreaRed(2).ne.0d0) PresB=(PresDom(2)+PresRed(2))/(AreaDom(2)+AreaRed(2))
   end subroutine

   ! --------------------------------------------------------------------------------
   !> Routine for determining the pressure and accumulated
   !> area in a domain broken by a plane
   !> @param PresAcum - vector with the pressures accumulated by the broken
   !>                   domain.
   !> @param AreaAcum - vector with the areas accumulated by the broken
   !>                   domain.
   !> @param   xsides - vector with the centers of mass of the sub-domains.
   !> @param       xI - point belonging to the transmissivity plane.
   !> @param  vnormal - normal vector defining the transmissivity plane.
   subroutine AcumDominio(PresAcum, AreaAcum, xsides, xI, vnormal)
      use mGlobaisEscalares, only: nrowsh, npint
      use mGlobaisArranjos,  only: mat
      use mMalha,    only: nsd, nen, x, numel
      use mMalha,    only: conecNodaisElem, normal, numelFratura
      use mMalha,    only: local
      use mPotencial,        only: estrutSistEqP
      use mFuncoesDeForma,   only: shlt, shlq, shlt3D, shlq3d
      use mFuncoesDeForma,   only: shlten, shlt, shlqen, shlq
      use mFuncoesDeForma,   only: shlten3D, shlt3D, shlqen3D, shlq3d
      use mFuncoesDeForma,   only: oneshg, shgq, shg3d
      use mFuncoesDeForma,   only: shgCodimOneSurface_points, shl1D_intPoints
      
      implicit none
      
      DOUBLE PRECISION, DIMENSION(nsd), INTENT( IN) :: xI, vnormal
      DOUBLE PRECISION, DIMENSION(2), INTENT(OUT) :: PresAcum, AreaAcum
      DOUBLE PRECISION, DIMENSION(nsd,2), INTENT(OUT) :: xsides
      
      INTEGER :: i, j, nel, noGlobal, side, dof, mat_ign
      DOUBLE PRECISION :: PressTemp, AreaTemp, wDet, dLPlano
      DOUBLE PRECISION :: xtemp(nsd)
      DOUBLE PRECISION :: w(npint), det(npint),  vetorParaNo(nsd)
      DOUBLE PRECISION :: shg(nrowsh,nen,npint), shl(nrowsh,nen,npint)
      DOUBLE PRECISION :: xl(nsd,nen), normalNoLocal(nsd,nen)
      DOUBLE PRECISION :: shlCentroElem(nrowsh,nen,1), wCentroElem(1), xCentroElem(nsd)
      LOGICAL :: flagPressaoMedia, quad
      
      ! Calculo dos pesos e da função de interpolação local
      w=0.0
      shl=0.0
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
         call shlq3d(shlCentroElem,wCentroElem,1,nen)
      end if
      
      ! Loop nos elementos
      PresAcum(:)=0d0
      AreaAcum(:)=0d0
      xsides(:,:)=0d0
      do nel=1,numel
         flagPressaoMedia=.False.
         call centroMassa(xtemp,x(:,conecNodaisElem(:,nel)), nsd, nen)
         dLPlano = distanciaPontoPlano(xtemp, vnormal, xI, nsd)
         
         ! Avaliação se pertence ao domínio A 
         if ((dLPlano.le.(0.e0)) .and. (dLPlano.ge.(-trans_limite))) then
            flagPressaoMedia=.TRUE.
            side=1
         end if
         
         ! Avaliação se pertence ao domínio B
         if ((dLPlano.ge.(0.e0)) .and. (dLPlano.le.(trans_limite))) then 
            flagPressaoMedia=.TRUE.
            side=2
         end if

         ! Lista dos materiais igonrados
         do mat_ign=1,Material_Ignore%nelm
            if (mat(nel).eq.Material_Ignore%list(mat_ign)) then
               flagPressaoMedia=.False.
               exit
            end if
         end do

         ! Caso estar no domínio A ou B, determinar pressão e area do elemento e acumular
         if (flagPressaoMedia) then
            ! LOCALIZE COORDINATes and Dirichlet b.c.
            call local(conecNodaisElem(1,nel), x, xl, nen, nsd, nsd)
            
            if (numelFratura > 0) then
               call local(conecNodaisElem(:,nel), normal, normalNoLocal, nen, nsd, nsd)
               do i=1,nsd
                  xCentroElem(i) = dot_product(xl(i,:), shlCentroElem(nrowsh,:,1))
               end do
            end if
            
            quad = .true.
            if (nen.eq.4.and.conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
            
            if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,1,nel,1)
            if(nen==3) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
            if(nen==4.and.nsd==2) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
            if(nen==4.and.nsd==3) call shg3d (xl,det,shl,shg,npint,nel,nen)
            if(nen==8) call shg3d (xl,det,shl,shg,npint,nel,nen)
            
            PressTemp = 0d0
            AreaTemp  = 0d0
            do i=1,npint
               wDet = w(i)*det(i)
               do j=1,nen
                  dof = 1
                  if(numelFratura>0)then
                     if ((estrutSistEqP%lm(2,J,nel) > 0) .or. (estrutSistEqP%lm(3,J,nel) > 0)) then
                        vetorParaNo(:) = xl(:,J) - xCentroElem(:)                   
                        if (dot_product(normalNoLocal(:,J), vetorParaNo(:)) < 0d0) then
                           dof = 2
                           ! print*, dof
                        end if
                     end if
                  end if
                  AreaTemp = AreaTemp + shg(nrowsh,j,i)*wDet
                  noGlobal=conecNodaisElem(j,nel)
                  PressTemp = PressTemp + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(dof,noGlobal)
               end do
            end do
            ! Atualização do centro de massa do sub-domínio
            do i=1,nsd
               xsides(i,side) = ( xsides(i,side)*AreaAcum(side) + xtemp(i)*AreaTemp)/(AreaAcum(side)+AreaTemp)
            end do
            ! Atualização da pressão no sub-domínio
            PresAcum(side) = PresAcum(side) + PressTemp
            AreaAcum(side) = AreaAcum(side) + AreaTemp
         end if
      end do

   end subroutine

   ! --------------------------------------------------------------------------------
   !> Routine for determining the pressure and accumulated area in a reduced
   !> domain broken by a plane
   !> @param PresAcum - vector with the pressures accumulated by the broken
   !>                   domain.
   !> @param AreaAcum - vector with the areas accumulated by the broken
   !>                   domain.
   !> @param   xsides - vector with the centers of mass of the sub-domains.
   !> @param       xI - point belonging to the transmissivity plane.
   !> @param  vnormal - normal vector defining the transmissivity plane.
   subroutine AcumDominioReduzido(PresAcum, AreaAcum, xsides, xI, vnormal)
      use mGlobaisEscalares, only: nrowsh, npintFratura
      use mGlobaisArranjos,  only: c, matFratura
      use mMalha,    only: nsd, x
      use mMalha,    only: numelFratura, nenFratura, conecNodaisFratura
      use mMalha,    only: local
      use mPotencial,        only: estrutSistEqP
      use mFuncoesDeForma,   only: shlt, shlq, shlt3D, shlq3d
      use mFuncoesDeForma,   only: shlten, shlt, shlqen, shlq
      use mFuncoesDeForma,   only: shlten3D, shlt3D, shlqen3D, shlq3d
      use mFuncoesDeForma,   only: oneshg, shgq, shg3d
      use mFuncoesDeForma,   only: shgCodimOneSurface_points, shl1D_intPoints
      use mMLayer,           only: LayerMaterial
      
      implicit none
      
      DOUBLE PRECISION, DIMENSION(nsd), INTENT( IN) :: xI, vnormal
      DOUBLE PRECISION, DIMENSION(2), INTENT(OUT) :: PresAcum, AreaAcum
      DOUBLE PRECISION, DIMENSION(nsd,2), INTENT(OUT) :: xsides

      INTEGER :: i, j, k, nel, noGlobal, side, ncamadas
      DOUBLE PRECISION :: wDet, dLPlano
      DOUBLE PRECISION :: PressTempA, PressTempB, AreaTempA, AreaTempB
      DOUBLE PRECISION :: xtemp(nsd)
      
      DOUBLE PRECISION :: xl(nsd,nenFratura)
      DOUBLE PRECISION :: w(npintFratura), det(npintFratura)
      DOUBLE PRECISION :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)

      DOUBLE PRECISION :: layDh(estrutSistEqP%ndof)

      LOGICAL :: flagPressaoMedia, flagCamada, quad
      
      ! Calculo dos pesos e da função de interpolação local
      w=0.0
      shl=0.0
      if(nenFratura==2) then
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then
         call shlq(shl,w,npintFratura,nenFratura)
      end if
   
      ! Loop nos elementos reduzidos
      PresAcum(:)=0d0
      AreaAcum(:)=0d0
      xsides(:,:)=0d0
      do nel=1,numelFratura
         flagPressaoMedia=.False.
         call centroMassa(xtemp,x(:,conecNodaisFratura(:,nel)), nsd, nenFratura)
         dLPlano = distanciaPontoPlano(xtemp, vnormal, xI, nsd)
         
         ! Avaliação sobre qual domínio pertence o elemento
         flagCamada=.FALSE.
         if((dLPlano.ge.(-errDistanciaPlano)).and.(dLPlano.le.(errDistanciaPlano))) then
            flagPressaoMedia=.TRUE.
            flagCamada=.TRUE.
         else if ((dLPlano.le.(-errDistanciaPlano)) .and. (dLPlano.ge.(-trans_limite))) then
            flagPressaoMedia=.TRUE.
            side=1
         else if ((dLPlano.ge.(errDistanciaPlano)) .and. (dLPlano.le.(trans_limite))) then
            flagPressaoMedia=.TRUE.
            side=2
         end if
         
         ! Caso o elemento se encontre sobre a região analisada,
         ! determinação das propriedades necessárias para calculo da pressão
         if (flagPressaoMedia) then
            call local(conecNodaisFratura(:,nel), x, xl, nenFratura, nsd, nsd)
            quad = .true.
            if (nenFratura.eq.4.and.conecNodaisFratura(3,nel).eq.conecNodaisFratura(4,nel)) quad = .false.
            call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)
            ncamadas = LayerMaterial(matFratura(nel),1)
            layDh = 0d0
            do k=1,ncamadas
               layDh(k) = c(2,LayerMaterial(matFratura(nel),k+1))
            end do
         end if
         
         ! Situação em que a descontinuidade se encontra coincidente com o plano da transmissividade
         if (flagCamada) then
            PressTempA = 0d0
            PressTempB = 0d0
            AreaTempA  = 0d0
            AreaTempB  = 0d0

            if (MOD(ncamadas,2).eq.1) then

               ! Situação em que o número de camadas é impar
               if (ncamadas.gt.1) then

                  ! Sub-domínio A
                  do k=1,(ncamadas-1)/2
                     do i=1,npintFratura
                        wDet = w(i)*det(i)
                        do j=1,nenFratura
                           AreaTempA = AreaTempA + shg(nrowsh,j,i)*wDet*layDh(k)
                           noGlobal=conecNodaisFratura(j,nel)
                           PressTempA = PressTempA + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(k+2,noGlobal)!/layDh(k)
                        end do
                     end do
                  end do
                  
                  ! Sub-domínio B
                  do k=(ncamadas-1)/2+2,ncamadas
                     do i=1,npintFratura
                        wDet = w(i)*det(i)
                        do j=1,nenFratura
                           AreaTempB = AreaTempB + shg(nrowsh,j,i)*wDet*layDh(k)
                           noGlobal=conecNodaisFratura(j,nel)
                           PressTempB = PressTempB + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(k+2,noGlobal)!/layDh(k)
                        end do
                     end do
                  end do
               end if
               
               k=(ncamadas-1)/2+1
               do i=1,npintFratura
                  wDet = w(i)*det(i)
                  do j=1,nenFratura
                     AreaTempA = AreaTempA + shg(nrowsh,j,i)*wDet*(layDh(k)/2.)
                     AreaTempB = AreaTempB + shg(nrowsh,j,i)*wDet*(layDh(k)/2.)
                     noGlobal=conecNodaisFratura(j,nel)
                     PressTempA = PressTempA + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(k+2,noGlobal)!*2/layDh(k)
                     PressTempB = PressTempB + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(k+2,noGlobal)!*2/layDh(k)
                  end do
               end do
               
               ! Atualização do centro de massa do sub-domínio
               do i=1,nsd
                  xsides(i,1) = ( xsides(i,1)*AreaAcum(1) + xtemp(i)*AreaTempA)/(AreaAcum(1)+AreaTempA)
                  xsides(i,2) = ( xsides(i,2)*AreaAcum(2) + xtemp(i)*AreaTempB)/(AreaAcum(2)+AreaTempB)
               end do
               ! Atualização da pressão no sub-domínio
               PresAcum(1) = PresAcum(1) + PressTempA
               PresAcum(2) = PresAcum(2) + PressTempB
               AreaAcum(1) = AreaAcum(1) + AreaTempA
               AreaAcum(2) = AreaAcum(2) + AreaTempB
   
            else
               ! Situação em que o número de camadas é par
               ! Sub-domínio A
               do k=1,(ncamadas-1)/2
                  do i=1,npintFratura
                     wDet = w(i)*det(i)
                     do j=1,nenFratura
                        AreaTempA = AreaTempA + shg(nrowsh,j,i)*wDet*layDh(k)
                        noGlobal=conecNodaisFratura(j,nel)
                        PressTempA = PressTempA + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(k+2,noGlobal)!/layDh(k)
                     end do
                  end do
               end do
               
               ! Sub-domínio B
               do k=(ncamadas-1)/2+1,ncamadas
                  do i=1,npintFratura
                     wDet = w(i)*det(i)
                     do j=1,nenFratura
                        AreaTempB = AreaTempB + shg(nrowsh,j,i)*wDet*layDh(k)
                        noGlobal=conecNodaisFratura(j,nel)
                        PressTempB = PressTempB + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(k+2,noGlobal)!/layDh(k)
                     end do
                  end do
               end do

               ! Atualização do centro de massa do sub-domínio
               do i=1,nsd
                  xsides(i,1) = ( xsides(i,1)*AreaAcum(1) + xtemp(i)*AreaTempA)/(AreaAcum(1)+AreaTempA)
                  xsides(i,2) = ( xsides(i,2)*AreaAcum(2) + xtemp(i)*AreaTempB)/(AreaAcum(2)+AreaTempB)
               end do
               ! Atualização da pressão no sub-domínio
               PresAcum(1) = PresAcum(1) + PressTempA
               PresAcum(2) = PresAcum(2) + PressTempB
               AreaAcum(1) = AreaAcum(1) + AreaTempA
               AreaAcum(2) = AreaAcum(2) + AreaTempB
            end if
         end if
   
         ! Situação em que a descontinuidade fora do plano da transmissividade
         if (flagPressaoMedia.and.(.not.flagCamada)) then
            PressTempA = 0d0
            AreaTempA  = 0d0
            do k=1,ncamadas
               do i=1,npintFratura
                  wDet = w(i)*det(i)
                  do j=1,nenFratura
                     AreaTempA = AreaTempA + shg(nrowsh,j,i)*wDet   
                     noGlobal=conecNodaisFratura(j,nel)
                     PressTempA = PressTempA + wDet*shg(nrowsh,j,i)*estrutSistEqP%u(k+2,noGlobal)!/layDh(k)
                  end do
               end do
            end do

            ! Atualização do centro de massa do sub-domínio
            do i=1,nsd
               xsides(i,side) = ( xsides(i,side)*AreaAcum(side) + xtemp(i)*AreaTempA)/(AreaAcum(side)+AreaTempA)
            end do
            ! Atualização da pressão no sub-domínio
            PresAcum(side) = PresAcum(side) + PressTempA
            AreaAcum(side) = AreaAcum(side) + AreaTempA
         endif
      end do
   end subroutine

   ! --------------------------------------------------------------------------------
   !> Routine for determining the normal vector to a line defined by two
   !> points considers itself vector perpendicular to the z-plane
   !> @param   vect - normal vector of the line.
   !> @param pointA - point A on the line.
   !> @param pointB - point B on the line.
   subroutine vecNormal(vect, pointA, pointB)

      implicit none

      DOUBLE PRECISION, DIMENSION(2), INTENT( IN) :: pointA, pointB
      DOUBLE PRECISION, DIMENSION(2), INTENT(OUT) :: vect
      DOUBLE PRECISION, DIMENSION(2) :: vecTang

      vecTang = pointB-pointA
      vect(1) =  vecTang(2)
      vect(2) = -vecTang(1)
      vect = vect/(sqrt(dot_product(vect,vect)))
   end subroutine

   ! --------------------------------------------------------------------------
   !> Routine for determining the distance between two points
   !> @param  dL - distance between A and B by L2 norm.
   !> @param  xA - point A.
   !> @param  xB - point B.
   !> @param nsd - space dimension.
   double precision function distancia(xA, xB, nsd) result(dL) 
      INTEGER, INTENT(IN) :: nsd
      DOUBLE PRECISION, INTENT( IN) :: xA(nsd), xB(nsd)
      DOUBLE PRECISION :: diffAeB(nsd)

      diffAeB = xA-xB
      dL = sqrt(dot_product(diffAeB,diffAeB))
   end function distancia

   ! --------------------------------------------------------------------------
   !> Routine for determining the distance between a point and the plane/line 
   !> @param      dL - distance between a point and the plane/line.
   !> @param      Pt - point value.
   !> @param  normal - normal vector of the plane
   !> @param PtPlano - point belonging to the plane.
   !> @param     nsd - space dimension.
   double precision function distanciaPontoPlano(Pt, normal, PtPlano, nsd) result(dL)
      INTEGER, INTENT(IN) :: nsd
      DOUBLE PRECISION, INTENT( IN) :: Pt(nsd), normal(nsd), PtPlano(nsd)
      INTEGER :: i

      dL =0.
      do i=1,nsd
         dL = dL + normal(i)*(Pt(i)-PtPlano(i))
      end do
      dL = dL/(sqrt(dot_product(normal,normal)))
   end function distanciaPontoPlano

end module mTransmissividade
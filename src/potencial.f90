!=================================================================================
!
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Garcia,        bidu@lncc.br [4]
!         Tuane Lopes,          tuane@lncc.br [5]
!
!         LNCC/MCT
!         Petropolis, 12.2015
!=================================================================================

module mPotencial
        
   use mEstruturasDadosSistEq
        
   type (estruturasArmazenamentoSistemaEq) :: estrutSistEqP
        
   contains 

   !==========================================================================================================
   subroutine montarSistEqAlgPotencial(estrutSistEqP_)

      use mMalha,               only: x, conecNodaisElem, numel, nen, nsd, numnp, numelFratura
      use mUtilSistemaEquacoes, only: dirichletConditions, load
      use mGlobaisEscalares,    only: transiente, passoTempo
      use mLeituraEscrita,      only: escreverArquivosSaida_Potencial

      implicit none

      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqP_
        
      real*8 :: t1, t2, t3

      IF(PASSOTEMPO==1) THEN
      
         if (estrutSistEqP_%nlvect.gt.0) call load (estrutSistEqP_%id,estrutSistEqP_%f,estrutSistEqP_%brhs,&
                  estrutSistEqP_%ndof,numnp,estrutSistEqP_%nlvect)

         if (estrutSistEqP_%nlvect.gt.0) call dirichletConditions (estrutSistEqP_%id,estrutSistEqP_%u,estrutSistEqP_%f, &
                  estrutSistEqP_%ndof,numnp,estrutSistEqP_%nlvect)
                                     
         call escreverArquivosSaida_Potencial(estrutSistEqP_, 0, 0.d0)
      ENDIF

      call calcCoefSistAlgPotencial (estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
      if (numelFratura > 0) call calcCoefSistAlgPotencial_Descontinuidade(estrutSistEqP_)
    
   end subroutine montarSistEqAlgPotencial

   !==========================================================================================================
   subroutine calcCoefSistAlgPotencial(estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
         
      use mGlobaisEscalares, only: nrowsh, npint, transiente, PASSOTEMPO, viscosidade
      use mGlobaisArranjos,  only: mat, c, grav, campoPermeabilidade
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq, shlt3D
      use mMalha,            only: local, numelFratura, normal
      use mGlobaisEscalares, only: deltaT, saltoPressao

      use mUtilSistemaEquacoes,  only: kdbc2

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre
         
      implicit none

      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqP_
      real*8,    intent(in)    :: x(nsd, numnp)
      integer*4, intent(in)    :: numnp, numel, nen, nsd
      integer*4, intent(in)    :: conecNodaisElem(nen,numel)
         
      real*8 :: xl(nsd,nen), pl(estrutSistEqP_%ndof,nen), plAnt(estrutSistEqP_%ndof,nen)
      real*8 :: plIterAnt(estrutSistEqP_%ndof,nen), plIterAnt_pontoInt
      real*8 :: shg(nrowsh,nen,npint), shl(nrowsh,nen,npint)
      real*8 :: det(npint), w(npint)
         
      integer*4 :: nee
      real*8    :: elresf(estrutSistEqP_%ndof*nen), eleffm(estrutSistEqP_%ndof*nen,estrutSistEqP_%ndof*nen)
         
      integer*4 :: nel, m, l, i, j, ni, nj, sd_i, sd_j, indexI, indexJ, noGlobal
      real*8  :: permeabilidadeX, permeabilidadeY, permeabilidadeZ
      real*8  :: condutividadeX, condutividadeY, condutividadeZ
      real*8  :: beta, coef_poisson, mod_young
      real*8  :: wDet, gf1, gf2, gf3
      real*8  :: djx, djy, djz, djn, dix, diy, diz, din
      real*8  :: UU, UUP(estrutSistEqP_%ndof), betat
      logical :: diag,zerodl,quad,lsym, direita, esquerda
         
      real*8  :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1), vetorParaNo(nsd), normalNoLocal(nsd,nen)
      real*8  :: di(nrowsh), dj(nrowsh)

      real*8 :: poisson, modYoung, bulkGrain, bulkOil, bulkRocha
      real*8 :: alpha, umPorM, porosidade
         
      nee  = nen*estrutSistEqP_%ndof
      diag = .false.
         
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

      do nel=1,numel

         ! clear stiffness matrix and force array
         eleffm=0.0
         elresf=0.0

         ! LOCALIZE COORDINATes and Dirichlet b.c.
         call local(conecNodaisElem(1,nel), x, xl, nen, nsd, nsd)
         call local(conecNodaisElem(1,nel), estrutSistEqP_%u, pl, nen, estrutSistEqP_%ndof, estrutSistEqP_%ndof)
         if(transiente.eqv..true.) then
            call local(conecNodaisElem(1,nel),estrutSistEqP_%uTempoAnt,plAnt,nen,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
         end if
            
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
         
         ! ....... form stiffness matrix
         !... Properties            
         m = mat(nel)
         porosidade = c(1,m)
         permeabilidadeX = campoPermeabilidade(1,nel) !c(3,m)
         permeabilidadeY = campoPermeabilidade(2,nel) !c(4,m)
         permeabilidadeZ = campoPermeabilidade(3,nel) !c(5,m)
         beta=c(6,m)

         condutividadeX = permeabilidadeX/viscosidade
         condutividadeY = permeabilidadeY/viscosidade
         condutividadeZ = permeabilidadeZ/viscosidade
            
         ! poisson=0.2
         ! modYoung=1.0e+09
         ! bulkGrain=7.5e+09
         ! bulkOil=1.62e+09
            
         ! bulkRocha=modYoung/(3.d0*(1.d0-2.d0*poisson))
         ! alpha=1-(bulkRocha/bulkGrain)
            
         ! umPorM=(porosidade/bulkOil)+((alpha-porosidade)/bulkGrain)
            
         ! beta=umPorM + (alpha*alpha)/bulkRocha

         ! mod_young = c(4,m)  ! c(4,m): Modulo de Young
         ! coef_poisson = c(5,m)  ! c(5,m): Coeficiente de Poisson       
         ! beta = 3d0*(1d0-2d0*coef_poisson)/mod_young

         !.... loop on integration points
         do l=1,npint
            wDet = w(l)*det(l)

            if(transiente.eqv..true.) then
               UUP=0.D0
               do j=1,nen
                  UUP(1)  = UUP(1) + SHG(nrowsh,J,L)*plAnt(1,J)
                  if(saltoPressao) UUP(2)  = UUP(2) + SHG(nrowsh,J,L)*plAnt(2,J)
               enddo
            endif

            do j=1,nen

               noGlobal=conecNodaisElem(j,nel)
               direita=.false.
                  
               dj(:) = shg(:,j,l)*wDet
               indexJ = (j-1)*estrutSistEqP_%ndof + 1
                  
               if ((numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) then
                  if ((estrutSistEqP_%lm(2,j,nel) > 0) .or. (estrutSistEqP_%lm(3,j,nel) > 0)) then
                     vetorParaNo(:) = xl(:,j) - xCentroElem(:)
                     if (dot_product(normalNoLocal(:,j), vetorParaNo(:)) < 0d0) then
                        indexJ = indexJ + 1
                        direita=.true.
                     end if
                  end if
               end if

               ! determinacao do termo transiente (dP/dt)
               if ( transiente .eqv. .true. ) then
                  if(saltoPressao) then
                     if(direita.eqv..true.)then
                        elresf(indexJ)=elresf(indexJ) + beta*UUP(2)*dj(nrowsh)
                     else
                        elresf(indexJ)=elresf(indexJ) + beta*UUP(1)*dj(nrowsh)
                     endif
                  else
                     elresf(indexJ)=elresf(indexJ) + beta*UUP(1)*dj(nrowsh)
                  endif
               end if  

               ! determinacao da parela espacial (dP/dx)
               do i=1,nen
                  di(:) = shg(:,i,l)
                  indexI = (i-1)*estrutSistEqP_%ndof + 1
                     
                  if ( (numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) then 
                     if ((estrutSistEqP_%lm(2,i,nel) > 0) .or. (estrutSistEqP_%lm(3,i,nel) > 0)) then
                        vetorParaNo(:) = xl(:,i) - xCentroElem(:)
                        if (dot_product(normalNoLocal(:,i), vetorParaNo(:)) < 0d0) then
                           indexI = indexI + 1
                        end if
                     end if
                  end if

                  if(transiente.eqv..true.) then
                     eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ beta*di(nrowsh)*dj(nrowsh)
                     eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ deltaT*condutividadeX*(di(1)*dj(1)) 
                     eleffm(indexJ, indexI)=eleffm(indexJ, indexI)+ deltaT*condutividadeY*(di(2)*dj(2)) 
                     if(nsd==3) eleffm(indexJ, indexI) = eleffm(indexJ, indexI) + deltaT*condutividadeZ*(di(3)*dj(3)) 
                  else
                     eleffm(indexJ, indexI) = eleffm(indexJ, indexI) + condutividadeX*(di(1)*dj(1))
                     eleffm(indexJ, indexI) = eleffm(indexJ, indexI) + condutividadeY*(di(2)*dj(2))
                     if(nsd==3) eleffm(indexJ, indexI) = eleffm(indexJ, indexI) + condutividadeZ*(di(3)*dj(3)) 
                  endif
               end do
            end do
         end do
            
         ! computation of Dirichlet b.c. contribution
         call ztest(pl,nee,zerodl)
            
         if(.not.zerodl) then
            call KDBC2(ELEFFM,ELRESF,pl,NEE,estrutSistEqP_%lm,NEL)
         endif
            
         !.... assemble element stifness matrix and force array into global
         !        left-hand-side matrix and right-hand side vector
         lsym=.true.
            
         if (passoTempo==1) then
            if (estrutSistEqP_%optSolver=='skyline')   then
               ! Write(*,*) " call addlhsN  (estrutSistEqP_, eleffm, nel) "
               call addlhsN    (estrutSistEqP_%alhs, eleffm, estrutSistEqP_%lm, estrutSistEqP_%idiag, nee, nel, diag, lsym)
            endif
            if (estrutSistEqP_%optSolver=='pardiso') then
               ! Write(*,*) " call addlhsCSR  (estrutSistEqP_, eleffm, nel) "
               call addlhsCSR  (estrutSistEqP_,eleffm, estrutSistEqP_%lm, nee, nel)
            endif
            if (estrutSistEqP_%optSolver=='hypre')   then
               call addnslHYPRE(A_HYPRE, eleffm, estrutSistEqP_%lm, nee, nel)
            endif
         endif

         call addrhsN (estrutSistEqP_%brhs, elresf, estrutSistEqP_%lm, nee, nel)
      end do

      return

   end subroutine calcCoefSistAlgPotencial

   !==========================================================================================================
   subroutine calcCoefSistAlgPotencial_fraturas(estrutSistEqP_)
      use mGlobaisEscalares, only: transiente
      use mGlobaisEscalares, only: npintFratura, viscosidade, nrowsh
      use mGlobaisArranjos,  only: c, grav, matFratura
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shgq, shlq, shg1D_points, shl1D_intPoints
      use mfuncoesDeForma,   only: shgCodimOneSurface_points, shgSurface_points
      use mMalha,            only: local, x, numnp, nenFratura, nsd, numelFratura, conecNodaisFratura, normal
      use mGlobaisEscalares, only: deltaT, passoTempo, carregamento, carregamento_ref, saltoPressao
      use mMLayer,           only: LayerMaterial

      use mUtilSistemaEquacoes,  only: kdbc,KDBC2

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre

      implicit none

      type (estruturasArmazenamentoSistemaEq), intent(inout) :: estrutSistEqP_

      real*8 :: xl(nsd,nenFratura), xlAux(nsd,nenFratura)
      real*8 :: pl(estrutSistEqP_%ndof,nenFratura), plAnt(estrutSistEqP_%ndof,nenFratura)
      real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)
      real*8 :: det(npintFratura), w(npintFratura)
      real*8 :: normalNoLocal(nsd,nenFratura)

      integer*4 :: nee

      real*8   :: elresf(estrutSistEqP_%ndof*nenFratura), eleffm(estrutSistEqP_%ndof*nenFratura,estrutSistEqP_%ndof*nenFratura)

      integer*4:: nel, m, l, i, j, ni, nj, k

      real*8  :: condutividade, abertura, beta, rigidezNormalInicial
      real*8  :: coef_poisson, mod_young
      real*8  :: wDet, gf1, gf2, gf3, permeabilidade
      real*8  :: djx, djy, djz, djn, dix, diy, diz, din
      logical :: diag,zerodl,quad,lsym
      real*8  :: UU, UUP(estrutSistEqP_%ndof)
      real*8  :: resistividade, di(nrowsh), dj(nrowsh), aux, xi
      integer*4 :: indexI, indexJ

      real*8 :: poisson, modYoung, bulkGrain, bulkOil, bulkRocha
      real*8 :: alpha, umPorM, porosidade
      real*8, dimension(:), allocatable :: layKN, layKT, layDh, layBeta
      real*8 :: coeff(9), coeffInter(4)
      integer*4:: ncamadas, posI1, posI2, posI3, posJ1, posJ2, posJ3

      nee  = nenFratura*estrutSistEqP_%ndof
      diag = .false.

      w=0.0
      shl=0.0
      xi=0.49999d0
      ! xi=3.d0/4.d0
      ! xi=1.d0
        
      if(nenFratura==2) then
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then
         call shlq(shl,w,npintFratura,nenFratura)
      end if
         
      do nel=1, numelFratura

         ! clear stiffness matrix and force array
         eleffm=0.0
         elresf=0.0

         ! LOCALIZE COORDINATes and Dirichlet b.c.
         call local(conecNodaisFratura(:,nel), x, xl, nenFratura, nsd, nsd)
         call local(conecNodaisFratura(:,nel), estrutSistEqP_%u, pl, nenFratura, estrutSistEqP_%ndof, estrutSistEqP_%ndof)
         if(transiente.eqv..true.) then
            call local(conecNodaisFratura(:,nel),estrutSistEqP_%uTempoAnt,plAnt,nenFratura,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
         endif

         quad = .true.
         if (nenFratura.eq.4.and.conecNodaisFratura(3,nel).eq.conecNodaisFratura(4,nel)) quad = .false.

         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)

         !....... form stiffness matrix
         !... Properties
         m = matFratura(nel)
         beta=c(6,m)
         ! Propriedades das camadas
         ncamadas = LayerMaterial(m,1)
         allocate(layDh(ncamadas))
         allocate(layKN(ncamadas))
         allocate(layKT(ncamadas))
         allocate(layBeta(ncamadas))
         do k=1, ncamadas
            layDh(k) = c(2,LayerMaterial(m,k+1))
            layKN(k) = (2.d0*c(3,LayerMaterial(m,k+1)))/(viscosidade*layDh(k))
            layKT(k) = c(4,LayerMaterial(m,k+1))/viscosidade
            layBeta(k) = c(6,LayerMaterial(m,k+1))
         enddo

         !.... loop on integration points
         do l=1,npintFratura
            wDet = w(l)*det(l)

            if(transiente.eqv..true.) then
               UUP=0.D00
               DO J=1,nenFratura
                  UUP(1)  = UUP(1) + SHG(nrowsh,J,L)*PLAnt(1,J) 
                  if(saltoPressao)  then
                     UUP(3)  = UUP(3) + SHG(nrowsh,J,L)*PLAnt(3,J) 
                  endif
               ENDDO
            endif

            do j=1,nenFratura
               dj(:) = shg(:,j,l)*wDet
               indexJ = (j-1)*estrutSistEqP_%ndof

               ! determinacao do termo transiente (dP/dt)
               if(transiente.eqv..true.) then
                  elresf(indexJ+3)=elresf(indexJ+3) + abertura*beta*UUP(3)*dj(nrowsh)
               endif

               ! determinacao da parela espacial (dP/dx)
               do i=1, nenFratura
                  di(:) = shg(:,i,l)
                  indexI = (i-1)*estrutSistEqP_%ndof

                  do k=1, ncamadas
                     ! Determinacao das posicoes da camada
                     posJ1 = indexJ+1+k
                     posJ2 = indexJ+3+k
                     posJ3 = indexJ+2+k
                     posI1 = indexI+1+k
                     posI2 = indexI+3+k
                     posI3 = indexI+2+k
                     
                     coeffInter(1)=1.d0;
                     coeffInter(2)=0.d0;
                     coeffInter(3)=0.d0;
                     coeffInter(4)=1.d0;
                     if(k.eq.1) then
                        posI1=indexI+1
                        posJ1=indexJ+1
                        !    coeffInter(1)=1.d0
                        !    coeffInter(2)=0.d0
                        ! else
                        !    coeffInter(1)=5.d-1
                        !    coeffInter(2)=5.d-1
                        !    ! coeffInter(1)=layKN(k-1)/(layKN(k-1)+layKN(k))
                        !    ! coeffInter(2)=layKN(  k)/(layKN(k-1)+layKN(k))
                     endif
                     if(k.eq.ncamadas) then
                        posI2=indexI+2
                        posJ2=indexJ+2
                        !    coeffInter(3)=0.d0
                        !    coeffInter(4)=1.d0
                        ! else
                        !    coeffInter(3)=5.d-1
                        !    coeffInter(4)=5.d-1
                        !    ! coeffInter(3)=layKN(  k)/(layKN(k+1)+layKN(k))
                        !    ! coeffInter(4)=layKN(k+1)/(layKN(k+1)+layKN(k))
                     endif
                     if (nel.eq.1) then
                        print*, "camada ->", k
                        print*, "Beta(1)=", coeffInter(1)
                        print*, "Beta(2)=", coeffInter(2)
                        print*, "Beta(3)=", coeffInter(3)
                        print*, "Beta(4)=", coeffInter(4)
                     endif
                     ! coeff(1) =  coeffInter(1)*(1-2*xi-ncamadas)
                     ! coeff(2) = -coeffInter(1)*(1-2*xi+ncamadas)
                     ! coeff(3) =  coeffInter(1)*(2*ncamadas)
                     ! coeff(4) = -coeffInter(4)*(1-2*xi+ncamadas)
                     ! coeff(5) =  coeffInter(4)*(1-2*xi-ncamadas)
                     ! coeff(6) =  coeffInter(4)*(2*ncamadas)
                     ! coeff(7) = -coeffInter(3)*(1-2*xi+ncamadas) &
                     !            +coeffInter(2)*(1-2*xi-ncamadas) +2*ncamadas
                     ! coeff(8) =  coeffInter(3)*(1-2*xi-ncamadas) &
                     !            -coeffInter(2)*(1-2*xi+ncamadas) +2*ncamadas                        
                     ! coeff(9) = (coeffInter(2)+coeffInter(3)-2)*2*ncamadas
                     
                     coeff(1) =  coeffInter(1)*coeffInter(1)*(1-2*xi-ncamadas)
                     coeff(2) = -coeffInter(1)*coeffInter(4)*(1-2*xi+ncamadas)
                     coeff(3) = -coeffInter(1)*coeffInter(3)*(1-2*xi+ncamadas) &
                     +coeffInter(1)*coeffInter(2)*(1-2*xi-ncamadas) &
                     +coeffInter(1)*(2*ncamadas)
                     
                     coeff(4) = -coeffInter(1)*coeffInter(4)*(1-2*xi+ncamadas)
                     coeff(5) =  coeffInter(4)*coeffInter(4)*(1-2*xi-ncamadas)
                     coeff(6) = -coeffInter(3)*coeffInter(4)*(1-2*xi-ncamadas) &
                     +coeffInter(2)*coeffInter(4)*(1-2*xi+ncamadas) &
                     +coeffInter(4)*(2*ncamadas)
                     
                     coeff(7) = -coeffInter(1)*coeffInter(3)*(1-2*xi+ncamadas) &
                     +coeffInter(1)*coeffInter(2)*(1-2*xi-ncamadas) &
                     +coeffInter(1)*(2*ncamadas)
                     coeff(8) = -coeffInter(3)*coeffInter(4)*(1-2*xi-ncamadas) &
                     +coeffInter(2)*coeffInter(4)*(1-2*xi+ncamadas) &
                     +coeffInter(4)*(2*ncamadas)
                     coeff(9) = (coeffInter(2)*coeffInter(2)+ &
                     coeffInter(3)*coeffInter(3))*(1-2*xi-ncamadas) &
                     -2*coeffInter(2)*coeffInter(3)*(1-2*xi+ncamadas) &
                     +(coeffInter(2)+coeffInter(3)-1)*(4*ncamadas)
                     
                     ! Termo de divergencia
                     eleffm(posJ3,posI3) = eleffm(posJ3,posI3) + layKT(k)*layDh(k)*dot_product(di(1:nsd),dj(1:nsd))
                     
                     ! Termos de saltos e medias
                     aux = ( layKN(k)/( dble(ncamadas)*(2.d0-4.d0*xi) ) )*di(nrowsh)*dj(nrowsh)
                     eleffm(posJ1,posI1) = eleffm(posJ1,posI1) + aux*coeff(1)
                     eleffm(posJ1,posI2) = eleffm(posJ1,posI2) + aux*coeff(2)
                     eleffm(posJ1,posI3) = eleffm(posJ1,posI3) + aux*coeff(3)
                     
                     eleffm(posJ2,posI1) = eleffm(posJ2,posI1) + aux*coeff(4)
                     eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + aux*coeff(5)
                     eleffm(posJ2,posI3) = eleffm(posJ2,posI3) + aux*coeff(6)
                     
                     eleffm(posJ3,posI1) = eleffm(posJ3,posI1) + aux*coeff(7)
                     eleffm(posJ3,posI2) = eleffm(posJ3,posI2) + aux*coeff(8)
                     eleffm(posJ3,posI3) = eleffm(posJ3,posI3) + aux*coeff(9)
                     
                     ! ! Termos de saltos e medias
                     ! aux = (layKN(k)/(2.d0*dble(ncamadas)))*di(nrowsh)*dj(nrowsh)
                     ! eleffm(posJ1,posI1) = eleffm(posJ1,posI1) + aux
                     ! eleffm(posJ1,posI2) = eleffm(posJ1,posI2) - aux
                     ! eleffm(posJ2,posI1) = eleffm(posJ2,posI1) - aux
                     ! eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + aux
                     
                     ! aux = (layKN(k)/(2.d0-4.d0*xi))*di(nrowsh)*dj(nrowsh)
                     ! eleffm(posJ1,posI1) = eleffm(posJ1,posI1) - aux
                     ! eleffm(posJ1,posI2) = eleffm(posJ1,posI2) - aux
                     ! eleffm(posJ1,posI3) = eleffm(posJ1,posI3) + aux*2.d0
                     
                     ! eleffm(posJ2,posI1) = eleffm(posJ2,posI1) - aux
                     ! eleffm(posJ2,posI2) = eleffm(posJ2,posI2) - aux
                     ! eleffm(posJ2,posI3) = eleffm(posJ2,posI3) + aux*2.d0
                     
                     ! eleffm(posJ3,posI1) = eleffm(posJ3,posI1) + aux*2.d0
                     ! eleffm(posJ3,posI2) = eleffm(posJ3,posI2) + aux*2.d0
                     ! eleffm(posJ3,posI3) = eleffm(posJ3,posI3) - aux*4.d0
                  end do
               end do
            end do
         end do

         ! Rotina para condensação de graus de liberdades das pontas das descontinuidades
         do i=1, nenFratura
            ! write(*,'("lmFract(",i5,"->",i5,",",i5,",",i5,",",i5,",",i5,")")') i, estrutSistEqP_%lmFratura(1,i,nel), estrutSistEqP_%lmFratura(2,i,nel), estrutSistEqP_%lmFratura(3,i,nel), estrutSistEqP_%lmFratura(4,i,nel), estrutSistEqP_%lmFratura(5,i,nel)
            if ((estrutSistEqP_%lmFratura(1,i,nel) .ne. 0) .and. (estrutSistEqP_%lmFratura(2,i,nel) == 0) &
            .and. (estrutSistEqP_%lmFratura(3,i,nel) == 0)) then
               indexI=(i-1)*estrutSistEqP_%ndof+1
               
               eleffm(:,indexI)=eleffm(:,indexI)+eleffm(:,indexI+1)+eleffm(:,indexI+2)
               eleffm(indexI,:)=eleffm(indexI,:)+eleffm(indexI+1,:)+eleffm(indexI+2,:)
               elresf(indexI)=elresf(indexI)+elresf(indexI+1)+elresf(indexI+2)
               
               eleffm(:,indexI+1)=0.d0
               eleffm(:,indexI+2)=0.d0
               eleffm(indexI+1,:)=0.d0
               eleffm(indexI+2,:)=0.d0
               elresf(indexI+1)=0.d0
               elresf(indexI+2)=0.d0                     
            endif
         end do

         ! computation of Dirichlet b.c. contribution
         call ztest(pl,nee,zerodl)

         if(.not.zerodl) then
            call KDBC2(ELEFFM,ELRESF,pl,NEE,estrutSistEqP_%lmFratura,NEL)
         endif

         !.... assemble element stifness matrix and force array into global
         !        left-hand-side matrix and right-hand side vector
         lsym=.true.

         if(passoTempo==1) then
            if (estrutSistEqP_%optSolver=='skyline')   then
               call addlhsN (estrutSistEqP_%alhs, eleffm, estrutSistEqP_%lmFratura, estrutSistEqP_%idiag, nee, nel, diag, lsym)
            endif
            if (estrutSistEqP_%optSolver=='pardiso') then
               call addlhsCSR  (estrutSistEqP_,eleffm, estrutSistEqP_%lmFratura, nee, nel)
            endif
               
            if (estrutSistEqP_%optSolver=='hypre')   then
               call addnslHYPRE(A_HYPRE, eleffm, estrutSistEqP_%lmFratura, nee, nel)
            endif
         endif

         call addrhsN(estrutSistEqP_%brhs, elresf, estrutSistEqP_%lmFratura, nee, nel)

         deallocate(layDh)
         deallocate(layKN)
         deallocate(layKT)
         deallocate(layBeta)
      end do

   end subroutine calcCoefSistAlgPotencial_fraturas



      !==========================================================================================================
   subroutine calcCoefSistAlgPotencial_Descontinuidade(estrutSistEqP_)
      ! Calcula e acumula a matriz do elemento de descontinuidade
      ! O calculo das pressões barra são realizados utilizando as velocidades
      ! nas interfaces das camadas.
      use mGlobaisEscalares, only: transiente
      use mGlobaisEscalares, only: npintFratura, viscosidade, nrowsh
      use mGlobaisArranjos,  only: c, grav, matFratura
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shgq, shlq, shg1D_points, shl1D_intPoints
      use mfuncoesDeForma,   only: shgCodimOneSurface_points, shgSurface_points
      use mMalha,            only: local, x, numnp, nenFratura, nsd, numelFratura, conecNodaisFratura, normal
      use mGlobaisEscalares, only: deltaT, passoTempo, carregamento, carregamento_ref, saltoPressao
      use mMLayer,           only: LayerMaterial

      use mUtilSistemaEquacoes,  only: kdbc,KDBC2

      use mSolverGaussSkyline,   only: addrhsN, addlhsN
      use mSolverPardiso,        only: addlhsCSR
      use mSolverHypre

      implicit none

      type (estruturasArmazenamentoSistemaEq), intent(inout) :: estrutSistEqP_

      real*8 :: xl(nsd,nenFratura), xlAux(nsd,nenFratura)
      real*8 :: pl(estrutSistEqP_%ndof,nenFratura), plAnt(estrutSistEqP_%ndof,nenFratura)
      real*8 :: shg(nrowsh,nenFratura,npintFratura), shl(nsd,nenFratura,npintFratura)
      real*8 :: det(npintFratura), w(npintFratura)
      real*8 :: normalNoLocal(nsd,nenFratura)

      integer*4 :: nee

      real*8   :: elresf(estrutSistEqP_%ndof*nenFratura), eleffm(estrutSistEqP_%ndof*nenFratura,estrutSistEqP_%ndof*nenFratura)

      integer*4:: nel, m, l, i, j, ni, nj, k

      real*8  :: condutividade, abertura, beta, rigidezNormalInicial
      real*8  :: coef_poisson, mod_young
      real*8  :: wDet, gf1, gf2, gf3, permeabilidade
      real*8  :: djx, djy, djz, djn, dix, diy, diz, din
      logical :: diag,zerodl,quad,lsym
      real*8  :: UU, UUP(estrutSistEqP_%ndof)
      real*8  :: resistividade, di(nrowsh), dj(nrowsh), aux, xi
      integer*4 :: indexI, indexJ

      real*8 :: poisson, modYoung, bulkGrain, bulkOil, bulkRocha
      real*8 :: umPorM, porosidade
      real*8, dimension(:), allocatable :: layKN, layKT, layDh, layBeta
      real*8 :: alpha(4), coeff(9), coeffInter(4)
      integer*4:: ncamadas, posI1, posI2, posI3, posJ1, posJ2, posJ3

      nee  = nenFratura*estrutSistEqP_%ndof
      diag = .false.

      w=0.0
      shl=0.0
      ! xi=0.49999d0
      ! xi=3.d0/4.d0
      ! xi=1.d0
        
      if(nenFratura==2) then
         call shl1D_intPoints(shl, w, npintFratura, nenFratura)
      else if(nenFratura==3) then
         call shlt(shl,w,npintFratura,nenFratura)
      else if(nenFratura==4) then
         call shlq(shl,w,npintFratura,nenFratura)
      end if
         
      do nel=1, numelFratura

         ! clear stiffness matrix and force array
         eleffm=0.0
         elresf=0.0

         ! LOCALIZE COORDINATes and Dirichlet b.c.
         call local(conecNodaisFratura(:,nel), x, xl, nenFratura, nsd, nsd)
         call local(conecNodaisFratura(:,nel), estrutSistEqP_%u, pl, nenFratura, estrutSistEqP_%ndof, estrutSistEqP_%ndof)
         if(transiente.eqv..true.) then
            call local(conecNodaisFratura(:,nel),estrutSistEqP_%uTempoAnt,plAnt,nenFratura,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
         endif

         ! write(*,*) (plAnt(k,1),k=1,5)

         quad = .true.
         if (nenFratura.eq.4.and.conecNodaisFratura(3,nel).eq.conecNodaisFratura(4,nel)) quad = .false.

         call shgCodimOneSurface_points(xl, det, shl, shg, nenFratura, npintFratura, nsd)

         !....... form stiffness matrix
         !... Properties
         m = matFratura(nel)
         beta=c(6,m)
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

            ! Atualizar as pressões nos graus de liberdade da descontinuidade.
            if(transiente.eqv..true.) then
               if(saltoPressao) then
                  UUP=0.D00
                  do J=1,nenFratura
                     do k=1,ncamadas
                        UUP(k+2) = UUP(k+2) + SHG(nrowsh,J,L)*plAnt(k+2,J)
                     enddo
                  enddo
               endif
            endif

            do j=1,nenFratura
               dj(:) = shg(:,j,l)*wDet
               indexJ = (j-1)*estrutSistEqP_%ndof

               ! determinacao do termo transiente (dP/dt)
               if(transiente.eqv..true.) then
                  do k=1,ncamadas
                     elresf(indexJ+k+2)=elresf(indexJ+k+2) + layDh(k)*layBeta(k)*UUP(k+2)*dj(nrowsh)
                  enddo
               endif

               ! determinacao da parela espacial (dP/dx)
               do i=1, nenFratura
                  di(:) = shg(:,i,l)
                  indexI = (i-1)*estrutSistEqP_%ndof

                  do k=1, ncamadas
                     ! Determinacao das posicoes da camada
                     posJ1 = indexJ+1+k
                     posJ2 = indexJ+2+k
                     posJ3 = indexJ+3+k
                     posI1 = indexI+1+k
                     posI2 = indexI+2+k
                     posI3 = indexI+3+k
                                        
                     if(k.eq.1) then
                        posI1=indexI+1
                        posJ1=indexJ+1
                     endif
                     if(k.eq.ncamadas) then
                        posI3=indexI+2
                        posJ3=indexJ+2
                     endif
                     
                     coeff=0.d0
                     if(k.eq.1) then
                        coeff(1) = (2.d0*layKN(k))/layDh(k)
                        coeff(4) = -coeff(1)
                        coeff(2) = -coeff(1)
                        coeff(5) =  coeff(1)
                     else
                        coeff(2) = -(2.d0*layKN(k-1)*layKN(k))/(layKN(k-1)*layDh(k)+layKN(k)*layDh(k-1))
                        coeff(5) = -coeff(2)
                     endif
                     if(k.eq.ncamadas) then
                        coeff(8) = -(2.d0*layKN(k))/layDh(k)
                        coeff(5) =  coeff(5)-coeff(8)
                        coeff(6) =  coeff(8)
                        coeff(9) = -coeff(8)
                     else
                        coeff(8) = -(2.d0*layKN(k+1)*layKN(k))/(layKN(k+1)*layDh(k)+layKN(k)*layDh(k+1))
                        coeff(5) = coeff(5)-coeff(8)
                     endif

                     if(transiente.eqv..true.) then
                        ! //TODO : Corrigir o termo di*dj 
                        ! Termo de divergencia
                        eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + layDh(k)*layBeta(k)*di(nrowsh)*dj(nrowsh)
                        ! eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + layBeta(k)*dot_product(di(1:nsd),dj(1:nsd))
                        eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + deltaT*layKT(k)*layDh(k)*dot_product(di(1:nsd),dj(1:nsd))
                     
                        ! Termos de saltos e medias
                        aux = deltaT*di(nrowsh)*dj(nrowsh)
                        eleffm(posJ1,posI1) = eleffm(posJ1,posI1) + aux*coeff(1)
                        eleffm(posJ1,posI2) = eleffm(posJ1,posI2) + aux*coeff(2)
                        eleffm(posJ1,posI3) = eleffm(posJ1,posI3) + aux*coeff(3)
                        
                        eleffm(posJ2,posI1) = eleffm(posJ2,posI1) + aux*coeff(4)
                        eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + aux*coeff(5)
                        eleffm(posJ2,posI3) = eleffm(posJ2,posI3) + aux*coeff(6)
                        
                        eleffm(posJ3,posI1) = eleffm(posJ3,posI1) + aux*coeff(7)
                        eleffm(posJ3,posI2) = eleffm(posJ3,posI2) + aux*coeff(8)
                        eleffm(posJ3,posI3) = eleffm(posJ3,posI3) + aux*coeff(9)
                     else
                        ! Termo de divergencia
                        eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + layKT(k)*layDh(k)*dot_product(di(1:nsd),dj(1:nsd))
                     
                        ! Termos de saltos e medias
                        aux = di(nrowsh)*dj(nrowsh)
                        eleffm(posJ1,posI1) = eleffm(posJ1,posI1) + aux*coeff(1)
                        eleffm(posJ1,posI2) = eleffm(posJ1,posI2) + aux*coeff(2)
                        eleffm(posJ1,posI3) = eleffm(posJ1,posI3) + aux*coeff(3)
                     
                        eleffm(posJ2,posI1) = eleffm(posJ2,posI1) + aux*coeff(4)
                        eleffm(posJ2,posI2) = eleffm(posJ2,posI2) + aux*coeff(5)
                        eleffm(posJ2,posI3) = eleffm(posJ2,posI3) + aux*coeff(6)
                     
                        eleffm(posJ3,posI1) = eleffm(posJ3,posI1) + aux*coeff(7)
                        eleffm(posJ3,posI2) = eleffm(posJ3,posI2) + aux*coeff(8)
                        eleffm(posJ3,posI3) = eleffm(posJ3,posI3) + aux*coeff(9)
                     endif
                  end do
               end do
            end do
         end do

         ! Rotina para condensação de graus de liberdades das pontas das descontinuidades
         do i=1, nenFratura
            ! write(*,'("lmFract(",i5,"->",i5,",",i5,",",i5,",",i5,",",i5,")")') i, estrutSistEqP_%lmFratura(1,i,nel), estrutSistEqP_%lmFratura(2,i,nel), estrutSistEqP_%lmFratura(3,i,nel), estrutSistEqP_%lmFratura(4,i,nel), estrutSistEqP_%lmFratura(5,i,nel)
            if ((estrutSistEqP_%lmFratura(1,i,nel) .ne. 0) .and. (estrutSistEqP_%lmFratura(2,i,nel) == 0) &
            .and. (estrutSistEqP_%lmFratura(3,i,nel) == 0)) then
               indexI=(i-1)*estrutSistEqP_%ndof+1
               
               eleffm(:,indexI)=eleffm(:,indexI)+eleffm(:,indexI+1)+eleffm(:,indexI+2)
               eleffm(indexI,:)=eleffm(indexI,:)+eleffm(indexI+1,:)+eleffm(indexI+2,:)
               elresf(indexI)=elresf(indexI)+elresf(indexI+1)+elresf(indexI+2)
               
               eleffm(:,indexI+1)=0.d0
               eleffm(:,indexI+2)=0.d0
               eleffm(indexI+1,:)=0.d0
               eleffm(indexI+2,:)=0.d0
               elresf(indexI+1)=0.d0
               elresf(indexI+2)=0.d0                     
            endif
         end do

         ! computation of Dirichlet b.c. contribution
         call ztest(pl,nee,zerodl)

         if(.not.zerodl) then
            call KDBC2(ELEFFM,ELRESF,pl,NEE,estrutSistEqP_%lmFratura,NEL)
         endif

         !.... assemble element stifness matrix and force array into global
         !        left-hand-side matrix and right-hand side vector
         lsym=.true.

         if(passoTempo==1) then
            if (estrutSistEqP_%optSolver=='skyline')   then
               call addlhsN (estrutSistEqP_%alhs, eleffm, estrutSistEqP_%lmFratura, estrutSistEqP_%idiag, nee, nel, diag, lsym)
            endif
            if (estrutSistEqP_%optSolver=='pardiso') then
               call addlhsCSR  (estrutSistEqP_,eleffm, estrutSistEqP_%lmFratura, nee, nel)
            endif
               
            if (estrutSistEqP_%optSolver=='hypre')   then
               call addnslHYPRE(A_HYPRE, eleffm, estrutSistEqP_%lmFratura, nee, nel)
            endif
         endif

         call addrhsN(estrutSistEqP_%brhs, elresf, estrutSistEqP_%lmFratura, nee, nel)

         deallocate(layDh)
         deallocate(layKN)
         deallocate(layKT)
         deallocate(layBeta)
      end do

   end subroutine calcCoefSistAlgPotencial_Descontinuidade


   !==========================================================================================================
   subroutine calcMatrizElemento_Descontinuidade(matElem, prop, w, det, shg, ndof, nen, npint,&
      ncamada, nrowsh)
      ! Rotina para montagem da matriz de elemento de descontinuidades com múltiplas

      ! Descrição das variáveis para a rotina
      !  matElem -> matriz de elemento
      !     prop -> matriz com as propriedades das camadas (espessura, alpha, permeabilidade, beta)
      ! prop(1,:)-> espessura por camada
      ! prop(2,:)-> coeficiente alpha por camada
      ! prop(3,:)-> permeabilidade por camada
      ! prop(4,:)-> compressibilidade por camada
      !        w -> matriz de pesos para integração do elemento
      !      det -> determinantes
      !      shg -> coeficientes de interpolação
      !     ndof -> número de graus de liberdade
      !  ncamada -> número de camadas
      !      nen -> número de nós para o elemento
      !    npint -> número de pontos de integração
      !   nrowsh -> número da quantidade de derivadas mais a função (nsd+1)
      
      ! Variáveis de entradas e saídas
      INTEGER, INTENT(IN) :: ndof, nen, npint, nrowsh, ncamada 
      DOUBLE PRECISION, INTENT(IN) :: det(npint), w(npint), shg(nrowsh,nen,npint)
      DOUBLE PRECISION, INTENT(IN) :: prop(4,ncamada)
      DOUBLE PRECISION, INTENT(INOUT) :: matElem(nen*ndof,nen*ndof)
      
      ! Variáveis internas
      INTEGER :: i, j, l, posI, posJ
      DOUBLE PRECISION :: xi, wDet, aux, esc, di(nrowsh), dj(nrowsh)
      ! LOGICAL, INTENT(IN) :: transiente
      ! INTEGER :: indexI, indexJ
      ! integer*4:: nel, m, ni, nj, k

      ! Forma do perfil de pressão dentro da camada
      xi=0.49999d0
      ! xi=3.d0/4.d0
      ! xi=1.d0

      !.... loop on integration points
      matElem = 0.d0
      do l=1,npint
         wDet = w(l)*det(l)

         ! if(transiente.eqv..true.) then
         !    UUP=0.D00
         !    DO J=1,nenFratura
         !       UUP(1)  = UUP(1) + SHG(nrowsh,J,L)*PLAnt(1,J) 
         !       if(saltoPressao)  then
         !          UUP(3)  = UUP(3) + SHG(nrowsh,J,L)*PLAnt(3,J) 
         !       endif
         !    ENDDO
         ! endif

         ! Loop nos nós dos elementos 
         do j=1,nen
            dj(:) = shg(:,j,l)*wDet
            
            ! TODO:Corrigir entrada do termo transiente
            ! Determinacao do termo transiente (dP/dt)
            ! if(transiente.eqv..true.) then
            !    elresf(indexJ+3)=elresf(indexJ+3) + abertura*beta*UUP(3)*dj(nrowsh)
            ! endif
            
            ! Determinacao da parela espacial (dP/dx)
            do i=1, nen
               di(:) = shg(:,i,l)
               posI = ndof*(i-1)
               posJ = ndof*(j-1)

               ! Loop nas camadas 
               do k=1, ncamada
                  matElem(posI+k,posJ+k) = matElem(posI+k,posJ+k) + &
                     prop(3,k)*prop(1,k)*dot_product(di(1:nsd),dj(1:nsd))
                  
                  aux = ( prop(2,k)/( dble(ncamada)*(2.d0-4.d0*xi) ) )*di(nrowsh)*dj(nrowsh)

                  matElem(posI+k,posJ+k) = matElem(posI+k,posJ+k) - 4.d0*camadas*aux

                  esc=2.d0*dble(ncamada)*aux
                  call BackwardSubs_Pres(matElem,esc,nen,ncamada,posI,posJ,j,k,.TRUE.)
                  call BackwardSubs_Pres(matElem,esc,nen,ncamada,posI,posJ,j,k-1,.FALSE.)
                  call BackwardSubs_W(matElem,esc,nen,ncamada,posI,posJ,i,k)
                  call BackwardSubs_W(matElem,esc,nen,ncamada,posI,posJ,i,k-1)
                  
                  esc=(1.d0-2.d0*xi-dble(ncamada))*aux

                  esc=-(1.d0-2.d0*xi+dble(ncamada))*aux

               end do
            end do
         end do
      end do

      ! Reordenação da matriz de elemento
      
   end subroutine

   !==========================================================================================================
   subroutine BackwardSubs_Pres(matElem,escalar,nen,ncamada,posI,posJ,nno,cam,flag)
      ! Rotina para backward substitution do termo da pressão

      ! Descrição das variáveis para a rotina
      !  matElem -> matriz de elemento
      !  escalar -> escalar multiplicado na substituição backward
      !      nen -> número de nós para o elemento
      !  ncamada -> número de camadas
      !     posI -> posição i na matriz de elemento
      !     posJ -> posição j na matriz de elemento
      !      cam -> camada trabalha
      !      nno -> número do nó avaliado
      !     flag -> Flag para indicar se esta analisado em j ou j+1

      INTEGER, INTENT(IN) :: nen, ncamada, cam
      INTEGER, INTENT(IN) :: posI, posJ
      DOUBLE PRECISION, INTENT(IN) :: escalar
      DOUBLE PRECISION, INTENT(INOUT) :: matElem(nen*(ncamada+2),nen*(ncamada+2))
      LOGICAL, INTENT(IN) :: flag

      INTEGER :: k

      if (cam.eq.1) then
         matElem(posI+1,posJ+2) = matElem(posI+1,posJ+2)+escalar
         return
      else
         matElem(posI+cam+1,posJ+cam+1) = matElem(posI,posJ)+2*escalar
         
         matElem((nno-1)*(ncamada+2)+1,posJ) = matElem((nno-1)*(ncamada+2)+1,posJ)+((-1)**(cam))*escalar
      endif
   

      if (cam.ne.1) then
         do k=1,cam-1
            matElem((nno-1)*(ncamada+2)+k+1,posJ) = matElem((nno-1)*(ncamada+2)+k+1,posJ)  &
            + 2*((-1)**(cam-k))*escalar
         end do
      endif
      
   end subroutine

   !==========================================================================================================
   subroutine BackwardSubs_W(matElem,escalar,nen,ncamada,posI,posJ,nno,cam)
      ! Rotina para backward substitution do termo da função teste (W)

      ! Descrição das variáveis para a rotina
      !  matElem -> matriz de elemento
      !  escalar -> escalar multiplicado na substituição backward
      !      nen -> número de nós para o elemento
      !  ncamada -> número de camadas
      !     posI -> posição i na matriz de elemento
      !     posJ -> posição j na matriz de elemento
      !      cam -> camada trabalha
      !      nno -> número do nó avaliado

      INTEGER, INTENT(IN) :: nen, ncamada, cam
      INTEGER, INTENT(IN) :: posI, posJ
      DOUBLE PRECISION, INTENT(IN) :: escalar
      DOUBLE PRECISION, INTENT(INOUT) :: matElem(nen*(ncamada+2),nen*(ncamada+2))

      INTEGER :: k
   
      matElem(posI,posJ) = matElem(posI,posJ)+2*escalar
      matElem(posI,(nno-1)*(ncamada+2)+1) = matElem(posI,(nno-1)*(ncamada+2)+1)+((-1)**(cam))*escalar

      if(cam.ne.1)then
         do k=1,cam-1
            matElem(posI,(nno-1)*(ncamada+2)+k+1) = matElem(posI,(nno-1)*(ncamada+2)+k+1)  &
            + 2*((-1)**(cam-k))*escalar
         end do
      endif
      
   end subroutine

   !==========================================================================================================
   subroutine BackwardSubs(matElem,escalar,nen,ncamada,posI,posJ,nno,cam)
      ! Rotina para backward substitution do termo da função teste (W) e da pressão (P)

      ! Descrição das variáveis para a rotina
      !  matElem -> matriz de elemento
      !  escalar -> escalar multiplicado na substituição backward
      !      nen -> número de nós para o elemento
      !  ncamada -> número de camadas
      !     posI -> posição i na matriz de elemento
      !     posJ -> posição j na matriz de elemento
      !      cam -> camada trabalha
      !      nno -> número do nó avaliado

      INTEGER, INTENT(IN) :: nen, ncamada, cam
      INTEGER, INTENT(IN) :: posI, posJ
      DOUBLE PRECISION, INTENT(IN) :: escalar
      DOUBLE PRECISION, INTENT(INOUT) :: matElem(nen*(ncamada+2),nen*(ncamada+2))

      INTEGER :: k
   
      matElem(posI,posJ) = matElem(posI,posJ)+2*escalar
      matElem(posI,(nno-1)*(ncamada+2)+1) = matElem(posI,(nno-1)*(ncamada+2)+1)+((-1)**(cam))*escalar

      if(cam.ne.1)then
         do k=1,cam-1
            matElem(posI,(nno-1)*(ncamada+2)+k+1) = matElem(posI,(nno-1)*(ncamada+2)+k+1)  &
            + 2*((-1)**(cam-k))*escalar
         end do
      endif
      
   end subroutine

end module
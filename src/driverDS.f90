!>=================================================================================
!> \mainpage
!>
!>         Programa de elementos finitos em fortran 90 ++
!>         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!>         + implementacoes de Abimael Loula
!>         + implementacoes de Tuane Lopes
!>
!>         novo projeto: 
!>         Eduardo da Silva Castro,		ecastro@lncc.br [5]
!>
!>         LNCC/MCT
!>         Petropolis, 2.2017
!>>
!>         Propriedades e Funcionalidades:
!>>
!>         0. calculo do potencial nos pontos nodais pelo metodo de galerkin
!>           e do fluxo nodal pela tecnica de posprocessamento global (Loula, 1995) 
!>         1. resultados identicos ao programa original
!>         2. mantem parte das rotinas originais do livro do hughes: 
!>             shlt, shlq, shgt, shgq, local, 
!>             addrhs, addlhs, kdbc, load, btod, ftod, fact, back
!>         3  eliminacao do vetor A estatico e da funcao mpoint
!>         4. alocacao dinamica de memoria: stack (automaticas) + heap (allocatable)
!>         5. novos nomes, + significativoes e maiores, de variaveis e procedimentos  
!>         6. utilizacao "module":
!>             em substituicao aos commons: exigencia fortran 90 em conhecer as interfaces dos procedimentos
!>             agrupamento de variaveis e procedimentos
!>             
!>         7. obrigatoriedade de declaracao de variaveis: implicit none
!>         8. composto de varios arquivos
!>               variaveisGlobais.F90 malha.F90 utilitarios.F90
!>               leituraEscrita.F90 mInputReader.F90
!>               funcoesDeForma.F90 estruturasDadosSistEq.F90
!>               potencial.F90 fluxo.F90
!>               solverGaussSkyline.F90 solverHypre.F90 solverPardisoCSR.F90
!>               utilSistemaEquacoes.F90 driverDS.F90
!>           e varios diretorios
!>               fontes, include, bin 
!>         9. extensao "F90" dos arquivos com maiuscula sendo importante para uso de opcao de preprocessamento
!>        10. 3 opcoes de solvers: 2 diretos (gauss original e pardiso) e 1 iterativo (HYPRE)
!>        11. leitura de dados usando palavras chaves (desenvolvido pela DeepSoft)
!>        10. utilizacao do comando make no linux
!>              porem com possibilidade de unir os arquivos em um unico
!>              colocando-os na ordem que aparecem no item 8 
!>              (exigencia por existencia de dependencia entre "module"
!>=================================================================================
program poisson
   use mGlobaisEscalares, only: exec
   use mLeituraEscrita,   only: abrirArquivos, fecharArquivos
   use mSolverHypre,      only: inicializarMPI, finalizarMPI
   use mSolverHypre,      only: myid, num_procs, mpi_comm
   use mVisualizacaoHdf5, only: hdf5_basic_structure, hdf5_close_structure
   use mMLayer,           only: MLayer_Finalizar

   implicit none

   ! Inicializacao
#ifdef withHYPRE
   call inicializarMPI(myid, num_procs, mpi_comm)
#endif
   call abrirArquivos()
   call hdf5_basic_structure()

   ! Preprocessamento
   print*, " call preprocessadorDS ()"
   call preprocessamentoDS()

   ! Processamento
   print*, " call processamento ()"
   call processamento()

   ! Encerramento de arquivos abertos
   call MLayer_Finalizar()
   call fecharArquivos()
   call hdf5_close_structure()     
#ifdef withHYPRE
   call finalizarMPI()
#endif

end program poisson
!=================================================================================
!> Routine for the preprocessement 
subroutine preprocessamentoDS ()
   use mInputReader,      only: readInputFileDS
   use mInputReader,      only: leituraGeracaoCoordenadasDS
   use mInputReader,      only: leituraCodigosCondContornoDS
   use mInputReader,      only: leituraValoresCondContornoDS, leituraValoresIniciaisDS
   ! use mInputReader,      only: leituraCodigosCondContornoDS_PressaoFratura
   ! use mInputReader,      only: leituraCodigosCondContornoDS_PressaoFratura_3D
   use mInputReader,      only: leituraCodigosCondContornoDS_PressaoDescontinuidade
   use mInputReader,      only: leituraCodigosCondContornoDS_PressaoDescontinuidade_3D
   use mInputReader,      only: leituraGeracaoConectividadesDS,leituraValoresCondContornoDS_PressaoFratura
   use mInputReader,      only: leituraMatrizInteiros, leituraMatrizInteiros2

   use mInputReader,      only: leituraCodigosCondContornoFluxoNormalDS
   
   use mGlobaisArranjos,  only: etime, title, mat, c, matFratura, tiposElemBordas, aberturasBorda, alturasBorda
   use mGlobaisArranjos,  only: listaSolverDisponivel, campoPermeabilidade, matDescontNaBorda, tipoBordaElemDescontNaBorda
   use mGlobaisArranjos,  only: NelmDescNoBordo
   use mGlobaisEscalares, only: exec, iprtin, npint,transiente
   use mGlobaisEscalares, only: deltaT, numPassosTempo, tempoTotalSimulacao, saltoPressao
   use mMalha,            only: numnp, numel, nen, nsd, numelFratura, numnpFratura, nenFratura, numelBordas
   use mMalha,            only: x, conecNodaisElem,conecNodaisFratura, conecNodaisReduzido, conecNodaisBordas
   use mMalha,            only: normal, normalElem, gerarNormais, gerarNormaisElemFratura, gerarNumnpFratura
   use mMalha,            only: gerarNormaisElemBordas, conectIntegroParaReduzido
   use mLeituraEscrita,   only: iin, iecho, icoords, echo
   use mLeituraEscrita,   only: prntel, iparaviewPerm, iparaviewMat, escreverArqParaviewVector, escreverArqParaviewMateriais
     
   use mPotencial,        only: estrutSistEqP
   use mFluxo,            only: estrutSistEqF, estrutSistEqFluxoNormal
   use mFluxo,            only: velocidadeBordas, velocidadeBordasFratura, calcularInterseccoesDesconBorda
   use mFluxo,            only: velocidadeBordasFalhas, vm_falha, vm_bloco, vm_integro, vm_reduzido

   use mVisualizacaoHdf5, only: hdf5_dataset_integer, hdf5_dataset_real, hdf5_dataset_double
   use mVisualizacaoHdf5, only: hdf5_dataset_dvec, hdf5_structure_addMeshFalha, hdf5_structure_addPropFalha
   use mVisualizacaoHdf5, only: name_struct, p_tree_hdf5, p_root_hdf5
   use mMLayer,           only: nTypesMLayers, LayerMaterial
   use mMLayer,           only: MLayer_Inic

   use mMalha,            only: normalElemBorda, numelDescontNaBorda, conecNodaisDesconNaBorda

   use mTransmissividade, only: trans_lm, trans_nelm, Material_Ignore



   implicit none

   character(len=50) keyword_name
   logical :: flagGerarConectividadesDesconNaBorda
   ! logical :: simetria
   integer :: i,j,k!, cont
   integer :: noF,noB, nelB, nelF, m, nel, aux
   integer :: itipoBordaElemDescontNaBorda
   real*8  :: abertura_ref, no(2)
   !---------------------------------------------------------------------------------
   !.... input phase       
   call readInputFileDS()
   call readSetupPhaseDS()

   ! Leitura dos ponteiros para os materiais das falhas
   call MLayer_Inic(estrutSistEqP%ndof)
   keyword_name = "lista_prop_desc"
   call leituraMatrizInteiros(keyword_name, LayerMaterial, nTypesMLayers, estrutSistEqP%ndof)
     
   ! saltoPressao= .false.
   saltoPressao= .true.            
      
   call identificaSolversDisponiveis(listaSolverDisponivel)
   call verificarSolver(estrutSistEqP%optSolver, listaSolverDisponivel)
   call verificarSolver(estrutSistEqFluxoNormal%optSolver, listaSolverDisponivel)
   
   etime = 0.0
   !---------------------------------------------------------------------------------
   !.... initialization phase
   !.... set memory pointers for static data arrays, and call associated input routines
   !.... Numero de graus de liberdade da analise, caso sem descontinuidade ndof=1
   if (.not.((numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) ) then
      estrutSistEqP%ndof=1
   end if
   estrutSistEqF%ndof=nsd
   estrutSistEqFluxoNormal%ndof=1

   !---------------------------------------------------------------------------------
   !.... input coordinate data
   allocate(          x(nsd,numnp));              x=0.0
   write(*,*) "call leituraGeracaoCoordenadasDS"
   call leituraGeracaoCoordenadasDS(x,nsd,numnp, icoords, iprtin)
   call hdf5_dataset_double(name_struct(1), x, nsd, numnp, p_root_hdf5(5))

   ! Impressao da variavel X (Posicao espacial)
   ! do i=1,numnp
   !    write(*,'("x(",I5,",",F8.5,",",F8.5,")")') i, x(1,i), x(2,i)
   ! enddo

   !---------------------------------------------------------------------------------
   !.....generation of conectivities
   allocate(mat(numel))
   allocate(conecNodaisElem(nen,numel))
   write(*,*) "call leituraGeracaoConectividadesDS", nen, numel
   keyword_name = "conectividades_nodais"
   call leituraGeracaoConectividadesDS(keyword_name, conecNodaisElem,mat,nen)
   call hdf5_dataset_integer(name_struct(1), conecNodaisElem-1, nen, numel, p_root_hdf5(6)) 

   ! Impressao da variavel conecNodais ( conectividades - macico rochoso )
   ! do i=1,numel
   !    write(*,'("elm(",i5,")=")',advance='no') i
   !    do j=1,nen
   !       write(*,'(i5)',advance='no') conecNodaisElem(j,i)
   !    enddo
   !    write(*,*)
   ! enddo
      
   !---------------------------------------------------------------------------------
   !.....generation of descontinuity conectivities
   if (numelFratura > 0) then
      allocate(matFratura(numelFratura))
      allocate(conecNodaisFratura(nenFratura,numelFratura))
      allocate(conecNodaisReduzido(nenFratura,numelFratura))
      write(*,*) "call leituraGeracaoConectividadesDS(fraturas)"
      keyword_name = "conectividades_nodais_fratura"
      call leituraGeracaoConectividadesDS(keyword_name, conecNodaisFratura,matFratura,nenFratura)
      call conectIntegroParaReduzido()  
      call hdf5_structure_addMeshFalha(name_struct(2))
      ! Impressao da variavel conecNodaisFratura ( conectividades - falha ou fratura )
      ! do i=1,numelFratura
      !    write(*,'("elm(",i5,")=")',advance='no') i
      !    do j=1,nenFratura
      !       write(*,'(i5)',advance='no') conecNodaisFratura(j,i)
      !    enddo
      !    write(*,*)
      ! enddo
   endif

   !---------------------------------------------------------------------------------
   !.....materal ignore list for transmissibility
   if (Material_Ignore%nelm.gt.0) then
      allocate(Material_Ignore%list(Material_Ignore%nelm))
      keyword_name = "transmissividade_mat_ignore"
      call leituraMatrizInteiros2(keyword_name, Material_Ignore%list, Material_Ignore%nelm, 1)
   end if

   !---------------------------------------------------------------------------------
   !.....generation of transmissibility conectivities
   allocate(trans_lm(trans_nelm,2))
   keyword_name = "transmissividade_conectividades"
   call leituraMatrizInteiros(keyword_name, trans_lm, trans_nelm, 2)
   
   ! do i=1,trans_nelm
   !    WRITE(*,'(i2,2x,2("LM(",i2,")=",i5,3x))') i, (j,trans_lm(i,j),j=1,2)
   ! end do

      
   !---------------------------------------------------------------------------------
   !.....generation of normal vectors pointwise
   call gerarNormais()
   if (numelFratura > 0) then
      ! generation of normal vectors elementwise
      call gerarNormaisElemFratura()
      ! generation of numnpFratura and nosFratura array
      call gerarNumnpFratura()
   end if

   ! Impressao da variavel conecNodaisFratura ( conectividades - falha ou fratura )
   ! do i=1,numnp
   !    write(*,'("normal(",i5,")=")',advance='no') i
   !    do j=1,nsd
   !       write(*,'(E10.3)',advance='no') normal(j,i)
   !    enddo
   !    write(*,*)
   ! enddo
   ! do i=1,numelFratura
   !    write(*,'("normalElem(",i5,")=")',advance='no') i
   !    do j=1,nsd
   !       write(*,'(E10.3)',advance='no') normalElem(j,i)
   !    enddo
   !    write(*,*)
   ! enddo
      
   !---------------------------------------------------------------------------------
   !.....reading of the boundary conditiong for the model      
   allocate(tiposElemBordas(numelBordas)); tiposElemBordas=0
   allocate(conecNodaisBordas(nenFratura,numelBordas))
   write(*,*) "call leituraGeracaoConectividadesDS(bordas)", numelBordas, nenFratura
   keyword_name = "conectividades_nodais_bordas"
   call leituraGeracaoConectividadesDS(keyword_name, conecNodaisBordas, tiposElemBordas, nenFratura)

   ! Impressao da variavel tiposElemBordas ( A qual borda o elemento pertence )
   ! do i=1,numelBordas
   !    write(*,'("borda(",i5,"->",i5,")")') i,tiposElemBordas(i)
   ! enddo
   ! Impressao da variavel conecNodaisBordas ( conectividades dos elementos das bordas )
   ! do i=1,numelBordas
   !    write(*,'("borda(",i5,")=")',advance='no') i
   !    do j=1,nenFratura
   !       write(*,'(i5)',advance='no') conecNodaisBordas(j,i)
   !    enddo
   !    write(*,*)
   ! enddo

   ! Acrescentar para o calculo da vazao consistente
   if (numelFratura > 0) then
      write(*,*) "call calcularInterseccoesDesconBorda", numelBordas, nenFratura
      if(nsd==3) call calcularInterseccoesDesconBorda(flagGerarConectividadesDesconNaBorda)
      if ( flagGerarConectividadesDesconNaBorda .eqv. .true. ) then            
         allocate(matDescontNaBorda(numelDescontNaBorda));  matDescontNaBorda=0
         allocate(conecNodaisDesconNaBorda(nenFratura-1, numelDescontNaBorda))            
         write(*,*) "call leituraGeracaoConectividadesDS(fraturasNaBorda)", numelDescontNaBorda, nenFratura-1      
         keyword_name = "conectividades_nodais_fraturas_borda"
         call leituraGeracaoConectividadesDS(keyword_name, conecNodaisDesconNaBorda(:,:), matDescontNaBorda(:), nenFratura-1) 
         allocate(tipoBordaElemDescontNaBorda(numelDescontNaBorda))
         itipoBordaElemDescontNaBorda = 50
         open(unit=itipoBordaElemDescontNaBorda,  file='tipoBordaElemDescontNaBorda.dat')!, status='old'
         do nel=1, numelDescontNaBorda
            read(itipoBordaElemDescontNaBorda, "(2I10)") aux, tipoBordaElemDescontNaBorda(nel)
         end do   
         close(itipoBordaElemDescontNaBorda)
      endif
   endif

   !---------------------------------------------------------------------------------
   !.....generation of the normal for the boundary conditiong
   call gerarNormaisElemBordas()
   if (iprtin.eq.0) call prntel(mat,conecNodaisElem,nen,numel,1_4)
   ! Impressao das normais no bordo
   ! do i=1,numelBordas
   !    write(*,'("normal(",i5,")=")',advance='no') i
   !    do j=1,nsd
   !       write(*,'(E11.3)',advance='no') normalElemBorda(j,i)
   !    enddo
   !    write(*,*)
   ! enddo

   !---------------------------------------------------------------------------------
   !.... input boundary condition data and establish equation numbers
   allocate(estrutSistEqP%u (estrutSistEqP%ndof,numnp));
   allocate(estrutSistEqP%uTempoAnt(estrutSistEqP%ndof,numnp));
   allocate(estrutSistEqP%id(estrutSistEqP%ndof,numnp));
   allocate(estrutSistEqP%listaElemPorNo(nen,numnp));
   estrutSistEqP%u=0.0
   estrutSistEqP%uTempoAnt=0.0
   estrutSistEqP%id=0
   estrutSistEqP%listaElemPorNo=0
   write(*,*) "call leituraCodigosCondContornoDS(idPotencial"
   keyword_name = "codigos_cond_contorno_potencial"
   if ((numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) then
      if(nsd==2) call leituraCodigosCondContornoDS_PressaoDescontinuidade(keyword_name, &
               estrutSistEqP%id, estrutSistEqP%ndof, numnp, estrutSistEqP%neq, iecho, iprtin)
      if(nsd==3) call leituraCodigosCondContornoDS_PressaoDescontinuidade_3D(keyword_name, estrutSistEqP%id, &
               estrutSistEqP%ndof, numnp, estrutSistEqP%neq, iecho, iprtin)
      ! //TODO: Corrigir rotina leituraCodigosCondContornoDS_PressaoDescontinuidade_3D para o caso de falha.
      ! ----------------------------------------------------------------------------------
      ! a alteração para identificar falhas multicamadas foi feita apenas para o caso 2d
      ! ----------------------------------------------------------------------------------
      ! if(nsd==2) call leituraCodigosCondContornoDS_PressaoFratura(keyword_name, estrutSistEqP%id, &
      !          estrutSistEqP%ndof, numnp, estrutSistEqP%neq, iecho, iprtin)
      ! if(nsd==3) call leituraCodigosCondContornoDS_PressaoFratura_3D(keyword_name, estrutSistEqP%id, &
      !          estrutSistEqP%ndof, numnp, estrutSistEqP%neq, iecho, iprtin)
   else                                    
      call leituraCodigosCondContornoDS(keyword_name, estrutSistEqP%id, estrutSistEqP%ndof, numnp, &
               estrutSistEqP%neq, iecho, iprtin)
   endif
   ! Impressao da variavel ID ( numero da equacao )
   ! do i=1,numnp
   !    write(*,'("id(",i5,")->")',advance='no') i
   !    do j=1,estrutSistEqP%ndof
   !       write(*,'(i5)',advance='no') estrutSistEqP%id(j,i)
   !    enddo
   !    WRITE(*,*) ''
   ! enddo

   !---------------------------------------------------------------------------------
   !.... input boundary condition data and establish equation numbers for consistent flux
   allocate(estrutSistEqFluxoNormal%u (estrutSistEqFluxoNormal%ndof,numnp));
   allocate(estrutSistEqFluxoNormal%uTempoAnt (estrutSistEqFluxoNormal%ndof,numnp));
   allocate(estrutSistEqFluxoNormal%id(estrutSistEqFluxoNormal%ndof,numnp));
   estrutSistEqFluxoNormal%u=0.0
   estrutSistEqFluxoNormal%uTempoAnt=0.0
   estrutSistEqFluxoNormal%id=0
   call leituraCodigosCondContornoFluxoNormalDS(estrutSistEqFluxoNormal%id, &
            estrutSistEqFluxoNormal%ndof, estrutSistEqP%id, estrutSistEqP%ndof, numnp, &
            estrutSistEqFluxoNormal%neq)
   ! Impressao da variavel ID ( numero da equacao )
   ! do i=1,numnp
   !    write(*,'("id(",i5,")->")',advance='no') i
   !    do j=1,estrutSistEqFluxoNormal%ndof
   !       write(*,'(i5)',advance='no') estrutSistEqFluxoNormal%id(j,i)
   !    enddo
   !    WRITE(*,*) ''
   ! enddo

   if (estrutSistEqP%nlvect.gt.0) then
      allocate(estrutSistEqP%f(estrutSistEqP%nlvect,estrutSistEqP%ndof,numnp))
      ! Alterado para fluxo consistente
      ! allocate(estrutSistEqP%f(estrutSistEqP%ndof, numnp, estrutSistEqP%nlvect))
      estrutSistEqP%f = 0.0
      write(*,*) "call leituraValoresCondContornoDS(fPotencial,ndofP,numnp,0,nlvectP,iprtin)"
      keyword_name = "valores_cond_contorno_potencial"
      if ((numelFratura > 0) .AND. (saltoPressao .EQV. .true.) ) then    
         call leituraValoresCondContornoDS_PressaoFratura (keyword_name, estrutSistEqP%f, &
                  estrutSistEqP%ndof, numnp, 1_4, estrutSistEqP%nlvect, iprtin)
      else  
         call leituraValoresCondContornoDS(keyword_name, estrutSistEqP%f, estrutSistEqP%ndof, &
                  numnp, 1_4, estrutSistEqP%nlvect, iprtin)       
      endif
   end if
   ! Impressao da variavel ID ( numero da equacao )
   ! do j=1,estrutSistEqP%nlvect
   !    do i=1,numnp
   !       WRITE(*,'("force(",i5,")")',advance='no') i
   !       do k=1,estrutSistEqP%ndof   
   !          write(*,'(E8.3)',advance='no') estrutSistEqP%f(j,k,i)
   !       enddo
   !       WRITE(*,*) ""
   !    enddo
   !    write (*,*)
   ! enddo

   !---------------------------------------------------------------------------------
   !.... Alocacao e leitura da estimativa inicial da pressao em cada ponto
   allocate(estrutSistEqP%uIterAnt(estrutSistEqP%ndof,numnp))
   estrutSistEqP%uIterAnt = 0.0d0
   write(*,*) "call leituraValoresIniciaisDS(pressaoIterAnt)"
   keyword_name = "valores_iteracao_inicial_potencial"
   call leituraValoresIniciaisDS(keyword_name, estrutSistEqP%uIterAnt, estrutSistEqP%ndof, numnp)

   ! Impressao da variavel "uIterAnt"
   ! do i=1,numnp
   !    WRITE(*,'("uIteAnt(",i5,")->")',advance='no') i
   !    do j=1,estrutSistEqP%ndof
   !       write(*,'(E11.3)',advance='no') estrutSistEqP%uIterAnt(j,i)
   !    enddo
   !    WRITE(*,*)
   ! enddo
      
   !---------------------------------------------------------------------------------
   !.... Alocacao e leitura da pressao de referencia em cada ponto
   allocate(estrutSistEqP%uRef(estrutSistEqP%ndof,numnp))
   estrutSistEqP%uRef = 0.0d0
   write(*,*) "call leituraValoresIniciaisDS(pressao_ref)"
   keyword_name = "valores_potencial_ref"
   call leituraValoresIniciaisDS(keyword_name, estrutSistEqP%uRef, estrutSistEqP%ndof, numnp)

   ! Impressao da variavel "uRef"
   ! do i=1,numnp
   !    WRITE(*,'("uRef(",i5,")")',advance='no') i
   !    do j=1,estrutSistEqP%ndof
   !       write(*,'(E11.3)',advance='no') estrutSistEqP%uRef(1,i)
   !    enddo
   !    WRITE(*,*)
   ! enddo

   !---------------------------------------------------------------------------------
   !.... input nodal force and prescribed kinematic boundary-value data
   allocate(estrutSistEqF%u (estrutSistEqF%ndof,numnp)); estrutSistEqF%u =0.0
   allocate(estrutSistEqF%id(estrutSistEqF%ndof,numnp)); estrutSistEqF%id=0
   write(*,*) "call leituraCodigosCondContornoDS(idFluxo,"
   keyword_name = "codigos_cond_contorno_fluxo"
   call leituraCodigosCondContornoDS(keyword_name, estrutSistEqF%id, estrutSistEqF%ndof, numnp, &
            estrutSistEqF%neq, iecho, iprtin)

   ! Impressao da variavel
   ! write(*,'("ndof(",i5,")")') estrutSistEqF%ndof
   ! do i=1,numnp
   !    WRITE(*,'("id(",i5,")->")',advance='no') i
   !    do j=1,estrutSistEqF%ndof
   !       write(*,'(i5)',advance='no') estrutSistEqF%id(j,i)
   !    enddo
   !    WRITE(*,*)
   ! enddo
   ! do i=1,numnp
   !    WRITE(*,'("u(",i5,")->")',advance='no') i
   !    do j=1,estrutSistEqF%ndof
   !       write(*,'(E11.3)',advance='no') estrutSistEqF%u(j,i)
   !    enddo
   !    WRITE(*,*)
   ! enddo

   if (estrutSistEqF%nlvect.gt.0)  then
      allocate(estrutSistEqF%f(estrutSistEqF%nlvect,estrutSistEqF%ndof,numnp))
      estrutSistEqF%f = 0.0
      write(*,*) 'call leituraValoresCondContornoDS(fFluxo,ndofF,numnp,0,nlvectF,iprtin)'
      keyword_name = "valores_cond_contorno_fluxo"
      call leituraValoresCondContornoDS(keyword_name, estrutSistEqF%f, estrutSistEqF%ndof, numnp, 1_4, &
               estrutSistEqF%nlvect, iprtin)
   end if

   ! Impressao da variavel 
   ! WRITE(*,*) "estrutSistEqF%nlvect->", estrutSistEqF%nlvect
   ! do j=1,estrutSistEqF%nlvect
   !    do i=1,numnp
   !       WRITE(*,'("force(",i5,")->")',advance='no') i
   !       do k=1,estrutSistEqF%ndof
   !          write(*,'(E11.3)',advance='no') estrutSistEqF%f(j,k,i)
   !       enddo
   !       WRITE(*,*)
   !    enddo
   !    write (*,*)
   ! enddo
     
   !---------------------------------------------------------------------------------
   !.... input element data
   write (*,*) "call leituraParamNumericosPropFisicaDS()"
   call leituraParamNumericosPropFisicaDS()
     
   !---------------------------------------------------------------------------------
   !.... entender melhor o que faz esse trecho do código, só está sendo utilizado no cálculo da vazão
   allocate(aberturasBorda(nenFratura,numelBordas))
   ! allocate(NelmDescNoBordo(nenFratura,numelBordas))
   allocate(alturasBorda(2,numelBordas))
   allocate(velocidadeBordas(nenFratura,numelBordas))
   allocate(velocidadeBordasFratura(nenFratura,numelBordas))
   allocate(velocidadeBordasFalhas(nenFratura,numelBordas,(estrutSistEqP%ndof-2))) ! 3=numero de camadas

   ! Para o fluxo consistente
   ! ! DESCOMENTAR PARA ANALITICO E SEMIANALITICO 
   ! write (*,*) "calculando novas aberturas"     
   ! call calcularAberturasFratura()   
   
   ! ! DESCOMENTAR PARA HIDRODINÂMICA
   ! write (*,*) "montando aberturasBorda"   
   ! if (numelFratura > 0 ) call montarAberturasBorda()

   allocate(vm_integro(nsd,numnp))
   allocate(vm_reduzido(nsd,numnpFratura,(estrutSistEqP%ndof-2)))

   allocate(vm_bloco(nen,numel,nsd,1))
   allocate(vm_falha(nenFratura,numelFratura,nsd,(estrutSistEqP%ndof-2))) ! 3=numero de camadas

   write(*,*) "valu", nenFratura,numelFratura,nsd,(estrutSistEqP%ndof-2)

   ! Substituido por função  acima ------------
   aberturasBorda = 0d0
   do nelF=1,numelFratura
      m = matFratura(nelF)
      abertura_ref=0.
      do j=1,LayerMaterial(m,1)
         abertura_ref = abertura_ref+c(2,LayerMaterial(m,j+1))
      enddo
           
      do i=1,nenFratura
         noF = conecNodaisFratura(i,nelF)
         do nelB=1,numelBordas
            do j=1,nenFratura
               noB = conecNodaisBordas(j,nelB)
               if (noF == noB) then
                  aberturasBorda(j,nelB) = abertura_ref
               end if
            end do
         end do
      end do
   end do
   ! --------------------------------

   !---------------------------------------------------------------------------------
   !.... Cria campo de permeabilidade
   print*, "Mapenado campos de permabilidade"
   allocate(campoPermeabilidade(nsd,numel))
   do nel=1, numel
      m=mat(nel)
      campoPermeabilidade(1,nel) = c(3,m)
      campoPermeabilidade(2,nel) = c(4,m)
      if(nsd==3)campoPermeabilidade(3,nel) = c(5,m)
   end do

   !---------------------------------------------------------------------------------
   !.... imprime campo com permeabilidade ( mapa de permeabilidade e id do material )
   print*, "escrevendo materiais e permeabilidade"
   open(unit=iparaviewPerm ,file= './out/resuladtoPerm.vtk')
   call escreverArqParaviewVector(campoPermeabilidade, nsd, numnp, nen, &
            conecNodaisElem, 1, trim('perm'), len(trim('perm')), iparaviewPerm)
                    
   open(unit=iparaviewPerm ,file= './out/resultadoMaterial.vtk')
   call escreverArqParaviewMateriais(mat, nen, conecNodaisElem, trim('material'), len(trim('material')), iparaviewPerm)
   print*, " apos escrevendo materiais e permeabilidade"

   ! escrita do arquivo hdf5
   call hdf5_dataset_integer( "0", mat, 1, numel, p_tree_hdf5(5))
   call hdf5_dataset_dvec(  "0", campoPermeabilidade, nsd, numel, p_tree_hdf5(6))
   call hdf5_structure_addPropFalha("0")

   !---------------------------------------------------------------------------------
   !.... Formatos para escrita
   !  1000 format(20a4)
   !  2000 format(4i10,i10,11i10)
   !  3000 format(5x,&
   !      ' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,&
   !      ' execution code  . . . . . . . . . . . . . . (exec  ) = ',i10//5x,&
   !      '    eq. 0, data check                                   ',   //5x,&
   !      '    eq. 1, execution                                    ',  //5x,&
   !      ' input data print code . . . . . . . . . . . (iprtin ) = ',i10//5x,&
   !      '    eq. 0, print nodal and element input data           ',   //5x,&
   !      '    eq. 1, do not print nodal and element input data    ',   //5x, &
   !      ' number of space dimensions  . . . . . . . . (nsd    ) = ',i10)
   !  4000 format(5x,&
   !      ' number of nodal points  . . . . . . . . . . (numnp  ) = ',i10//5x,&
   !      ' number of elements      . . . . . . . . . . (numel  ) = ',i10//5x,&
   !      ' number of potencial load vectors  . . . . . (nlvectP) = ',i10//5x,&
   !      ' number of fluxos   load vectors   . . . . . (nlvectF) = ',i10//5x)
end subroutine preprocessamentoDS

!==========================================================================================================    

subroutine montarAberturasBorda()
    
   use mGlobaisArranjos,  only: aberturasBorda, c, matFratura
   use mLeituraescrita,   only: iaberturasBorda, iaberturasBordaOut         
   use mMalha,            only: numelFratura, nenFratura, numelBordas, conecNodaisFratura, conecNodaisBordas
   
   use mMLayer,           only: LayerMaterial
   
   implicit none

   integer :: nelF, m, noF, nelB, noB, i, j
   real*8  :: abertura_ref
 
   aberturasBorda = 0d0

   open(unit=iaberturasBorda,  file='aberturasBorda.dat', status='old', err=100) 
   print*, "Encontrado arquivo aberturasBorda.dat"
   
   do nelB=1, numelBordas
      if (nenFratura==2) read(iaberturasBorda, '(2es20.8)')  aberturasBorda(:,nelB)
      if (nenFratura==3) read(iaberturasBorda, '(3es20.8)')  aberturasBorda(:,nelB)
      if (nenFratura==4) read(iaberturasBorda, '(4es20.8)')  aberturasBorda(:,nelB)
   end do   
   
   close(iaberturasBorda)
   
   return
   
   100 do nelF=1,numelFratura
      m = matFratura(nelF)
      abertura_ref=0.
      do j=1,LayerMaterial(m,1)
         abertura_ref = abertura_ref+c(2,LayerMaterial(m,j+1))
      enddo
      ! write(*,*) 'abertura_ref',abertura_ref
        
      do i=1,nenFratura
         noF = conecNodaisFratura(i,nelF)
         do nelB=1,numelBordas
            do j=1,nenFratura
               noB = conecNodaisBordas(j,nelB)
               if (noF == noB) then
                  aberturasBorda(j,nelB) = abertura_ref
               end if
            end do
         end do
      end do         
   end do
   
   print*, "numelBordas", numelBordas
   open(unit=iaberturasBordaOut,  file='aberturasBorda.dat', status='new', err=200)
   do nelB=1, numelBordas
       if (nenFratura==2) write(iaberturasBordaOut,'(2es20.8)') aberturasBorda(:,nelB)
       if (nenFratura==3) write(iaberturasBordaOut,'(3es20.8)') aberturasBorda(:,nelB)
       if (nenFratura==4) write(iaberturasBordaOut,'(4es20.8)') aberturasBorda(:,nelB)
   end do 
   close(iaberturasBordaOut)
   
   return
   
   200 write(*,*) ' arquivo: aberturasBorda.dat NAO pode ser criado'      
   ! do i=1,numelBordas
   !    write(950, '(3es20.8)') aberturasBorda(:, i)
   ! end do
end subroutine montarAberturasBorda

!==========================================================================================================    

subroutine mapearPermeabilidadeODA(conecNodaisElem, numel, numnp, nsd, nen)
   use mMalha,            only: x
   use mGlobaisEscalares, only: nrowsh
   use mGlobaisArranjos,  only: c, mat, tiposElemBordas, campoPermeabilidade
   use mMalha,            only: local, normal, normalElemBorda, numelBordas, nenFratura, conecNodaisBordas
   use mFuncoesDeForma,   only: shlt, shlq, shlq3d, shlt3D
   use mPotencial,        only: estrutSistEqP

   implicit none
    
   integer, intent(in) :: numel, numnp, nsd, nen
   integer, intent(in) :: conecNodaisElem(nen, numel)
    
   real*8 :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1)
   real*8 :: normalNoLocal(nsd,nen)
   real*8 :: xl(nsd,nen), domX, domY, domZ
   integer :: numParticoes, tamIntervalos, i, j, nel, m
   real*8 :: pontoInicialX, pontoInicialY, pontoInicialZ
   real*8 :: px, py, pz
   integer :: particoesEmX, particoesEmY, particoesEmZ, numnpParticoes, cont, cont2, tamanhoDaMatriz
   real*8, allocatable :: ptsIntervalos(:,:), permODA(:,:), permODATemp(:,:), permODATemp2(:,:)
      

    
   tamIntervalos=10 !2 para Arapua, 10 para cenario 3
   ! tamIntervalos=2 
      
   ! Calculando o tamanho do dominio
   domX=maxval(x(1,:))-minval(x(1,:))
   domY=maxval(x(2,:))-minval(x(2,:))
   if(nsd==3)domZ=maxval(x(3,:))-minval(x(3,:))
   print*, "Dominio", domX, domY, domZ
      
   ! Calculando o número de partições em cada direção
   particoesEmX=int(domX/tamIntervalos)
   particoesEmY=int(domY/tamIntervalos)
   if(nsd==3)particoesEmZ=int(domZ/tamIntervalos)
   print*, particoesEmX,particoesEmY !100, 87
      
   ! Número total de partições
   numParticoes=particoesEmX*particoesEmY
   if(nsd==3) numParticoes=numParticoes*particoesEmZ
   print*, "numParticoes", numParticoes

   ! Define os pontos das janelas
   allocate(ptsIntervalos(4,numParticoes))
   cont=0
   pontoInicialX=minval(x(1,:))
   pontoInicialY=minval(x(2,:))
   do j=1, particoesEmY
      py=pontoInicialY+(j-1)*tamIntervalos
      do i=1, particoesEmX
         cont=cont+1
         px=pontoInicialX+(i-1)*tamIntervalos
         ptsIntervalos(1,cont)=px
         ptsIntervalos(2,cont)=py
         ptsIntervalos(3,cont)=px+tamIntervalos
         ptsIntervalos(4,cont)=py+tamIntervalos
      end do
   end do

   ! Leitura do arquivo de permeabilidades ODA vindos do Matlab
   tamanhoDaMatriz=numParticoes*2
   allocate(permODATemp(2,tamanhoDaMatriz))
   allocate(permODATemp2(2,200))
   print*, "tamanhoDaMatriz", tamanhoDaMatriz
   ! stop
      
   open(75,file="permeabilidadeODA.inc")
   do j=1, 2
      read(75,*) (permODATemp(j,i),i=1,tamanhoDaMatriz)
   end do
   close(75)

   ! esse trecho é para o Arapua
   ! open(75,file="permeabilidadeODA_2.inc")
   ! do j=1, 2
   !    read(75,*) (permODATemp2(j,i),i=1,200)
   ! end do
   ! close(75)
      
   ! Separação das permeabilidades de intersse, Kx e Ky, descarte de Kxy e Kyx
   allocate(permODA(2,numParticoes))
   cont =1
   do i=1, tamanhoDaMatriz, 2
      permODA(1,cont)=permODATemp(1,i) 
      permODA(2,cont)=permODATemp(2,i+1) 
      cont =cont +1
   enddo
   deallocate(permODATemp)
      

   ! Usei essa parte para o Arapua.
   ! cont=1
   ! do i=1, numParticoes
   !    if(mod(i,100)==0) then
   !       permODA(1,i)=permODATemp2(1,cont)
   !       permODA(2,i)=permODATemp2(2,cont+1)
   !       cont=cont+2
   !    endif
   ! end do
   ! stop
      
   ! Ordenação. No matlab numera de cima para baixo e nós de baixo para cima.
   allocate(permODATemp(2,numParticoes))
   permODATemp=permODA
      
   cont=1
   do j=1, particoesEmY
      cont2=numParticoes-(particoesEmX*j)+1
      do i=1, particoesEmX
         permODA(1,cont)=permODATemp(1,cont2)
         permODA(2,cont)=permODATemp(2,cont2)
         if((permODA(1,cont).ne.0d0).or.(permODA(2,cont).ne.0d0)) write(252,'(i5, 2e15.5)') cont, permODA(:,cont)
         cont=cont+1
         cont2=cont2+1
      end do      
   end do
   deallocate(permODATemp)
      
   ! do i=1, numParticoes
   !    print*, i , permODA(1,i), permODA(2,i)
   ! end do      
      
    
   ! Usado para localizar o centro do elemento
   if(nen==3) then
      call shlt(shlCentroElem,wCentroElem,1,nen)
   else if(nen==4.and.nsd==2) then
      call shlq(shlCentroElem,wCentroElem,1,nen)
   else if(nen==4.and.nsd==3) then
      call shlt3D(shlCentroElem,wCentroElem,1,nen)
   else if(nen==8) then
      call shlq3d(shlCentroElem,wCentroElem,1,nen)
   end if  
    
   DO NEL=1,NUMEL
      !....    LOCALIZE UNKNOWNS AND COORDINATES
      CALL LOCAL(conecNodaisElem(1,NEL),x,xl,NEN,NSD,NSD)
      CALL LOCAL(conecNodaisElem(:,nel), normal, normalNoLocal, nen, nsd, nsd)

      do i=1,nsd
         xCentroElem(i) = dot_product(xl(i,:), shlCentroElem(nrowsh,:,1))
      end do
         
      ! Se for elemento de Matriz (rocha intacta):
      ! Pesquisa se o centro do elemento pertence a janela, e se sim  atualiza a permeabilidade
      if(mat(nel)==1) then  !no Arapua é 28
         do i=1, numParticoes
            if ((xCentroElem(1)>=ptsIntervalos(1,i)).and.(xCentroElem(1)<=ptsIntervalos(3,i))) then
               if ((xCentroElem(2)>=ptsIntervalos(2,i)).and.(xCentroElem(2)<=ptsIntervalos(4,i))) then
                  campoPermeabilidade(1,nel)=campoPermeabilidade(1,nel)+permODA(1,i)
                  campoPermeabilidade(2,nel)=campoPermeabilidade(2,nel)+permODA(2,i)
               endif
            endif
         end do
      endif         
   ENDDO
     
end subroutine

!==========================================================================================================    
subroutine mapearPermeabilidadeODA3D(conecNodaisElem, numel, numnp, nsd, nen)
   use mMalha,            only: x
   use mGlobaisEscalares, only: nrowsh
   use mGlobaisArranjos,  only: c, mat, tiposElemBordas, campoPermeabilidade
   use mMalha,            only: local, normal, normalElemBorda, numelBordas, nenFratura, conecNodaisBordas
   use mFuncoesDeForma,   only: shlt, shlq, shlq3d, shlt3D
   use mPotencial,        only: estrutSistEqP
    
   implicit none
    
   integer, intent(in) :: numel, nsd, nen, numnp
   integer, intent(in) :: conecNodaisElem(nen, numel)
    
   real*8 :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1)
   real*8 :: normalNoLocal(nsd,nen)
   real*8 :: xl(nsd,nen), domX, domY, domZ
   integer :: numParticoes, tamIntervalos, i, j, nel!, m
   real*8 :: pontoInicialX, pontoInicialY!, pontoInicialZ
   real*8 :: px, py!, pz
   integer :: particoesEmX, particoesEmY, particoesEmZ!, numnpParticoes
   integer :: cont, cont2, tamanhoDaMatriz
   real*8, allocatable :: ptsIntervalos(:,:), permODA(:,:,:), permODATemp(:,:,:)
   character(len=24) :: nomeArq
   integer :: celula
   real*8 :: menor, maior
      
   tamIntervalos=10 !2 para Arapua, 10 para cenario 3
   ! tamIntervalos=2 
      
   ! Calculando o tamanho do dominio
   domX=maxval(x(1,:))-minval(x(1,:))
   domY=maxval(x(2,:))-minval(x(2,:))
   if(nsd==3)domZ=maxval(x(3,:))-minval(x(3,:))
   print*, "Dominio", domX, domY, domZ
      
   ! Calculando o número de partições em cada direção
   particoesEmX=int(domX/tamIntervalos)
   particoesEmY=int(domY/tamIntervalos)
   if(nsd==3)particoesEmZ=int(domZ/tamIntervalos)
   print*, particoesEmX,particoesEmY !100, 87
      
   ! Número total de partições
   numParticoes=particoesEmX*particoesEmY
   print*, "numParticoes", numParticoes

   ! Define os pontos das janelas
   allocate(ptsIntervalos(4,numParticoes))
   cont=0
   pontoInicialX=minval(x(1,:))
   pontoInicialY=minval(x(2,:))
   do j=1, particoesEmY
      py=pontoInicialY+(j-1)*tamIntervalos
      do i=1, particoesEmX
         cont=cont+1
         px=pontoInicialX+(i-1)*tamIntervalos
         ptsIntervalos(1,cont)=px
         ptsIntervalos(2,cont)=py
         ptsIntervalos(3,cont)=px+tamIntervalos
         ptsIntervalos(4,cont)=py+tamIntervalos
      end do
   end do

   ! Leitura do arquivo de permeabilidades ODA vindos do Matlab
   tamanhoDaMatriz=numParticoes*2
   allocate(permODATemp(20,2,tamanhoDaMatriz))
      
   print*, "tamanhoDaMatriz", tamanhoDaMatriz
      
   do celula=1, 20
      call gerarLabelArqs(nomeArq, celula)
         
      open(75,file=nomeArq)
      do j=1, 2
         read(75,*) (permODATemp(celula,j,i),i=1,tamanhoDaMatriz)
      end do
      close(75)
   end do
      
   print*, "depois das leituras"
   ! Separação das permeabilidades de intersse, Kx e Ky, descarte de Kxy e Kyx
   allocate(permODA(20,2,numParticoes))
      
   do celula=1, 20
      cont =1
      do i=1, tamanhoDaMatriz, 2
         permODA(celula,1,cont)=permODATemp(celula,1,i) 
         permODA(celula,2,cont)=permODATemp(celula,2,i+1) 
         cont =cont +1
      enddo     
   end do
   deallocate(permODATemp)

   ! Ordenação. No matlab numera de cima para baixo e nós de baixo para cima.
   allocate(permODATemp(20,2,numParticoes))
   permODATemp=permODA
      
   do celula=1,20
      cont=1
      do j=1, particoesEmY
         cont2=numParticoes-(particoesEmX*j)+1
         do i=1, particoesEmX
            permODA(celula,1,cont)=permODATemp(celula,1,cont2)
            permODA(celula,2,cont)=permODATemp(celula,2,cont2)
            if((permODA(celula,1,cont).ne.0d0).or.(permODA(celula,2,cont).ne.0d0)) then
               write(252,'(i5, 2e15.5)') cont, permODA(celula,:,cont)
            endif
            cont=cont+1
            cont2=cont2+1
         end do      
      end do
   enddo
   deallocate(permODATemp)
      
   ! Usado para localizar o centro do elemento
   if(nen==3) then
      call shlt(shlCentroElem,wCentroElem,1,nen)
   else if(nen==4.and.nsd==2) then
      call shlq(shlCentroElem,wCentroElem,1,nen)
   else if(nen==4.and.nsd==3) then
      call shlt3D(shlCentroElem,wCentroElem,1,nen)
   else if(nen==8) then
      call shlq3d(shlCentroElem,wCentroElem,1,nen)
   end if  
    
   DO NEL=1,NUMEL
      !....    LOCALIZE UNKNOWNS AND COORDINATES
      CALL LOCAL(conecNodaisElem(1,NEL),x,xl,NEN,NSD,NSD)
      CALL LOCAL(conecNodaisElem(:,nel), normal, normalNoLocal, nen, nsd, nsd)

      do i=1,nsd
         xCentroElem(i) = dot_product(xl(i,:), shlCentroElem(nrowsh,:,1))
      end do
         
      ! Se for elemento de Matriz (rocha intacta):
      ! Pesquisa se o centro do elemento pertence a janela, e se sim  atualiza a permeabilidade
      if(mat(nel)==1) then  !no Arapua é 28
         do i=1, numParticoes
            if ((xCentroElem(1)>=ptsIntervalos(1,i)).and.(xCentroElem(1)<=ptsIntervalos(3,i))) then
               if ((xCentroElem(2)>=ptsIntervalos(2,i)).and.(xCentroElem(2)<=ptsIntervalos(4,i))) then
                  do celula=1, 20
                     menor=(celula-1)*1.d0
                     maior=menor+1.d0
                     if((xCentroElem(3)>=menor).and.(xCentroElem(3)<maior)) then
                        campoPermeabilidade(1,nel)=campoPermeabilidade(1,nel)+permODA(celula,1,i)
                        campoPermeabilidade(2,nel)=campoPermeabilidade(2,nel)+permODA(celula,2,i)
                     endif
                  end do
               endif
            endif
         end do
      endif
   ENDDO

end subroutine

!==========================================================================================================    
subroutine gerarLabelArqs(label, iteracao)

   implicit none

   character(LEN=24), intent(out) :: label
   integer, intent(in) :: iteracao

   character(LEN=30) :: labelAux, num     
   character(LEN=4) :: extensao
   integer :: i
            
   write(num,'(i5)') iteracao
   ! labelAux="t="//ADJUSTL(num)
   labelAux="permeabilidadeODA_"//ADJUSTL(num)
   extensao=".inc"
   label = trim(labelAux)//".inc"

end subroutine     

!==========================================================================================================    
subroutine mapearMaterialCarste3D(conecNodaisElem, numel, numnp, nsd, nen)
   use mMalha,            only: x
   use mGlobaisEscalares, only: nrowsh
   use mGlobaisArranjos,  only: c, mat, tiposElemBordas, campoPermeabilidade
   use mMalha,            only: local, normal, normalElemBorda, numelBordas, nenFratura, conecNodaisBordas
   use mFuncoesDeForma,   only: shlt, shlq, shlq3d, shlt3D
   use mPotencial,        only: estrutSistEqP
    
   implicit none
    
   integer, intent(in) :: numel, numnp, nsd, nen
   integer, intent(in) :: conecNodaisElem(nen, numel)
    
   real*8 :: xCentroElem(nsd), shlCentroElem(nrowsh,nen,1), wCentroElem(1)
   real*8 :: normalNoLocal(nsd,nen)
   real*8 :: xl(nsd,nen), domX, domY, domZ
   integer :: i, j, nel, m
   real*8 :: tamBrechaCaotica, tamBrechaMosaica, tamBrechaCraquelada
   real*8 :: intervaloBrechaCaotica(2), intervaloBrechaMosaica(2), intervaloBrechaCraquelada(2)
    
      
   ! Calculando o tamanho do dominio
   domX=maxval(x(1,:))-minval(x(1,:))
   domY=maxval(x(2,:))-minval(x(2,:))
   domZ=maxval(x(3,:))-minval(x(3,:))
   print*, "Dominio", domX, domY, domZ
      
   tamBrechaCaotica   =(domZ*40)/100 !regra de 3
   tamBrechaMosaica   =(domZ*30)/100 !regra de 3
   tamBrechaCraquelada=(domZ*30)/100 !regra de 3
      
   intervaloBrechaCaotica(1)=minval(x(3,:))
   intervaloBrechaCaotica(2)=intervaloBrechaCaotica(1)+tamBrechaCaotica

   intervaloBrechaMosaica(1)=intervaloBrechaCaotica(2)
   intervaloBrechaMosaica(2)=intervaloBrechaMosaica(1)+tamBrechaMosaica
      
   intervaloBrechaCraquelada(1)=intervaloBrechaMosaica(2)
   intervaloBrechaCraquelada(2)=intervaloBrechaCraquelada(1)+tamBrechaCraquelada

   ! print*, intervaloBrechaCaotica
   ! print*, intervaloBrechaMosaica
   ! print*, intervaloBrechaCraquelada
   ! stop

   ! Usado para localizar o centro do elemento
   if(nen==3) then
      call shlt(shlCentroElem,wCentroElem,1,nen)
   else if(nen==4.and.nsd==2) then
      call shlq(shlCentroElem,wCentroElem,1,nen)
   else if(nen==4.and.nsd==3) then
      call shlt3D(shlCentroElem,wCentroElem,1,nen)
   else if(nen==8) then
      call shlq3d(shlCentroElem,wCentroElem,1,nen)
   end if  
    
   DO NEL=1,NUMEL
      !....    LOCALIZE UNKNOWNS AND COORDINATES
      CALL LOCAL(conecNodaisElem(1,NEL),x,xl,NEN,NSD,NSD)
      CALL LOCAL(conecNodaisElem(:,nel), normal, normalNoLocal, nen, nsd, nsd)

      do i=1,nsd
         xCentroElem(i) = dot_product(xl(i,:), shlCentroElem(nrowsh,:,1))
      end do
         
      ! Se for elemento de Carste
      ! Pesquisa se o centro do elemento pertence ao intervalo da região do carste, se sim, altera o tipo de material
      if(mat(nel)==2) then
         if ((xCentroElem(3)>=intervaloBrechaCaotica(1)).and.(xCentroElem(3)<=intervaloBrechaCaotica(2))) then
            mat(nel)=2
         endif
         if ((xCentroElem(3)>=intervaloBrechaMosaica(1)).and.(xCentroElem(3)<=intervaloBrechaMosaica(2))) then
            mat(nel)=3
         endif
         if ((xCentroElem(3)>=intervaloBrechaCraquelada(1)).and.(xCentroElem(3)<=intervaloBrechaCraquelada(2))) then
            mat(nel)=4
         endif
         m=mat(nel)
         campoPermeabilidade(1,nel) = c(3,m)
         campoPermeabilidade(2,nel) = c(4,m)
         campoPermeabilidade(3,nel) = c(5,m)
      endif
   ENDDO

end subroutine

!==========================================================================================================    
subroutine processamento()
   use mLeituraEscrita,    only: ignuplotPotencial, ignuplotFluxo, iparaviewPotencial, iparaviewFluxo
   use mLeituraEscrita,    only: escreverArquivosSaida_Potencial, escreverArquivosSaida_Fluxo
   use mLeituraEscrita,    only: escreverValoresReferencia, escreverArquivoParaview

   use mGlobaisEscalares,  only: mpi_comm, transiente, tempo, passoTempo, simulandoReferencia, pi
   use mGlobaisEscalares,  only: tempoTotalSimulacao, numPassosTempo, deltaT, relax_param, viscosidade
   use mGlobaisEscalares,  only: max_iter_pressao, tolerancia_relativa, saltoPressao
   use mGlobaisEscalares,  only: flagTransmib, corrTempo
   use mGlobaisArranjos,   only: matFratura, tiposElemBordas, c, campoPermeabilidade
   use mGlobaisArranjos,   only: vazao

   use mUtilSistemaEquacoes,  only: montarEstrutDadosSistEqAlq

   use mMalha,    only: x, numel, nsd, nen, numnp
   use mMalha,    only: nenFratura, numnpFratura, numelFratura
   use mMalha,    only: numConexoesPorElem, numelBordas
   use mMalha,    only: conecNodaisBordas, conecNodaisElem, conecNodaisFratura

   use mPotencial,   only: estrutSistEqP
   use mPotencial,   only: montarSistEqAlgPotencial,estrutSistEqP

   use mFluxo,    only: estrutSistEqF, estrutSistEqFluxoNormal
   use mFluxo,    only: calcularVazao, calcularPropsEfetivas, velocidadeBordas
   use mFluxo,    only: montarSistEqAlgVelocidadeNormal
   use mFluxo,    only: calcularVelocidade,calcularVazaoNew

   use mTransmissividade, only: transmissivity

   use mSolverHypre,    only: inicializarMPI, finalizarMPI
   use mSolverHypre,    only: myid, num_procs

   use mVisualizacaoHdf5, only: hdf5_dataset_integer, hdf5_dataset_real, hdf5_dataset_double
   use mVisualizacaoHdf5, only: name_struct, p_tree_hdf5, p_root_hdf5
   use mVisualizacaoHdf5, only: hdf5_structure_addPressure, hdf5_structure_addVelocity
   use mVisualizacaoHdf5, only: hdf5_structure_addPressureFalha, hdf5_structure_addVelocFalha
   use mVisualizacaoHdf5, only: xdmf_structure_bloco, xdmf_structure_falha

   implicit none

   !---------------------------------------------------------------------------------
   !.... solution driver program
      
   real*8 :: t1, t2
   integer :: i, dof, j, k!, iteracaoAtual, no, noBordas, tipoBorda
   logical :: simetria
   character(LEN=12) :: label
   character(len=8) :: passoTempoStr
   integer, parameter :: ioutPotencial=21, ioutFluxo=22
   real*8 :: Lx, Ly, Lz
   logical :: alterouDeltaT
   integer :: nTempo
   real*8, allocatable :: time_array(:)
   real*8 :: vazLin(nsd)
   real*8 :: trans, vazBord(2)

   estrutSistEqP%u=estrutSistEqP%uInicial
   if(transiente) estrutSistEqP%uTempoAnt=estrutSistEqP%uInicial

   ! //TODO corrigir o calculo da transmissibilidade para funcionar na análise transiente.
   flagTransmib=.true.
   if(transiente) flagTransmib=.false.

   !---------------------------------------------------------------------------------
   !.... initialization phase
   !.... Determinacao para analise estatica ou transiente
   tempo = 0.d0
   passoTempo = 0
   numPassosTempo=1
   alterouDeltaT=.false.
   
   if(transiente.eqv..true.)then
      numPassosTempo=floor(tempoTotalSimulacao/deltaT)
      print*, "DELTAT=", deltaT
      ! print*, "numPassosTempo=", numPassosTempo

      ! Armazenar produção do poco 
      ! print*, "allocating space para o vazaoAnt" !XXXXX
      ! allocate(prodAcum(nsd*2+1+2))
      ! allocate(vazaoAnt(nsd*2+1+2))
      ! prodAcum(:) = 0.d0
      ! vazaoAnt(:) = 0.d0
   endif

   ! escrita da serie temporal no arquivo hdf5
   if(numPassosTempo.eq.1)then
      nTempo = 2
   else
      nTempo = numPassosTempo+1
   endif
   allocate(time_array(nTempo))
   do i=1,nTempo
      time_array(i)=tempo+(deltaT/corrTempo)*dble(i-1)
   end do
   time_array(nTempo)=tempoTotalSimulacao/corrTempo
   call hdf5_dataset_double("TimeSerie", time_array, 1, nTempo, p_root_hdf5(4))
   deallocate(time_array)

   ! Impressao da variavel
   ! write(*,'("NumPasso(",i5,")")') numPassosTempo
   ! write(*,'("Tf(",E8.3,")")') tempoTotalSimulacao
   ! write(*,'("Dt(",E8.3,")")') deltaT

   !---------------------------------------------------------------------------------
   !....  Montagen das estruturas de matrizes
   numConexoesPorElem=nen !  usado para criar lista de vizinhos dos nOs
   if (estrutSistEqP%optSolver=='pardiso') then
      if(nsd==2)estrutSistEqP%numCoefPorLinha=27    ! usado para criar LMStencilEq
      if(nsd==3)estrutSistEqP%numCoefPorLinha=81
   end if
   call montarEstrutDadosSistEqAlq(estrutSistEqP)
   ! Para o fluxo consistente
   if (flagTransmib) call montarEstrutDadosSistEqAlq(estrutSistEqFluxoNormal)
      
   ! !    Limpa matriz
   ! do i=1, estrutSistEqP%nalhs
   !    estrutSistEqP%alhs(i) = 0.0
   ! end do

   ! do i=1, estrutSistEqF%nalhs
   !    estrutSistEqF%alhs(i) = 0.0
   ! end do

   ! Impressao da variavel
   ! do j=1,estrutSistEqP%nlvect
   !    do i=1,numnp
   !       write(*,'("force(",i5,")=")',advance='no') i
   !       do k=1,estrutSistEqP%ndof
   !          write(*,'(E11.3)',advance='no') estrutSistEqP%f(j,k,i)
   !       enddo
   !       WRITE(*,*)
   !    enddo
   !    write (*,*)
   ! enddo
   ! do i=1, estrutSistEqP%nalhs
   !    write(*,'("estrutSistEqP(",i5,")=",E11.3)') i, estrutSistEqP%alhs(i)
   ! end do

   ! stop

   

   !---------------------------------------------------------------------------------
   !....  Loop para solucao para cada incremento de tempo. Caso estatico apenas uma iteracao
   do passoTempo=1, numPassosTempo
      ! Teste para saida do loop       
      if(tempo>tempoTotalSimulacao) return

      ! Determina o tempo analisado
      if(transiente.eqv..true.) then
         tempo = tempo + deltat
         write(*,*) "" ! Newline
         write(*,*) "Tempo", tempo, "Passo de tempo", passoTempo
      else
         tempo = tempoTotalSimulacao
      endif

      ! Limpa vetor
      estrutSistEqP%brhs(:) = 0d0
      ! do i=1,estrutSistEqP%neq
      !    estrutSistEqP%brhs(i) = 0d0
      ! end do

      ! Fluxo consistente
      if (flagTransmib) estrutSistEqFluxoNormal%brhs(:)= 0d0

      ! Montagem das equacoes
      call timing(t1)
      print*, "Calculando Pressao"
      call montarSistEqAlgPotencial(estrutSistEqP)
      call timing(t2)
      print*, "Tempo de montagem do sistema de equacoes", t2-t1

      ! Escrevendo o vetores iniciais no arquivo hdf5.
      if (passoTempo.eq.1) then
         write(passoTempoStr,'(i0.8)') 0
         call hdf5_structure_addPressure(passoTempoStr, p_tree_hdf5(7))
         call hdf5_structure_addVelocity(passoTempoStr, p_tree_hdf5(8))
         if (numelFratura>0) then
            call hdf5_structure_addPressureFalha(passoTempoStr, p_tree_hdf5(11))
            call hdf5_structure_addVelocFalha(passoTempoStr, p_tree_hdf5(12))
         endif
      endif
      ! Impressao da variavel (Matriz de forca)
      ! do j=1,estrutSistEqP%nlvect 
      !    do i=1,numnp
      !       WRITE(*,'("force(",i5,")=")',advance='no') i
      !       do k=1,estrutSistEqP%ndof
      !          write(*,'(E10.4)',advance='no') estrutSistEqP%f(j,k,i)
      !       enddo
      !       WRITE(*,*)
      !    enddo
      !    write (*,*)
      ! enddo
      ! Impressao da variavel (Matriz de ID)
      ! do i=1,numnp 
      !    WRITE(*,'("id(",i5,")=")',advance='no') i
      !    do k=1,estrutSistEqP%ndof
      !       write(*,'(i5)',advance='no') estrutSistEqP%id(k,i)
      !    enddo
      !    WRITE(*,*)
      ! enddo
      ! Impressao da variavel (Matriz brhs)
      ! do i=1,estrutSistEqP%neq 
      !    write(*,'("brhs(",i5,"->",E10.4,")")') i, estrutSistEqP%brhs(i)
      ! enddo
      ! Impressao da variavel (Vetor u)
      ! do i=1,numnp 
      !    WRITE(*,'("u(",i5,")")',advance='no') i
      !    do k=1,estrutSistEqP%ndof
      !       write(*,'(E12.4,"|")',advance='no') estrutSistEqP%u(k,i)
      !    enddo
      !    WRITE(*,*)
      ! enddo

      ! Resolvendo o problema (K u = F)
      call timing(t1)
      label='potencial'
      print*, "Resolvendo pressao"
      print*, "Numéro de Graus de Liberdade ->", estrutSistEqP%neq
      call solver(estrutSistEqP, label)
      call timing(t2)
      print*, "Tempo de solver:", t2-t1
      print*, "pressoes", maxval(estrutSistEqP%u(1,:)),minval(estrutSistEqP%u(1,:))
      ! Impressao da variavel (Vetor u)
      ! do i=1,numnp
      !    WRITE(*,'("u(",i5,")")',advance='no') i
      !    do k=1,estrutSistEqP%ndof
      !       write(*,'(E12.4,"|")',advance='no') estrutSistEqP%u(k,i)
      !    enddo
      !    WRITE(*,*)
      ! enddo

      ! Para todos os nós os que nao sao contorno e de falha, copia o 1 dof para os demais.
      if( (numelFratura>0) .AND. (saltoPressao .EQV. .true.) ) then    
         ! Para que o calculo do salto faca sentido temos que atribuir valores
         ! para os graus de liberdade adicionais das extremidades da fratura
         do i=1,numnp
            do dof=2,estrutSistEqP%ndof
               if (estrutSistEqP%id(dof,i) == 0) then
                  estrutSistEqP%u(dof,i) = estrutSistEqP%u(1,i)
               end if
            end do
         end do
      endif
      write(passoTempoStr,'(i0.8)') passoTempo
      call hdf5_structure_addPressure(passoTempoStr, p_tree_hdf5(7))
      if (numelFratura>0) then
         call hdf5_structure_addPressureFalha(passoTempoStr, p_tree_hdf5(11))
      endif

      
      ! Determinação da Velocidade pelo pos-processamento
      ! Comentado para tirar velocidade - pag 100 hughes.
      call timing(t1)
      print*, "Resolvendo Velocidade"
      call calcularVelocidade() 
      call timing(t2)
      print*, "Tempo para o computo da velocidade:", t2-t1
      
      ! escrita do arquivo hdf5
      call hdf5_structure_addVelocity(passoTempoStr, p_tree_hdf5(8))
      ! Possivel inclusão de rotina com dependencia cíclica.
      ! //TODO : Corrigir estrutura hdf5_structure_addVelocFalha
      ! if (numelFratura>0) then
      !    !    call hdf5_structure_addVelocBlocoFalha("VM", vm_mediaF, (estrutSistEqP%ndof-2), p_tree_hdf5(12))
      !    call hdf5_structure_addVelocFalha("1", p_tree_hdf5(12))
      ! endif
      
      ! -------------------------------------------
      ! Trabalhando Aqui!!!!!!!!
      ! computo velocidade e vazao usando a formulação do artigo do Gresho
      if (flagTransmib) then
         call montarSistEqAlgVelocidadeNormal()
         label='velocidade'
         print*, "Resolvendo velocidade normal nas bordas Dirichlet"!             
         call solver(estrutSistEqFluxoNormal, label)
         print*, "Calcular vazao"             
         call calcularVazaoNew(estrutSistEqP,estrutSistEqFluxoNormal, Lx, Ly, Lz, passoTempo)
         print*, "---separação entre metodos de determinação da vazão---"
         vazBord(1)=abs(vazao(2))
         vazBord(2)=abs(vazao(4))
      endif
      ! -------------------------------------------

      print*, "Calculando a Vazao e props efetivas"         
      call calcularVazao(estrutSistEqP, Lx, Ly, Lz)      
      call calcularPropsEfetivas()

      ! -------------------------------------------Trabalhando Aqui!!!!!!!!
      ! computo da transmissividade
      if (flagTransmib) call transmissivity(trans,vazBord)
      ! print*, "-------------Transmissividade-------------"
      ! print*, "Transmissividade ->", trans
      ! do i=1, trans_nelm
      !    WRITE(*,*) "transLM(", trans_lm(i,1), ",", trans_lm(i,2), ")"
      ! end do

      do i=1,nsd
         vazLin(i)=0.
      enddo
      ! LM -> conect
      ! veloc -> velocidades nos nós
      ! nelm -> numero de elementos na linha
      ! call vazaoLinha(vazLin, x, LM, veloc, nelm, nsd, numnp)
      
      !salva iteração atual
      estrutSistEqP%uTempoAnt = estrutSistEqP%u  
      
   end do !tempo

   ! Criação do arquivo de vizualisação
   call xdmf_structure_bloco(numel, numnp, nen, nsd, numPassosTempo+1)
   call xdmf_structure_falha(numelFratura, numnpFratura, nenFratura, nsd, 2)
   return

   contains

   !==========================================================================================================
   subroutine solver(estrutSistEq_, label_)         
      use mSolverGaussSkyline,      only: solverGaussSkyline, FACTOR, BACK
      use mSolverPardiso,           only: solverPardisoEsparso
      use mSolverHypre
      use mUtilSistemaEquacoes,     only: btod
      use mMalha,                   only: x, numnp
      use mEstruturasDadosSistEq,   only: estruturasArmazenamentoSistemaEq

      implicit none
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*), intent(in) :: label_
      
      character(LEN=12) :: etapa
      integer * 4 ::  num_iterations, printSol = 2
      real*8 ::  final_res_norm, elapsedT, tol, solutionNorm
      integer*4   :: equacao!, i, dof

      ! write(*,*) 'solver ', optSolver_ ,', ',  label_

      if (estrutSistEq_%optSolver=='skyline') then
         ! call solverGaussSkyline(estrutSistEq_)
         IF(passoTempo==1)call factor(estrutSistEq_%alhs,estrutSistEq_%idiag,estrutSistEq_%nalhs,estrutSistEq_%neq)            
         call back  (estrutSistEq_%alhs,estrutSistEq_%brhs,estrutSistEq_%idiag,estrutSistEq_%neq)

         ! deallocate(estrutSistEq_%alhs) ! pode ser desnecessario para problemas
      endif
          
      if (estrutSistEq_%optSolver=='pardiso') then
         simetria=.true.
         ! if(PASSOTEMPO==1) THEN
         !    etapa='full';
         ! else
         !    etapa='back';
         ! endif
         etapa='full';
         call solverPardisoEsparso(estrutSistEq_, simetria, etapa, label_)
      endif
          
      if(estrutSistEq_%optSolver=='hypre') then
         write(*,'(2a)') ' iterativo ', estrutSistEq_%optSolver
            
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

   !==========================================================================================================    
   subroutine escreverSistema_MTX(optSolver_, estrutSistEq_, nomeArq_)
      use mMalha,                   only: numnp, numel, nen
      use mSolverPardiso,           only: escreverSistemaAlgCSRemMTX
      use mSolverGaussSkyline,      only: escreverSistemaSkylineEmMTX
      use mEstruturasDadosSistEq,   only: estruturasArmazenamentoSistemaEq 

      character(LEN=*), intent(in)  :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*) :: nomeArq_
      integer*4 :: nee

      nee = estrutSistEq_%ndof * nen
            
      if (optSolver_=='skyline') then
         call escreverSistemaSkylineEmMTX (estrutSistEq_%alhs, estrutSistEq_%brhs, estrutSistEq_%idiag, &
                  estrutSistEq_%lm,estrutSistEq_%id, conecNodaisElem,nee,nen,numel,numnp,&
                  estrutSistEq_%neq, nomeArq_)
      end if


      if (optSolver_=='PardisoEsparso') then
         call escreverSistemaAlgCSRemMTX    (estrutSistEq_%Alhs, estrutSistEq_%brhs, estrutSistEq_%Ap, estrutSistEq_%Ai, &
         estrutSistEq_%nAlhs, estrutSistEq_%neq, nomeArq_)
      end if
   end subroutine escreverSistema_MTX

   !==========================================================================================================    
   subroutine escreverSolSistema_MTX(optSolver_, estrutSistEq_, nomeArq_)
      use mMalha,                   only: numnp, numel, nen
      use mSolverPardiso,           only: escreverSistemaAlgCSRemMTX
      use mSolverGaussSkyline,      only: escreverSistemaSkylineEmMTX
      use mEstruturasDadosSistEq,   only: estruturasArmazenamentoSistemaEq

      character(LEN=*), intent(in)  :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*) :: nomeArq_
      integer*4 :: nee

      call escreverBRHS (estrutSistEq_%brhs, estrutSistEq_%neq, nomeArq_)

   end subroutine escreverSolSistema_MTX

   !==========================================================================================================    
   subroutine escreverBRHS (brhs_, neq_, nomeArq_)
      real*8, pointer :: brhs_(:)
      integer :: neq_
      character (LEN=*) :: nomeArq_

      integer :: ii

      open (unit=1111, file=trim(nomeArq_))
      write(1111, *) 1, neq_
      do ii = 1, neq_
         write(1111, *) ii, brhs_(ii)
      end do
      close(1111)
   end subroutine
end subroutine processamento

!==========================================================================================================    
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     
!==========================================================================================================    
subroutine leituraParamNumericosPropFisica ()
   use mGlobaisEscalares, only: iprtin
   use mGlobaisEscalares, only: numParElem, numProps, nrowsh, numat, npint, nicode
   use mGlobaisArranjos,  only: npar,c, mat, grav
   use mMalha,            only: numel,nsd,numnp,nen,conecNodaisElem
   use mLeituraEscrita,   only: iin, iecho

   !.... program to read, generate and write element data
   implicit none

   integer*4 :: m, n, i
   character(len=80) :: formatoLeitura

   nrowsh = 3
   if (nsd==3) nrowsh=nrowsh+1

   allocate(npar(numParElem))
   allocate(grav(3))

   formatoLeitura='(16I10)'
   read(iin, formatoLeitura) (npar(i),i=1,numParElem)

   nicode = npint
   if (nicode.eq.0) nicode=nen

   numat  = npar( 1)
   allocate(c(7,numat))

   write(iecho,1000) numel,numat,nen,npint

   ! read material properties
   do 400 n=1,numat
      if (mod(n,50).eq.1) write(iecho,4000) numat
      read(iin,5000) m,(c(i,m),i=1,numProps)
      write(iecho,6000) m,(c(i,m),i=1,numProps)
   400 continue

   ! constant body forces
   read (iin,7000) (grav(i),i=1,3)
   write (iecho,8000) (grav(i),i=1,3)

   return

   1000 format(//,&
            ' two/three-n o d e    e l e m e n t s ',//,5x,&
            ' number of elements  . . . . . . . . . . . (numel ) = ',i8,//5x,&
            ' number of element material sets . . . . . (numat ) = ',i5,//5x,&
            ' number of element nodes . . . . . . . . . (nen   ) = ',i5,//5x,&
            ' number of integration points. . . . . . . (npint  ) = ',i5)
   ! 2000  format(16i5)
   4000  format(///,&
            ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
            ' number of material sets . . . . . . . . . . (numat ) = ', i5///,&
            2x,'set',4x,'Kx ', 10x,'Ky',10x,'Kz')
   5000  format(i10,6f10.0)
   6000  format(2x,i3,1x,6(1x,1pe14.4))
   7000 format(8f10.0)
   8000 format(///,&
            ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
            ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
            ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
            ' exemplo 3............................ = ',      1pe15.8,//)
   ! 9000  format(i5)
end subroutine leituraParamNumericosPropFisica

!==========================================================================================================    
subroutine readSetupPhaseDS
   !> Efetua a leitura completa dos dados da etapa de setup.
   use mInputReader,      only: readStringKeywordValue, readIntegerKeywordValue, readRealKeywordValue, readLogicalKeywordValue
   use mGlobaisArranjos,  only: title, bordaPrescrita
   use mGlobaisEscalares
   use mFuncoesDeForma,   only: gerarPtosIntegracaoLaterais2D
   use mMalha,            only: nsd, numnp, numnpFratura, numel, nen
   use mMalha,            only: numelFratura, nenFratura, calcularNenFratura, numelBordas
   use mPotencial,        only: estrutSistEqP
   use mFluxo,            only: estrutSistEqF, estrutSistEqFluxoNormal
   use mLeituraEscrita,   only: tipo_arq_saida, passosPorImpressao
   use mMLayer,           only: nTypesMLayers
   use mTransmissividade, only: trans_nelm, trans_limite, trans_lm, Material_Ignore

   implicit none
   character(len=50) keyword_name
   character(len=20) auxTempo
   integer temp, i

   numPassosTempo=1

   keyword_name = "title"
   call readStringKeywordValue(keyword_name, title, 'simulacao sem titulo')

   keyword_name = "exec"
   call readIntegerKeywordValue(keyword_name, exec, 0_4)

   keyword_name = "iprtin"
   call readIntegerKeywordValue(keyword_name, iprtin, 0_4)

   keyword_name = "tipo_arq_saida"
   call readIntegerKeywordValue(keyword_name, tipo_arq_saida, 0)

   keyword_name = "nsd"
   call readIntegerKeywordValue(keyword_name, nsd, 0_4)

   keyword_name = "ndof"
   call readIntegerKeywordValue(keyword_name, estrutSistEqP%ndof, 3)

   keyword_name = "numnp"
   call readIntegerKeywordValue(keyword_name, numnp, 0_4)

   keyword_name = "numnp_fratura"
   call readIntegerKeywordValue(keyword_name, numnpFratura, 0_4)

   keyword_name = "numel"
   call readIntegerKeywordValue(keyword_name, numel, 0_4)

   keyword_name = "numel_fratura"
   call readIntegerKeywordValue(keyword_name, numelFratura, 0_4)

   keyword_name = "numel_bordas"
   call readIntegerKeywordValue(keyword_name, numelBordas, 0_4)


   keyword_name = "transmissividade_nelm"
   call readIntegerKeywordValue(keyword_name, trans_nelm, 0_4)

   keyword_name = "transmissividade_limite"
   call readRealKeywordValue(keyword_name, trans_limite, 0d0)

   keyword_name = "transmissividade_mat_ign_nelm"
   call readIntegerKeywordValue(keyword_name, Material_Ignore%nelm, 0)

   ! keyword_name = "numel_bordas_dirichlet"
   ! call readIntegerKeywordValue(keyword_name, numelBordasDirichlet, 0_4)

   ! keyword_name = "numelDescontinuidadeNaBorda"
   ! call readIntegerKeywordValue(keyword_name, numelDescontinuidadeNaBorda, 0)
   ! print*, "numelDescontinuidadeNaBorda", numelDescontinuidadeNaBorda

   keyword_name = "num_type_desc"
   call readIntegerKeywordValue(keyword_name, nTypesMLayers, 0_4)

   keyword_name = "nen"
   call readIntegerKeywordValue(keyword_name, nen, 0_4)
   call calcularNenFratura

   keyword_name = "npint"
   call readIntegerKeywordValue(keyword_name, npint, 0_4)
   call gerarPtosIntegracaoLaterais2D

   keyword_name = "nlvectP"
   call readIntegerKeywordValue(keyword_name, estrutSistEqP%nlvect, 0_4)

   keyword_name = "nlvectF"
   call readIntegerKeywordValue(keyword_name, estrutSistEqF%nlvect, 0_4)

   keyword_name = "deltaT"
   call readRealKeywordValue(keyword_name, deltaT, 0d0)
   if(deltaT>0.d0) then
      transiente=.true.
      deltaT=deltaT*corrTempo
   endif

   keyword_name = "passos_por_impressao"
   call readIntegerKeywordValue(keyword_name, passosPorImpressao, 1_4)

   keyword_name = "tempoTotalSimulacao"
   call readRealKeywordValue(keyword_name, tempoTotalSimulacao, 0d0)
   tempoTotalSimulacao = tempoTotalSimulacao*corrTempo

   keyword_name = "relax_param"
   call readRealKeywordValue(keyword_name, relax_param, 1d0)

   keyword_name = "solver_pressao"
   call readStringKeywordValue(keyword_name, estrutSistEqP%optSolver, 'skyline')

   keyword_name = "solver_velocidade"
   call readStringKeywordValue(keyword_name, estrutSistEqFluxoNormal%optSolver, 'skyline')

   keyword_name = "carregamento"
   call readRealKeywordValue(keyword_name, carregamento, 0d0)

   keyword_name = "carregamento_ref"
   call readRealKeywordValue(keyword_name, carregamento_ref, 0d0)

   ! keyword_name = "pressaoMedia"
   ! call readRealKeywordValue(keyword_name, pressaoMedia, 0d0)
        
   ! keyword_name = "pressaoMedia_Ref"
   ! call readRealKeywordValue(keyword_name, pressaoMedia_Ref, 0d0)

   keyword_name = "pressao_inicial"
   call readRealKeywordValue(keyword_name, estrutSistEqP%uInicial, 1d0)
   simulandoReferencia = .true.

   keyword_name = "perm_y_efetiva_ref"
   call readRealKeywordValue(keyword_name, perm_y_efetiva_ref, 0d0)
   if (perm_y_efetiva_ref > 1d-21) then
      simulandoReferencia = .false.
   end if

   keyword_name = "porosidade_celula_ref"
   call readRealKeywordValue(keyword_name, porosidadeCelula_ref, 0d0)
   if (perm_y_efetiva_ref > 1d-10) then
      simulandoReferencia = .false.
   end if

   keyword_name = "tolerancia_relativa"
   call readRealKeywordValue(keyword_name, tolerancia_relativa, 1d-7)

   keyword_name = "max_iter_pressao"
   call readIntegerKeywordValue(keyword_name, max_iter_pressao, 200)

   keyword_name = "bordaPrescrita_1"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(1), .false.)    
   keyword_name = "bordaPrescrita_2"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(2), .false.)    
   keyword_name = "bordaPrescrita_3"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(3), .false.)    
   keyword_name = "bordaPrescrita_4"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(4), .false.)    
   keyword_name = "bordaPrescrita_5"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(5), .false.)    
   keyword_name = "bordaPrescrita_6"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(6), .false.)    
   keyword_name = "bordaPrescrita_7"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(7), .false.)    
   keyword_name = "bordaPrescrita_11"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(11), .false.)    
   keyword_name = "bordaPrescrita_12"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(12), .false.)    
   keyword_name = "bordaPrescrita_13"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(13), .false.)    
   keyword_name = "bordaPrescrita_14"
   call readLogicalKeywordValue(keyword_name, bordaPrescrita(14), .false.)

end subroutine readSetupPhaseDS 

!==========================================================================================================    
subroutine leituraParamNumericosPropFisicaDS()
   !> Faz a leitura de constant body forces.
   !> @param keyword_name  Keyword especifica para constant body forces.

   use mGlobaisEscalares, only: iprtin, numProps
   use mGlobaisEscalares, only: numParElem, nrowsh, numat, npint, nicode
   use mGlobaisArranjos,  only: npar,c, mat, grav
   use mMalha,            only: numel,nsd,numnp,nen,conecNodaisElem
   use mLeituraEscrita,   only: iecho
   use mInputReader,      only: findKeyword, file_lines

   !.... program to read, generate and write element data
   implicit none

   integer*4 :: m, n, i, keyword_line, nLinhaArqInput
   character(len=80) :: formatoLeitura
   character(len=50) keyword_name

   keyword_name = "nummat"
   keyword_line = findKeyword(keyword_name)
   nLinhaArqInput = keyword_line
   if (keyword_line.eq.-1) return

   nrowsh = 3
   if (nsd==3) nrowsh=nrowsh+1

   allocate(npar(numParElem))
   allocate(grav(3))

   formatoLeitura='(16I10)'
   ! read(iin, formatoLeitura) (npar(i),i=1,numParElem)
   read(file_lines(nLinhaArqInput:), formatoLeitura) (npar(i),i=1,numParElem)
   nLinhaArqInput = nLinhaArqInput + 1

   nicode = npint
   if (nicode.eq.0) nicode=nen

   numat  = npar( 1)
   allocate(c(numProps,numat))

   write(iecho,1000) numel,numat,nen,npint

   ! read material properties
   keyword_name = "prop_fisica_meio"
   keyword_line = findKeyword(keyword_name)
   nLinhaArqInput = keyword_line
   do n=1,numat
      if (mod(n,50).eq.1) write(iecho,4000) numat
      ! c(1,m): porosidade
      ! c(2,m): abertura da fratura (m)
      ! c(3,m): permeabilidade (m^2)
      read(file_lines(nLinhaArqInput:),5000) m,(c(i,m),i=1,numProps)
      ! write(*,5000) m,(c(i,m),i=1,numProps)
      nLinhaArqInput = nLinhaArqInput + 1
      write(iecho,6000) m,(c(i,m),i=1,numProps)
   end do

   !     constant body forces
   keyword_name = "grav"
   keyword_line = findKeyword(keyword_name)
   nLinhaArqInput = keyword_line
   read (file_lines(nLinhaArqInput:),7000) (grav(i),i=1,3)
   nLinhaArqInput = nLinhaArqInput + 1
   write (iecho,8000) (grav(i),i=1,3)

   return

   1000 format(//,&
            ' two/three-n o d e    e l e m e n t s ',//,5x,&
            ' number of elements  . . . . . . . . . . . (numel ) = ',i8,//5x,&
            ' number of element material sets . . . . . (numat ) = ',i5,//5x,&
            ' number of element nodes . . . . . . . . . (nen   ) = ',i5,//5x,&
            ' number of integration points. . . . . . . (npint  ) = ',i5)
   ! 2000 format(16i5)
   4000 format(///,&
            ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
            ' number of material sets . . . . . . . . . . (numat ) = ', i5///,&
            2x,'set',4x,'Kx ', 10x,'Ky',10x,'Kz')
   5000 format(i10,6f15.0)
   6000 format(2x,i3,1x,6(1x,1pe14.4))
   7000 format(8f10.0)
   8000 format(///,&
            ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
            ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
            ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
            ' exemplo 3............................ = ',      1pe15.8,//)
   ! 9000  format(i5)
end subroutine leituraParamNumericosPropFisicaDS

!==========================================================================================================    
subroutine verificarSolver(optSolver_, listaSolverDisponivel_) 

   character(len=7), intent(in) :: optSolver_
   logical, intent(in) :: listaSolverDisponivel_(*)

   if(optSolver_=='pardiso')then
      if (listaSolverDisponivel_(2).eqv..false.) then
         print*, "O Solver escolhido ...,  ", optSolver_,",  não está disponível"
         stop 100
      endif
   endif

   if(optSolver_=='hypre') then
      if(listaSolverDisponivel_(3).eqv..false.) then
         print*, "O Solver escolhido ...,  ", optSolver_, ",  não está disponível"
         stop 100
      endif
   endif
    
end subroutine      

!==========================================================================================================    
subroutine identificaSolversDisponiveis(listaSolverDisponivel_)
   logical, intent(inout) :: listaSolverDisponivel_(*)
         
   listaSolverDisponivel_(1)=.true. !skyline
#ifdef withPardiso
   listaSolverDisponivel_(2)=.true. !pardiso
#else
   listaSolverDisponivel_(2)=.false. !pardiso
#endif
#ifdef withHYPRE
   listaSolverDisponivel_(3)=.true. !hypre
#else
   listaSolverDisponivel_(3)=.false. !hypre
#endif

   print*, "Solvers disponiveis:"
   print*, "                      SKYLINE"
   if(listaSolverDisponivel_(2).eqv..true.) print*, "                      PARDISO"
   if(listaSolverDisponivel_(3).eqv..true.) print*, "                      HYPRE"

end subroutine

!==========================================================================================================    

subroutine lerPropriedadesCarstesElipses()
   use mMalha,            only: numel, conecNodaisElem, x
   use mPotencial,        only: estrutSistEqP

   implicit none
          
   ! real*8 :: domX, domY
   ! real*8 :: coordX, coordY
   ! integer :: i, j, cont
   ! integer :: contEsq, contDir, contSup, contInf

   integer :: no, I, nel, matMatriz=1
   integer :: inFile, numObjetos
   real*8  :: eps=1.d-10, eps2=5.d-1, x_, y_
          
   real*8,   allocatable ::  objetoX(:),  objetoY(:),  objetoW(:),  objetoH(:), objetoCor(:)
   integer,  allocatable ::  objetoNum(:), objetoMat(:)

   ! Depois, quando a rotina for para seu lugar, tirar essas linha
   !       integer, allocatable :: CCDir(:), CCEsq(:), CCSup(:), CCInf(:)
   !       integer :: numNosDir, numNosEsq, numNosSup, numNosInf

   write(*,*) "Inside lerPropriedadesCarstesElipses"
          
   inFile = 30;

   ! Borda direita
   open(inFile, file='propriedadesElipses.dat')
   read(inFile, '(I10)') numObjetos 
   write(*,*) "numObjetos", numObjetos

   allocate(objetoX(numObjetos));   objetoX=  0.0;
   allocate(objetoY(numObjetos));   objetoY=  0.0;
   allocate(objetoW(numObjetos));   objetoW=  0.0;
   allocate(objetoH(numObjetos));   objetoH=  0.0;
   allocate(objetoCor(numObjetos)); objetoCor=0.0;
   allocate(objetoNum(numObjetos)); objetoNum=0;
   allocate(objetoMat(numObjetos)); objetoMat=0;

   DO I=1,numObjetos
      read(inFile, 100)  objetoNum(I), objetoX(I), objetoY(i), objetoW(I), objetoH(I),  objetoCor(I), objetoMat(I)
      ! write(*,*) "linha", I
      write(*,100)      objetoNum(I), objetoX(I), objetoY(i), objetoW(I), objetoH(I),  objetoCor(I), objetoMat(I)
   END DO

   ! estrutSistEqP%u(1,no)=0.55           
   ! write(*,*) "elem, no1", nel, no
      
   do nel=1,numel
      matMatriz = 1
      no = conecNodaisElem(1,nel)
      estrutSistEqP%u(1,no)=1.0           
      ! write(*,*) "elem, no1", nel, no
            
      DO I=1,numObjetos
         x_= x(1,no)-objetoX(I)
         y_= x(2,no)-objetoY(I)            
         ! write(*,*) "nel, no, x, y, elipse", nel, no, x(1,no), x(2,no)
         ! write(*,*) "x_*x_", x_*x_
         ! write(*,*) "y_*y_", y_*y_
         ! write(*,*) "objetoW(I)*objetoW(I)", objetoW(I)*objetoW(I)
         ! write(*,*) "objetoH(I)*objetoH(I)", objetoH(I)*objetoH(I)
         ! write(*,*) "resultado",( ((x_*x_)/(objetoW(I)*objetoW(I))) + ((y_*y_)/(objetoH(I)*objetoH(I))) )
         if ( ( ((x_*x_)/(objetoW(I)*objetoW(I))) + ((y_*y_)/(objetoH(I)*objetoH(I))) ) .LE. 1.d0) then
            estrutSistEqP%u(1,no)=objetoCor(I)
            if ( x(2,no) .LE. objetoY(I) ) then
               estrutSistEqP%u(1,no) = 0.2
               matMatriz = 2
            else if ( x(2,no) .LE. (objetoY(I) + (objetoH(I)/2.0)) ) then
               estrutSistEqP%u(1,no) = 0.4
               matMatriz = 3
            else     
               estrutSistEqP%u(1,no) = 0.6
               matMatriz = 4                         
            endif    
            ! estrutSistEqP%u(1,no)=matMatriz !objetoCor(I)
            ! write(*,*) "nel, no, x, y, elipse", nel, no, x(1,no), x(2,no)
            ! write(*,*) "x_*x_", x_*x_
            ! write(*,*) "y_*y_", y_*y_
            ! write(*,*) "objetoW(I)*objetoW(I)", objetoW(I)*objetoW(I)
            ! write(*,*) "objetoH(I)*objetoH(I)", objetoH(I)*objetoH(I)
            ! write(*,*) "resultado",( ((x_*x_)/(objetoW(I)*objetoW(I))) + ((y_*y_)/(objetoH(I)*objetoH(I))) )
            exit;
         else 
            matMatriz = 1
            ! if ( (abs(x(1,no)-objetoX(I)) .LT. eps2) .AND. ( (x(2,no) .GT. objetoY(I)) .AND. (x(2,no) .LE. (objetoY(I)+objetoH(I)) ) )    ) then
            !    estrutSistEqP%u(1,no)=objetoCor(I)
            !    estrutSistEqP%u(1,conecNodaisElem(2,nel))=objetoCor(I)
            !    estrutSistEqP%u(1,conecNodaisElem(3,nel))=objetoCor(I)
            !    estrutSistEqP%u(1,conecNodaisElem(4,nel))=objetoCor(I)
            ! end if          
         end if                                
      END DO
      write(100,200) nel, matmatriz, conecNodaisElem(1,nel), conecNodaisElem(2,nel), &
               conecNodaisElem(3,nel), conecNodaisElem(4,nel)                        
   end do
          
   close(inFile)

   100  format(I10, f10.2, f10.2, f10.2, f10.2, f10.2, I10)
   200  format(I10, I10, I10, I10, I10, 1I10, I10)
end subroutine

!==========================================================================================================    
subroutine alteraPermeabilidadeDeElementosPorCoordenada()
   use mMalha,           only: conecNodaisElem, x, nsd, nen, numel
   use mGlobaisArranjos, only: campoPermeabilidade
    
   implicit none
      
   real*8  :: xCarste, yCarste, zCarste    
   real*8  :: permCarsteX, permCarsteY, permCarsteZ
   integer :: no1, no2, no3, no4, no5, no6, no7, no8
   integer :: i, nel

   ! Localiza coordenadas de Carste e define permeabilidade
    
   do i=1, 1
      xCarste=0.025
      yCarste=0.025
      zCarste=0.025
      permCarsteX=10.0
      permCarsteY=10.0
      permCarsteZ=10.0
         
      if(nsd==2 .and. nen==4) then
         do nel=1, numel
            no1=conecNodaisElem(1,nel)
            no2=conecNodaisElem(2,nel)
            no3=conecNodaisElem(3,nel)
            no4=conecNodaisElem(4,nel)
            if(x(1,no1)<=xCarste .and. x(1,no2)>=xCarste) then
               if(x(2,no1)<=yCarste .and. x(2,no4)>=yCarste) then
                  campoPermeabilidade(1,nel) = permCarsteX
                  campoPermeabilidade(2,nel) = permCarsteY
                  campoPermeabilidade(3,nel) = permCarsteZ
               endif
            endif
         end do
      endif

      if(nsd==2 .and. nen==3) then
         do nel=1, numel
            no1=conecNodaisElem(1,nel)
            no2=conecNodaisElem(2,nel)
            no3=conecNodaisElem(3,nel)
            if(x(1,no1)<=xCarste .and. x(1,no2)>=xCarste) then
               if(x(2,no1)<=yCarste .and. x(2,no3)>=yCarste) then
                  campoPermeabilidade(1,nel) = permCarsteX
                  campoPermeabilidade(2,nel) = permCarsteY
                  campoPermeabilidade(3,nel) = permCarsteZ
               endif
            endif
         end do
      endif

      if(nsd==3 .and. nen==8) then
         do nel=1, numel
            no1=conecNodaisElem(1,nel)
            no2=conecNodaisElem(2,nel)
            no3=conecNodaisElem(3,nel)
            no4=conecNodaisElem(4,nel)
            no5=conecNodaisElem(5,nel)
            no6=conecNodaisElem(6,nel)
            no7=conecNodaisElem(7,nel)
            no8=conecNodaisElem(8,nel)
            if(x(1,no1)<=xCarste .and. x(1,no2)>=xCarste) then
               if(x(2,no1)<=yCarste .and. x(2,no4)>=yCarste) then
                  if(x(3,no1)<=zCarste .and. x(3,no5)>=zCarste) then
                     campoPermeabilidade(1,nel) = permCarsteX
                     campoPermeabilidade(2,nel) = permCarsteY
                     campoPermeabilidade(3,nel) = permCarsteZ
                  endif
               endif
            endif
         end do
      endif   

      if(nsd==3 .and. nen==4) then
         do nel=1, numel
            no1=conecNodaisElem(1,nel)
            no2=conecNodaisElem(2,nel)
            no3=conecNodaisElem(3,nel)
            no4=conecNodaisElem(4,nel)
            
            ! elemento1
            if(x(1,no1)<=xCarste .and. x(1,no2)>=xCarste) then
               if(x(2,no2)<=yCarste .and. x(2,no3)>=yCarste) then
                  if(x(3,no1)<=zCarste .and. x(3,no4)>=zCarste) then
                     campoPermeabilidade(1,nel) = permCarsteX
                     campoPermeabilidade(2,nel) = permCarsteY
                     campoPermeabilidade(3,nel) = permCarsteZ   
                  endif
               endif
            endif
               
            ! elemento 2
            if(x(1,no3)<=xCarste .and. x(1,no1)>=xCarste) then
               if(x(2,no1)<=yCarste .and. x(2,no2)>=yCarste) then
                  if(x(3,no1)<=zCarste .and. x(3,no4)>=zCarste) then
                     campoPermeabilidade(1,nel) = permCarsteX
                     campoPermeabilidade(2,nel) = permCarsteY
                     campoPermeabilidade(3,nel) = permCarsteZ
                  endif
               endif
            endif
            
            ! elemento 3
            if(x(1,no3)<=xCarste .and. x(1,no4)>=xCarste) then
               if(x(2,no4)<=yCarste .and. x(2,no2)>=yCarste) then
                  if(x(3,no1)<=zCarste .and. x(3,no2)>=zCarste) then
                     campoPermeabilidade(1,nel) = permCarsteX
                     campoPermeabilidade(2,nel) = permCarsteY
                     campoPermeabilidade(3,nel) = permCarsteZ
                  endif
               endif
            endif
               
            ! elemento 4
            if(x(1,no2)<=xCarste .and. x(1,no4)>=xCarste) then
               if(x(2,no3)<=yCarste .and. x(2,no2)>=yCarste) then
                  if(x(3,no1)<=zCarste .and. x(3,no4)>=zCarste) then
                     campoPermeabilidade(1,nel) = permCarsteX
                     campoPermeabilidade(2,nel) = permCarsteY
                     campoPermeabilidade(3,nel) = permCarsteZ
                  endif
               endif
            endif
            
            ! elemento 5
            if(x(1,no2)<=xCarste .and. x(1,no1)>=xCarste) then
               if(x(2,no3)<=yCarste .and. x(2,no4)>=yCarste) then
                  if(x(3,no2)<=zCarste .and. x(3,no4)>=zCarste) then
                     campoPermeabilidade(1,nel) = permCarsteX
                     campoPermeabilidade(2,nel) = permCarsteY
                     campoPermeabilidade(3,nel) = permCarsteZ
                  endif
               endif
            endif     
               
            ! elemento 6
            if(x(1,no3)<=xCarste .and. x(1,no2)>=xCarste) then
               if(x(2,no1)<=yCarste .and. x(2,no3)>=yCarste) then
                  if(x(3,no1)<=zCarste .and. x(3,no4)>=zCarste) then
                     campoPermeabilidade(1,nel) = permCarsteX
                     campoPermeabilidade(2,nel) = permCarsteY
                     campoPermeabilidade(3,nel) = permCarsteZ
                  endif
               endif
            endif   
         end do
      endif 
   end do
end subroutine
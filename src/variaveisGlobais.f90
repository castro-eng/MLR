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
module mGlobaisArranjos
   
   integer*4, allocatable :: npar(:)
   real*8                 :: etime(6)
   character(len=20)   :: title
   integer*4, allocatable :: mat(:), matFratura(:), tiposElemBordas(:)
   real*8,    allocatable :: grav(:), bf(:,:), c(:,:)
   integer,   allocatable :: CCDir(:), CCEsq(:), CCSup(:), CCInf(:)
   real*8,    allocatable :: CCDir_fratura(:), CCEsq_fratura(:), CCSup_fratura(:), CCInf_fratura(:)
   real*8,    allocatable :: aberturasBorda(:,:), alturasBorda(:,:)
   logical :: listaSolverDisponivel(3)
   real*8,    allocatable :: campoPermeabilidade(:,:)
   ! real*8 :: vazao(6) ! Versão inicial
   real*8 :: vazao(14)     ! 6 bordas no 3d mais uma borda de poço no caso de WI 

   ! Contador para numero de elementos de descontinuidades nos nós do bordo. 
   integer*4, allocatable :: NelmDescNoBordo(:,:)

   ! Fluxo consistente
   integer*4, allocatable :: matPoco(:), matDescontNaBorda(:)
   integer*4, allocatable :: tipoBordaElemDescontNaBorda(:), matBordasDirichlet(:)
   logical :: bordaPrescrita(14)

end module ! mGlobaisArranjos


module mGlobaisEscalares
   
   integer*4:: exec,iprtin
   integer*4:: numParElem=15
   integer*4:: ntype,numat,nrowsh,nicode
   integer*4:: npint, npintFratura, npintFraturasNaBorda, nrowshFratNaBorda
   integer, parameter :: numProps = 6

   real*8, parameter  :: zero=0.0d0, pt25=0.25d0, pt5=0.5d0
   real*8, parameter  :: four=4.0d0, five=5.0d0, six=6.0d0
   real*8, parameter  :: pt1667=0.1666666666666667d0
   real*8, parameter  :: pt8=0.8d0
   real*8, parameter  :: pt45=0.45d0
   real*8, parameter  :: one=1.0d0, two=2.0d0, three=3.0d0
   real*8, parameter  :: pi=3.14159265359
   real*8, parameter  :: corrTempo = 3600. ! Variável mudar a escala de tempo de simulação, no SI (Ex. segundo->hora)
   integer*4:: id0
   integer*4 :: mpi_comm
        
   real*8   :: relax_param, tolerancia_relativa
   INTEGER  :: numPassosTempo, passoTempo, max_iter_pressao
   real*8   :: tempo, tempoTotalSimulacao, deltaT
        
   real*8   :: permX, permY, carregamento, carregamento_ref
   logical  :: transiente, simulandoReferencia
        
   real*8   :: areaCelula, areaFluido
   integer :: numNosDir, numNosEsq, numNosSup, numNosInf
   
   logical :: saltoPressao
   logical :: fazerWellIndex

   real*8, parameter :: perm_x_efetiva_ref = 1.3542196345814103d-15
   ! Permeabilidade de referencia
   real*8 :: perm_y_efetiva_ref
   ! Porosidade de referencia
   real*8 :: porosidadeCelula_ref

   ! real*8, parameter :: viscosidade = 1.0e-03 ! Visc. Agua
   real*8, parameter :: viscosidade = 2.0e-02
   ! real*8, parameter :: viscosidade = 1.0
   
   real*8 :: BIOTCOEF, BIOTMOD, BULKROCK, BULKGRAIN, BULKOIL

   logical :: flagTransmib ! Flag para calculo de transmissibilidade

end module mGlobaisEscalares



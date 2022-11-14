!=================================================================================
!
!         Rotinas para escrita de arquivo hdf5 para o programa de
!         elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!         usa implementações da HDF5 Library version 1.8
!
!         novo projeto: 
!         Eduardo Castro,        ecastro@lncc.br [4]
!
!         LNCC/MCT
!         Petropolis, 07.2019
!=================================================================================
module mVisualizacaoHdf5
   
   use hdf5
   implicit none

   ! Declaracao de variaveis visiveis
   character(len=16), parameter :: filename = "Experimento"
   character(len=16), allocatable :: attr_names(:), name_struct(:)

   integer(HID_T), dimension(8) :: p_root_hdf5
   integer(HID_T), allocatable :: p_tree_hdf5(:)

   integer :: flag_err ! identificador de erro
   integer :: num_result, num_tipo_elem, npoint_estr
   ! integer, dimension(3) :: attr_types = (/1,7,-6/)

   ! Funcoes e Subrotinas
   private :: hdf5_inic, hdf5_close
   public :: hdf5_basic_structure, hdf5_close_structure
   public :: hdf5_dataset_integer, hdf5_dataset_real, hdf5_dataset_double
   public :: xdmf_basic_structure
   public :: hdf5_dataset_dvec, hdf5_structure_addPressure, hdf5_structure_addVelocity
   public :: hdf5_structure_addMeshFalha, hdf5_structure_addPropFalha
   public :: hdf5_structure_addVelocBlocoElem, hdf5_structure_addVelocBlocoFalha
   public :: xdmf_structure_bloco

   public :: hdf5_structure_addPressureFalha, hdf5_structure_addVelocFalha

   ! Escrita do arquivo xdmf
   ! CALL xdmf_basic_structure(attr_names, attr_types, 3, filename)

   ! Declaracao das rotinas do modulo
   contains
   !=======================================================================
   !-------------- ROTINAS DE CONTROLE DO MODULO --------------------------
   !=======================================================================
   subroutine hdf5_inic()

      integer i,j,k

      num_result=4
      num_tipo_elem=2
      npoint_estr=num_result*(num_tipo_elem+1)

      allocate(attr_names(num_result))
      allocate(name_struct(num_tipo_elem))
      allocate(p_tree_hdf5(npoint_estr))

      name_struct(1)="struct1"
      name_struct(2)="struct2"

      attr_names(1)="material"
      attr_names(2)="permeability"
      attr_names(3)="pressure"
      attr_names(4)="velocity"

      call h5gcreate_f(p_root_hdf5(8),attr_names(1), p_tree_hdf5(1), flag_err)
      call h5gcreate_f(p_root_hdf5(8),attr_names(2), p_tree_hdf5(2), flag_err)
      call h5gcreate_f(p_root_hdf5(7),attr_names(3), p_tree_hdf5(3), flag_err)
      call h5gcreate_f(p_root_hdf5(7),attr_names(4), p_tree_hdf5(4), flag_err)
      do i=1,2
         do j=1,4
            k= num_result*(i-1)+j+4
            call h5gcreate_f(p_tree_hdf5(j),name_struct(i), p_tree_hdf5(k), flag_err)
         end do
      end do
      return

   end subroutine
   !=======================================================================
   subroutine hdf5_close()

      integer i
      do i=npoint_estr,1,-1
         call h5gclose_f(p_tree_hdf5(i), flag_err)
      end do

      deallocate(attr_names)
      deallocate(name_struct)
      deallocate(p_tree_hdf5)
      return

   end subroutine

   !=======================================================================
   !---------------- ROTINAS BASICAS PARA O HDF5 --------------------------
   !=======================================================================
   subroutine hdf5_basic_structure()
      ! Rotina para criação das estruturas de pasta básica para o arquivo hdf5 
      ! Cria uma estrutura de 3 pastas e retornas os ponteiros de cada uma.
      ! As pastas criadas tem as seguintes funções:
      ! -    TIME - armazena em um vetor a lista dos valores tempos utilizados
      !             na análise.
      ! -    MESH - armazena informação das coordenadas nodais e das
      !             conectividades dos elementos.
      ! - RESULTS - armazena informações relativas aos resultados da análise.
      ! Podendo estes valores estarem situados nos nós ou nos 
      ! elementos.

      implicit none
      ! integer(hid_t), dimension(*), intent(out) :: p_dados
      ! integer, intent(in):: n_pointer
      ! integer :: flag_err ! identificador de erro

      ! Initializando a interface para o FORTRAN.
      call h5open_f(flag_err)
      ! Criação de um novo arquivo utilizando propiedades padrões.
      call h5fcreate_f('./out/'//trim(filename)//'.h5', h5f_acc_trunc_f, p_root_hdf5(1), flag_err)
      ! Criação das pastas de grupo.
      call h5gcreate_f(p_root_hdf5(1),     "Mesh", p_root_hdf5(2), flag_err)
      call h5gcreate_f(p_root_hdf5(1),  "Results", p_root_hdf5(3), flag_err)
      call h5gcreate_f(p_root_hdf5(1),     "Time", p_root_hdf5(4), flag_err)
      call h5gcreate_f(p_root_hdf5(2), "Geometry", p_root_hdf5(5), flag_err)
      call h5gcreate_f(p_root_hdf5(2), "Topology", p_root_hdf5(6), flag_err)
      call h5gcreate_f(p_root_hdf5(3),     "Node", p_root_hdf5(7), flag_err)
      call h5gcreate_f(p_root_hdf5(3),     "Cell", p_root_hdf5(8), flag_err)

      ! Criação de pastas para resultados.
      call hdf5_inic()
      return
   end subroutine

   !=======================================================================
   subroutine hdf5_close_structure()
      ! Rotina para fechar os espaços indicados pelos ponteiros no arquivo hdf5
      implicit none
      integer :: i

      call hdf5_close()
      do i=8,2,-1
         call h5gclose_f(p_root_hdf5(i), flag_err)
      end do
      call h5fclose_f(p_root_hdf5(1), flag_err)
      return
   end subroutine

   !=======================================================================
   !---------------- ROTINAS PARA CRIAR DATASET ---------------------------
   !=======================================================================
   subroutine hdf5_dataset_integer(vector_name, vector_value, m, n, hdf5_pointer)
      ! Rotina para criação de dataset. 
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"

      character(len=*), intent(in) :: vector_name
      integer, dimension(*), intent(in) :: vector_value
      integer(HID_T), intent(in) :: hdf5_pointer
      integer, intent(in):: m, n

      integer(HSIZE_T), dimension(2) :: dims
      integer(HID_T) :: p_space, p_valores
      integer :: error

      dims = (/m,n/)

      ! Cria o espaço para o dataspace.
      call h5screate_simple_f(2, dims, p_space, error)

      ! Define as propriedades e posição do dataspace.
      call h5dcreate_f(hdf5_pointer, trim(vector_name), H5T_NATIVE_INTEGER, &
               p_space, p_valores, error)

      ! Copia o valores do vetor/matriz para o dataspace.
      call h5dwrite_f(p_valores, H5T_NATIVE_INTEGER, vector_value, dims, error)

      ! Fecha o acesso ao dataspace.
      call h5dclose_f(p_valores, error)
      call h5sclose_f(p_space, error)
      return
   end subroutine

   !=======================================================================
   SUBROUTINE hdf5_dataset_real(vector_name, vector_value, m, n, hdf5_pointer)
      ! Rotina para criação de dataset. 
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de real com o tamanho do vetor/matriz "vector_value"

      character(len=*), intent(in) :: vector_name
      real, dimension(*), intent(in) :: vector_value
      integer(HID_T), intent(in) :: hdf5_pointer
      integer, intent(in):: m, n

      integer(HSIZE_T), dimension(2) :: dims
      integer(HID_T) :: p_space, p_valores
      integer :: error

      dims = (/m,n/)

      ! Cria o espaço para o dataspace.
      call h5screate_simple_f(2, dims, p_space, error)

      ! Define as propriedades e posição do dataspace.
      call h5dcreate_f(hdf5_pointer, trim(vector_name), H5T_NATIVE_REAL, &
               p_space, p_valores, error)

      ! Copia o valores do vetor/matriz para o dataspace.
      call h5dwrite_f(p_valores, H5T_NATIVE_REAL, vector_value, dims, error)

      ! Fecha o acesso ao dataspace.
      call h5dclose_f(p_valores, error)
      call h5sclose_f(p_space, error)
      return
   end subroutine

   !=======================================================================
   SUBROUTINE hdf5_dataset_double(vector_name, vector_value, m, n, hdf5_pointer)
      ! Rotina para criação de dataset. 
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de double com o tamanho do vetor/matriz "vector_value"

      character(len=*), intent(in) :: vector_name
      double precision, dimension(*), intent(in) :: vector_value
      integer(HID_T), intent(in) :: hdf5_pointer
      integer, intent(in):: m, n

      integer(HSIZE_T), dimension(2) :: dims
      integer(HID_T) :: p_space, p_valores
      integer :: error

      dims = (/m,n/)

      ! cria o espaço para o dataspace.
      call h5screate_simple_f(2, dims, p_space, error)

      ! define as propriedades e posição do dataspace.
      call h5dcreate_f(hdf5_pointer, trim(vector_name), H5T_NATIVE_DOUBLE, &
               p_space, p_valores, error)

      ! copia o valores do vetor/matriz para o dataspace.
      call h5dwrite_f(p_valores, H5T_NATIVE_DOUBLE, vector_value, dims, error)

      ! fecha o acesso ao dataspace.
      call h5dclose_f(p_valores, error)
      call h5sclose_f(p_space, error)
      return
   end subroutine

   !=======================================================================
   !------------------ ROTINAS ESPECIALIZADAS -----------------------------
   !=======================================================================
   SUBROUTINE hdf5_dataset_dvec(vector_name, vector_value, m, n, hdf5_pointer)
      ! Rotina para criação de dataset. 
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de double com o tamanho do vetor/matriz "vector_value"

      character(len=*), intent(in) :: vector_name
      double precision, dimension(m, n), intent(in) :: vector_value
      integer(HID_T), intent(in) :: hdf5_pointer
      integer, intent(in):: m, n

      integer :: i,j
      real*8, allocatable :: array_value(:,:)

      if (m.gt.3) return
      if (m.lt.1) return
      if (m.eq.3) then
         call hdf5_dataset_double(vector_name, vector_value, m, n, hdf5_pointer)
         return
      end if

      allocate(array_value(3, n))
      do j=1,n
         do i=1,m
            array_value(i,j)=vector_value(i,j)
         end do
         do i=m+1,3
            array_value(i,j)=0.D0
         end do
      end do
      call hdf5_dataset_double(vector_name, array_value, 3, n, hdf5_pointer)
      deallocate(array_value)

      return
   end subroutine

   !=======================================================================
   subroutine hdf5_structure_addPressure(vector_name, hdf5_pointer)
      ! Rotina para inclusão do contorno na solução do problema. 
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mPotencial,        only: estrutSistEqP
      use mMalha,            only: numnp!, numel

      character(len=*), intent(in) :: vector_name
      integer(HID_T), intent(in) :: hdf5_pointer

      character(len=16):: name_data
      integer :: i,j
      real*8, allocatable :: array_value(:)

      allocate(array_value(numnp))
      do i=1,estrutSistEqP%ndof
         write (name_data,"(I5,2A)") i, "T", vector_name
         do j=1,numnp
            array_value(j)=estrutSistEqP%u(i,j)
         end do
         call hdf5_dataset_double(adjustl(name_data), array_value, 1, numnp, hdf5_pointer)
      end do
      deallocate(array_value)
      return
   end subroutine

   !=======================================================================
   subroutine hdf5_structure_addVelocity(vector_name, hdf5_pointer)
      ! Rotina para inclusão do contorno na solução do problema. 
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mFluxo,            only: estrutSistEqF
      use mMalha,            only: numnp, nsd!, numel

      character(len=*), intent(in) :: vector_name
      integer(HID_T), intent(in) :: hdf5_pointer

      integer :: i,j
      real*8, allocatable :: array_value(:,:)

      allocate(array_value(3,numnp))
      do i=1,nsd
         do j=1,numnp
            array_value(i,j)=estrutSistEqF%u(i,j)
         end do
      end do
      if (nsd/=3) then
         do i=nsd+1,3
            do j=1,numnp
               array_value(i,j)=0.
            end do
         end do
      end if
      call hdf5_dataset_double(adjustl(vector_name), array_value, 3, numnp, hdf5_pointer)
      deallocate(array_value)

      return
   end subroutine

   !=======================================================================
   subroutine hdf5_structure_addMeshFalha(vector_name)
      ! Rotina para inclusão das coordenadas e das conectividades da descontinuidade
      ! Cria os dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mMalha,    only: x, nsd, conecNodaisFratura, numelFratura, nenFratura!, numnpFratura, numnp

      character(len=*), intent(in) :: vector_name

      logical :: flag
      integer :: i,j,k,pos_conc,pos_array, tam_array
      integer, allocatable :: array_conect(:,:), array_mapeamento(:)
      real*8, allocatable :: array_coord(:,:)

      ! Montagem do mapeamento
      tam_array=0 
      allocate(array_mapeamento(numelFratura*nenFratura))
      do i=1,numelFratura
         do j=1,nenFratura
            pos_conc = conecNodaisFratura(j,i)
            flag=.false.
            k=1
            do 
               if (flag.eqv..true.) exit
               if (pos_conc.eq.array_mapeamento(k)) flag=.true.
               if (k.eq.(numelFratura*nenFratura)) flag=.true.
               if (k.gt.tam_array) then
                  array_mapeamento(k)=pos_conc
                  tam_array=k
                  flag=.true.
               end if
               k=k+1
            end do
         end do
      end do
      ! numnpFratura = tam_array

      ! Montagem e escrita das coordenadas
      allocate(array_coord(nsd,tam_array))
      do i=1,tam_array
         do j=1,nsd
            array_coord(j,i)=x(j,array_mapeamento(i))
         end do
      end do
      call hdf5_dataset_double(adjustl(vector_name), array_coord, nsd, tam_array, p_root_hdf5(5))
      deallocate(array_coord)

      ! Montagem e escrita das conectividades
      pos_array=0 
      allocate(array_conect(nenFratura,numelFratura))
      do i=1,numelFratura
         do j=1,nenFratura
            pos_conc = conecNodaisFratura(j,i)
            flag=.false.
            k=1
            do 
               if (flag.eqv..true.) exit
               if (k.gt.tam_array) flag=.true.
               if (pos_conc.eq.array_mapeamento(k)) then
                  flag=.true.
                  array_conect(j,i)=k-1
               end if
               k=k+1
            end do
         end do
      end do
      call hdf5_dataset_integer(adjustl(vector_name), array_conect, nenFratura, numelFratura, p_root_hdf5(6))
      deallocate(array_conect)
      deallocate(array_mapeamento)

      return
   end subroutine

   !=======================================================================
   subroutine hdf5_structure_addPropFalha(vector_name)
      ! Rotina para inclusão das propriedades (material e permeabilidade) para descontinuidade
      ! Cria os dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mMalha,            only: nsd, numelFratura
      use mGlobaisArranjos,  only: c, matFratura

      character(len=*), intent(in) :: vector_name

      integer :: i,j,m
      integer, allocatable :: array_material(:)
      real*8, allocatable :: array_permeab(:,:)

      allocate(array_material(numelFratura))
      allocate(array_permeab(nsd,numelFratura))
      do i=1, numelFratura
         m=matFratura(i)
         array_material(i)=m
         do j=1,nsd
            array_permeab(j,i) = c(j+2,m)
         end do
      end do
      call hdf5_dataset_integer( "0", array_material, 1, numelFratura, p_tree_hdf5(9))
      call hdf5_dataset_dvec(  "0", array_permeab, nsd, numelFratura, p_tree_hdf5(10))
      deallocate(array_material)
      deallocate(array_permeab)

   end subroutine

   !=======================================================================
   subroutine hdf5_structure_addPressureFalha(vector_name, hdf5_pointer)
      ! Rotina para inclusão da pressão na descontinuidade
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mPotencial,        only: estrutSistEqP
      ! use mMalha,            only: numnp, numel, nsd
      use mMalha,            only: conecNodaisFratura, numelFratura, numnpFratura, nenFratura

      character(len=*), intent(in) :: vector_name
      integer(HID_T), intent(in) :: hdf5_pointer

      character(len=16):: name_data
      logical :: flag
      integer :: i,j,k,pos_conc, tam_array!, pos_array
      integer, allocatable :: array_mapeamento(:)
      real*8, allocatable :: array_value(:)

      ! Montagem do mapeamento
      tam_array=0 
      allocate(array_mapeamento(numelFratura*nenFratura))
      do i=1,numelFratura
         do j=1,nenFratura
            pos_conc = conecNodaisFratura(j,i)
            flag=.false.
            k=1
            do 
               if (flag.eqv..true.) exit
               if (pos_conc.eq.array_mapeamento(k)) flag=.true.
               if (k.eq.(numelFratura*nenFratura)) flag=.true.
               if (k.gt.tam_array) then
                  array_mapeamento(k)=pos_conc
                  tam_array=k
                  flag=.true.
               end if
               k=k+1
            end do
         end do
      end do

      allocate(array_value(numnpFratura))
      do i=1,estrutSistEqP%ndof
         write (name_data,"(I5,2A)") i, "T", vector_name
         do j=1,numnpFratura
            array_value(j)=estrutSistEqP%u(i,array_mapeamento(j))
         end do
         call hdf5_dataset_double(adjustl(name_data), array_value, 1, numnpFratura, hdf5_pointer)
      end do
      deallocate(array_mapeamento)
      deallocate(array_value)

      return
   end subroutine

      !=======================================================================
   subroutine hdf5_structure_addVelocFalha(vector_name, hdf5_pointer)
      ! Rotina para inclusão da velocidade na descontinuidade
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mPotencial,        only: estrutSistEqP
      use mMalha,            only: nsd
      use mMalha,            only: conecNodaisFratura, numelFratura, numnpFratura, nenFratura
      use mFluxo,            only: vm_reduzido

      character(len=*), intent(in) :: vector_name
      integer(HID_T), intent(in) :: hdf5_pointer

      character(len=16):: name_data
      integer :: i,j
      real*8, allocatable :: array_value(:,:)

      allocate(array_value(3,numnpFratura))
      do i=1,estrutSistEqP%ndof
         write (name_data,"(I5,2A)") i, "T", vector_name
         do j=1,numnpFratura
            array_value(1,j)=vm_reduzido(1,j,i)
            array_value(2,j)=vm_reduzido(2,j,i)
            if(nsd==2) then
               array_value(3,j)=0d0
            else
               array_value(3,j)=vm_reduzido(3,j,i)
            endif
         end do
         call hdf5_dataset_double(adjustl(name_data), array_value, 3, numnpFratura, hdf5_pointer)
      end do
      deallocate(array_value)

   end subroutine

   !=======================================================================
   subroutine hdf5_structure_addVelocBlocoElem(vector_name, velc, hdf5_pointer)
      ! Rotina para inclusão da velocidade média da célula
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mMalha,            only: numel, nen, nsd
      implicit none

      character(len=*), intent(in) :: vector_name
      integer(HID_T), intent(in) :: hdf5_pointer
      double precision, dimension(nen,numel,nsd,1), intent(in) :: velc

      integer :: i, j, k
      double precision, allocatable :: array_value(:,:)

      allocate(array_value(3,numel))
      do i=1,3
         do j=1,numel
            array_value(i,j)=0.D0
         end do
      end do
      do i=1,nsd
         do j=1,numel
            do k=1,nen
               array_value(i,j)=array_value(i,j)+velc(k,j,i,1)
            end do
            array_value(i,j)=array_value(i,j)/dble(nen)
         end do
      end do
      call hdf5_dataset_double(adjustl(vector_name), array_value, 3, numel, hdf5_pointer)
      deallocate(array_value)
      return
   end subroutine

   !=======================================================================
   subroutine hdf5_structure_addVelocBlocoFalha(vector_name, velc, nlayer, hdf5_pointer)
      ! Rotina para inclusão da velocidade média da célula
      ! Cria um dataset com as seguintes propriedades:
      ! - nome do dataset dado na variável "vector_name".
      ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
      ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"
      use mMalha,            only: nsd!,numel, nen
      use mMalha,            only: numnpFratura!, numelFratura, nenFratura
      implicit none

      character(len=*), intent(in) :: vector_name
      integer(HID_T), intent(in) :: hdf5_pointer
      integer, intent(in) :: nlayer
      double precision, dimension(nsd,numnpFratura,nlayer), intent(in) :: velc

      character(len=16):: name_data!, temp1, temp2
      integer :: i, j, lay!, k
      double precision, allocatable :: array_value(:,:)

      allocate(array_value(3,numnpFratura))
      do lay=1,nlayer
         do j=1,numnpFratura
            do i=1,nsd
               array_value(i,j)=velc(i,j,lay)
            end do
            if (nsd.lt.3) then
               do i=nsd+1,3
                  array_value(i,j)=0.D0
               end do
            end if
         end do
         write (name_data,"(A, I2.2)") vector_name, lay
         call hdf5_dataset_double(adjustl(name_data), array_value, 3, numnpFratura, hdf5_pointer)
      end do
      deallocate(array_value)
      return

   end subroutine

   !=======================================================================
   ! subroutine hdf5_structure_addVelocBlocoFalha(vector_name, velc, nlayer, hdf5_pointer)
   !    ! Rotina para inclusão da velocidade média da célula
   !    ! Cria um dataset com as seguintes propriedades:
   !    ! - nome do dataset dado na variável "vector_name".
   !    ! - o dataset é criado no local indicado pelo pointeiro "hdf5_pointer"
   !    ! - o dataset é de inteiros com o tamanho do vetor/matriz "vector_value"

   !    use mMalha,            only: numel, nen, nsd
   !    use mMalha,            only: numelFratura, nenFratura, numnpFratura
   !    implicit none

   !    character(len=*), intent(in) :: vector_name
   !    integer(HID_T), intent(in) :: hdf5_pointer
   !    integer, intent(in) :: nlayer
   !    double precision, dimension(nenFratura,numelFratura,nsd,nlayer), intent(in) :: velc

   !    character(len=16):: name_data, temp1, temp2
   !    integer :: i, j, k, lay
   !    double precision, allocatable :: array_value(:,:)

   !    allocate(array_value(3,numelFratura))
   !    do i=1,3
   !       do j=1,numelFratura
   !          array_value(i,j)=0.D0
   !       end do
   !    end do
   !    do lay=1,nlayer
   !       do i=1,nsd
   !          do j=1,numelFratura
   !             do k=1,nenFratura
   !                array_value(i,j)=array_value(i,j)+velc(k,j,i,lay)
   !             end do
   !             array_value(i,j)=array_value(i,j)/dble(nen)
   !          end do
   !       end do
   !       write (name_data,"(A, I2.2)") vector_name, lay
   !       call hdf5_dataset_double(adjustl(name_data), array_value, 3, numelFratura, hdf5_pointer)
   !    end do
   !    deallocate(array_value)
   !    return
   ! end subroutine

   !=======================================================================
   !---------------- ROTINAS PARA O XDMF ----------------------------------
   !=======================================================================
   SUBROUTINE xdmf_structure_bloco(nelm, nnos, nen, nsd, npassos)
      ! Rotina para criação do arquivo XDML. 
      ! Cria um arquivo xdml com as seguintes propriedades:
      ! - a malha do arquivo não se altera para cada incremento de tempo.
      ! - a malha é composta por apenas 1 tipo de elemento.
      ! - os resultados são dados para as células ou para os nós.
      use mPotencial,        only: estrutSistEqP

      INTEGER, INTENT(IN) :: nelm, nnos, nen, nsd, npassos

      CHARACTER(LEN=8) :: passo_tempo!, stringTemp, posit_data, data_valueA, data_valueB
      ! CHARACTER(LEN=25) :: type_data
      INTEGER :: idFile, i, j

      ! Abrindo arquivo
      idFile=2000
      OPEN(UNIT=idFile,FILE='./out/'//TRIM(filename)//'Bloco.xdmf',STATUS='replace')

      ! Escrita do cabeçalho do arquivo.
      WRITE(idFile,"(     A)") '<?xml version="1.0"?>'
      WRITE(idFile,"(     A)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      WRITE(idFile,"(     A)") '<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
      WRITE(idFile,"( 2X, A)") '<Domain>'
      WRITE(idFile,"( 4X, A)") '<Grid Name="ReservatorioBloco" GridType="Collection" CollectionType="Temporal">'

      ! Listagem dos tempos
      WRITE(idFile,"( 6X, A)") '<Time TimeType="List">'
      WRITE(idFile,"( 8X,A,I5,A)") '<DataItem Format="HDF" Dimensions="',npassos,'">'
      WRITE(idFile,"(10X,2A)") TRIM(filename), '.h5:/Time/TimeSerie'
      WRITE(idFile,"( 8X, A)") '</DataItem>'
      WRITE(idFile,"( 6X, A)") '</Time>'
      DO i=1,npassos
         ! Inicialização da estrutura (coord, conec, valores)
         WRITE(passo_tempo, "(I0.8)") i-1
         WRITE(idFile,"( 6X, 3A)") '<Grid Name="BlocoT',adjustl(passo_tempo),'" GridType="Uniform">'

         ! Listagem da topologia
         WRITE(idFile,"( 8X,A,I5,A,I5,A)") '<Topology NumberOfElements="',nelm, &
                  '" TopologyType="Triangle" NodesPerElement="',nen,'">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nelm,nen,'" NumberType="UInt" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Mesh/Topology/struct1'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Topology>'

         ! Listagem da geometria
         WRITE(idFile,"( 8X,A)") '<Geometry GeometryType="XY">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nnos,nsd,'" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Mesh/Geometry/struct1'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Geometry>'

         ! Listagem de material
         WRITE(idFile,"( 8X,A)") '<Attribute Name="material" AttributeType="Scalar" Center="Cell">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nelm,1,'" NumberType="UInt" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Results/Cell/material/struct1/0'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Attribute>'

         ! Listagem de permeabilidade
         WRITE(idFile,"( 8X,A)") '<Attribute Name="permeabilidade" AttributeType="Vector" Center="Cell">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nelm,3,'" NumberType="Float" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Results/Cell/permeability/struct1/0'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Attribute>'

         ! Listagem de pressão
         DO j=1,estrutSistEqP%ndof
            WRITE(idFile,"( 8X,A,I1.1,A)") '<Attribute Name="pressure_dof', j, '" AttributeType="Scalar" Center="Node">'
            WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nnos,1,'" NumberType="Float" Format="HDF">'
            WRITE(idFile,"(12X,2A,I1.1,2A)") TRIM(filename), '.h5:/Results/Node/pressure/struct1/',j,"T", adjustl(passo_tempo)
            WRITE(idFile,"(10X, A)") '</DataItem>'
            WRITE(idFile,"( 8X, A)") '</Attribute>'
         END DO

         ! Listagem de velocidade
         WRITE(idFile,"( 8X,A)") '<Attribute Name="velocity" AttributeType="Vector" Center="Node">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nnos,3,'" NumberType="Float" Format="HDF">'
         WRITE(idFile,"(12X,3A)") TRIM(filename), '.h5:/Results/Node/velocity/struct1/', adjustl(passo_tempo)
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Attribute>'

         WRITE(idFile,"( 6X,A)") '</Grid>'
      END DO
      ! fechamento do arquivo
      WRITE(idFile,"( 4X,A)") '</Grid>'
      WRITE(idFile,"( 2X,A)") '</Domain>'
      WRITE(idFile,"(    A)") '</Xdmf>'

      CLOSE(idFile)
      RETURN  
   END SUBROUTINE

   !=======================================================================
   SUBROUTINE xdmf_structure_falha(nelm, nnos, nen, nsd, npassos)
      ! Rotina para criação do arquivo XDML. 
      ! Cria um arquivo xdml com as seguintes propriedades:
      ! - a malha do arquivo não se altera para cada incremento de tempo.
      ! - a malha é composta por apenas 1 tipo de elemento.
      ! - os resultados são dados para as células ou para os nós.
      use mPotencial,        only: estrutSistEqP

      INTEGER, INTENT(IN) :: nelm, nnos, nen, nsd, npassos

      CHARACTER(LEN=6) :: passo_tempo!, stringTemp, posit_data, data_valueA, data_valueB
      ! CHARACTER(LEN=25) :: type_data
      INTEGER :: idFile, i, j

      ! Abrindo arquivo
      idFile=2000
      OPEN(UNIT=idFile,FILE='./out/'//TRIM(filename)//'Falha.xdmf',STATUS='replace')

      ! Escrita do cabeçalho do arquivo.
      WRITE(idFile,"(     A)") '<?xml version="1.0"?>'
      WRITE(idFile,"(     A)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      WRITE(idFile,"(     A)") '<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
      WRITE(idFile,"( 2X, A)") '<Domain>'
      WRITE(idFile,"( 4X, A)") '<Grid Name="ReservatorioFalha" GridType="Collection" CollectionType="Temporal">'

      ! Listagem dos tempos
      WRITE(idFile,"( 6X, A)") '<Time TimeType="List">'
      WRITE(idFile,"( 8X,A,I5,A)") '<DataItem Format="HDF" Dimensions="',npassos,'">'
      WRITE(idFile,"(10X,2A)") TRIM(filename), '.h5:/Time/TimeSerie'
      WRITE(idFile,"( 8X, A)") '</DataItem>'
      WRITE(idFile,"( 6X, A)") '</Time>'
      DO i=1,npassos
         ! Inicialização da estrutura (coord, conec, valores)
         WRITE(passo_tempo, "(I6)") i-1
         WRITE(idFile,"( 6X, 3A)") '<Grid Name="FalhaT',adjustl(passo_tempo),'" GridType="Uniform">'

         ! Listagem da topologia
         WRITE(idFile,"( 8X,A,I5,A,I5,A)") '<Topology NumberOfElements="',nelm,&
                  '" TopologyType="Polyline" NodesPerElement="',nen,'">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nelm,nen,'" NumberType="UInt" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Mesh/Topology/struct2'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Topology>'

         ! Listagem da geometria
         WRITE(idFile,"( 8X,A)") '<Geometry GeometryType="XY">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nnos,nsd,'" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Mesh/Geometry/struct2'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Geometry>'

         ! Listagem de material
         WRITE(idFile,"( 8X,A)") '<Attribute Name="material" AttributeType="Scalar" Center="Cell">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nelm,1,'" NumberType="UInt" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Results/Cell/material/struct2/0'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Attribute>'

         ! Listagem de permeabilidade
         WRITE(idFile,"( 8X,A)") '<Attribute Name="permeabilidade" AttributeType="Vector" Center="Cell">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nelm,3,'" NumberType="Float" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Results/Cell/permeability/struct2/0'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Attribute>'

         ! Listagem de pressão
         DO j=1,estrutSistEqP%ndof
            WRITE(idFile,"( 8X,A,I1.1,A)") '<Attribute Name="pres_dof_', j, '" AttributeType="Scalar" Center="Node">'
            WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nnos,1,'" NumberType="Float" Format="HDF">'
            WRITE(idFile,"(12X,2A,I1.1,2A)") TRIM(filename), '.h5:/Results/Node/pressure/struct2/',j,"T", adjustl(passo_tempo)
            WRITE(idFile,"(10X, A)") '</DataItem>'
            WRITE(idFile,"( 8X, A)") '</Attribute>'
         END DO

         ! Listagem de velocidade
         DO j=1,(estrutSistEqP%ndof-2)
            WRITE(idFile,"( 8X,A,I2.2,A)") '<Attribute Name="vel_lay_',j,'" AttributeType="Vector" Center="Node">'
            WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nnos,3,'" NumberType="Float" Format="HDF">'
            WRITE(idFile,"(12X,2A,I2.2)") TRIM(filename), '.h5:/Results/Node/velocity/struct2/VM', j
            WRITE(idFile,"(10X, A)") '</DataItem>'
            WRITE(idFile,"( 8X, A)") '</Attribute>'
         END DO


         WRITE(idFile,"( 6X,A)") '</Grid>'
      END DO
      ! fechamento do arquivo
      WRITE(idFile,"( 4X,A)") '</Grid>'
      WRITE(idFile,"( 2X,A)") '</Domain>'
      WRITE(idFile,"(    A)") '</Xdmf>'

      CLOSE(idFile)
      RETURN  
   END SUBROUTINE

   !=======================================================================
   SUBROUTINE xdmf_basic_structure(var_names, var_types, num_names)
      ! Rotina para criação do arquivo XDML. 
      ! Cria um arquivo xdml com as seguintes propriedades:
      ! - a malha do arquivo não se altera para cada incremento de tempo.
      ! - a malha é composta por apenas 1 tipo de elemento.
      ! - os resultados são dados para as células ou para os nós.
      ! CHARACTER(LEN=*), INTENT(IN) :: filename
      CHARACTER(LEN=*), DIMENSION(*), INTENT(IN) :: var_names
      INTEGER, INTENT(IN) :: num_names
      INTEGER, DIMENSION(*), INTENT(IN) :: var_types

      CHARACTER(LEN=6) :: passo_tempo, posit_data, data_valueA, data_valueB
      CHARACTER(LEN=25) :: type_data
      INTEGER :: idFile, i, j, div
      INTEGER :: npassos =10,nelm =60,nnos=42,nnoelm =3, ndim=2 

      idFile=2000

      OPEN(UNIT=idFile,FILE=TRIM(filename)//'.xdmf',STATUS='replace')

      ! Escrita do cabeçalho do arquivo.
      WRITE(idFile,"(     A)") '<?xml version="1.0"?>'
      WRITE(idFile,"(     A)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      WRITE(idFile,"(     A)") '<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
      WRITE(idFile,"( 2X, A)") '<Domain>'

      ! Inicialização da estrutura da série temporal
      WRITE(idFile,"( 4X, A)") '<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
      ! Listagem dos tempos
      WRITE(idFile,"( 6X, A)") '<Time TimeType="List">'
      WRITE(idFile,"( 8X,A,I5,A)") '<DataItem Format="HDF" Dimensions="',npassos,'">'
      WRITE(idFile,"(10X,2A)") TRIM(filename), '.h5:/Time/TimeSerie'
      WRITE(idFile,"( 8X, A)") '</DataItem>'
      WRITE(idFile,"( 6X, A)") '</Time>'

      ! Listagem da estrutura - 1 para cada passo de tempo
      DO i=1,npassos
         WRITE(passo_tempo, "(I6)") i-1
         passo_tempo=adjustl(passo_tempo)
         WRITE(idFile,"( 6X, A)") '<Grid Name="mesh" GridType="Uniform">'

         ! Listagem da topologia
         WRITE(idFile,"( 8X,A,I5,A,I5,A)") '<Topology NumberOfElements="',nelm,&
                  '" TopologyType="Triangle" NodesPerElement="',nnoelm,'">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nelm,nnoelm,'" NumberType="UInt" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Mesh/topology'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Topology>'

         ! Listagem da geometria
         WRITE(idFile,"( 8X,A)") '<Geometry GeometryType="XY">'
         WRITE(idFile,"(10X,A,2I5,A)") '<DataItem Dimensions="',nnos,ndim,'" Format="HDF">'
         WRITE(idFile,"(12X,2A)") TRIM(filename), '.h5:/Mesh/geometry'
         WRITE(idFile,"(10X, A)") '</DataItem>'
         WRITE(idFile,"( 8X, A)") '</Geometry>'

         ! Listagem do atributo - repitido para cada atributo
         DO j=1,num_names
            ! Teste para verificar a posição que propriedade esta armazenada 
            div=ABS(var_types(j)-1)/5
            SELECT CASE (div)
            ! Atributo nos nós
            CASE (0)
               WRITE(data_valueA,"(I6)") nnos
               WRITE(posit_data,"(A)") 'Node'
            ! Atributo nas celulas
            CASE (1)
               WRITE(data_valueA,"(I6)") nelm
               WRITE(posit_data,"(A)") 'Cell'
            ! Atributo nas arrestas - Em Construção
            CASE (2)
               WRITE(data_valueA,"(I6)") 2
               WRITE(posit_data,"(A)") 'Edge'
            ! Atributo nas faces - Em Construção
            CASE (3)
               WRITE(data_valueA,"(I6)") 2
               WRITE(posit_data,"(A)") 'Face'
            ! Atributo na grade - Em Construção
            CASE (4)
               WRITE(data_valueA,"(I6)") 2
               WRITE(posit_data,"(A)") 'Grid'
            ! Atributo em outra posição - consultar manual - Em Construção
            CASE (5)
               WRITE(data_valueA,"(I6)") 2
               WRITE(posit_data,"(A)") 'Other'
            ! ERROR
            CASE DEFAULT
            END SELECT

            ! Teste para verificar a dimensão da propriedade armazenada 
            SELECT CASE (MOD(ABS(var_types(j)),5))
            ! Atributo escalar
            CASE (1)
               WRITE(data_valueB,"(I6)") 1
               WRITE(type_data,"(A)") '" AttributeType="Scalar"'
            ! Atributo vetorial
            CASE (2)
               WRITE(data_valueB,"(I6)") 3
               WRITE(type_data,"(A)") '" AttributeType="Vector"'
            ! Atributo tensorial
            CASE (3)
               WRITE(data_valueB,"(I6)") 9
               WRITE(type_data,"(A)") '" AttributeType="Tensor"'
            ! Atributo tensorial simétrico
            CASE (4)
               WRITE(data_valueB,"(I6)") 6
               WRITE(type_data,"(A)") '" AttributeType="Tensor6"'
            ! Atributo matriz genérica
            CASE (0)
               WRITE(data_valueB,"(I6)") 1
               WRITE(type_data,"(A)") '" AttributeType="Matrix"'
            ! ERROR
            CASE DEFAULT
            END SELECT

            ! Escrita dos dados
            WRITE(idFile,"( 8X,6A)") '<Attribute Name="',TRIM(var_names(j)), &
                                 TRIM(type_data),' Center="', &
                                 TRIM(posit_data),'">'
            WRITE(idFile,"(10X,4A)") '<DataItem Dimensions="',TRIM(data_valueA), &
                                 TRIM(data_valueB),'" Format="HDF">'
            ! Teste para verificar se a propriedade varia com o tempo
            IF (var_types(j).LT.0) THEN
               WRITE(idFile,"(12X,5A)") TRIM(filename), '.h5:/Results/', &
                                    TRIM(posit_data),'/',TRIM(var_names(j))
            ELSE
               WRITE(idFile,"(12X,7A)") TRIM(filename), '.h5:/Results/', &
                                    TRIM(posit_data),'/',TRIM(var_names(j)), &
                                    '/', TRIM(passo_tempo)
            END IF

            WRITE(idFile,"(10X, A)") '</DataItem>'
            WRITE(idFile,"( 8X, A)") '</Attribute>'
         END DO

         WRITE(idFile,"( 6X,A)") '</Grid>'
      END DO

      ! fechamento do arquivo
      WRITE(idFile,"( 4X,A)") '</Grid>'
      WRITE(idFile,"( 2X,A)") '</Domain>'
      WRITE(idFile,"(    A)") '</Xdmf>'

      CLOSE(idFile)
      RETURN  
   END SUBROUTINE

end module



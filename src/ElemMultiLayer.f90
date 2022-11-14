!=================================================================================
!
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Castro,        ecastro@lncc.br
!
!         LNCC/MCT
!         Petropolis, 03.2020
!=================================================================================
module mMLayer

   implicit none

   INTEGER :: nTypesMLayers
   INTEGER, POINTER :: LayerMaterial(:,:)=>null()

   PUBLIC   :: nTypesMLayers, LayerMaterial
   PUBLIC   :: MLayer_Inic, MLayer_Finalizar

   contains

   !=======================================================================
   subroutine MLayer_Inic(ndof)
      ! Rotina para alocar a memória para os ponteiros dos materiais
      INTEGER, INTENT(IN):: ndof
      ALLOCATE(LayerMaterial(nTypesMLayers,ndof))
   end subroutine

   !=======================================================================
   subroutine MLayer_Finalizar()
      ! Rotina para desalocar a memória consumida pelo modulo
      DEALLOCATE(LayerMaterial)
   end subroutine


end module
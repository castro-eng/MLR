!!>
!!         programa de elementos finitos em fortran 90
!!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!!
!!         Pilato Jr , Vinicius Pessoa  {pilato@deepsoft.com.br, pessoa@deepsoft.com.br}
!!
!!         Desenvolvido por DeepSoft para LNCC/MCT
!!         Rio de Janeiro, 11.2013-04.2014
!!
!!=================================================================================

module mInputReader
   !> Modulo responsavel por reunir subrotinas para leitura do arquivo de entrada.

   !> Armazena as linhas do arquivo de input.
   character(len=200), allocatable :: file_lines(:)
   !> Armazena o numero de linhas no arquivo.
   integer*4 number_of_lines

   contains

   subroutine leituraValoresCondContornoDS(keyword_name, f,ndof,numnp,j, nlvect,iprtin)
      !> Leitura e geracao de valores de condicoes de contorno.
    
      !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
      !! @param f             f.
      !! @param ndof          ndof
      !! @param numnp         numnp
      !! @param j             j
      !! @param nlvect        nlvect
      !! @param iprtin        iprtin

      use mleituraEscrita, only: iecho, printf, printd
      !  use mInputReader,    only: findKeyword, genflDS

      implicit none

      integer*4 :: ndof, numnp, j, nlvect, iprtin, keyword_line
      character(len=50) :: keyword_name
      real*8 f(ndof,numnp,nlvect)

      logical lzero
      integer*4 nlv
      character(len=35) :: rotulo

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.-1) return

      f = 0d0

      do 100 nlv=1,nlvect
         ! call genflDS(f(1,1,nlv),ndof,keyword_line)
         call genflDS_ptj(f(1,1,nlv),ndof,keyword_line)
         call ztest(f(1,1,nlv),ndof*numnp,lzero)
         if (iprtin.eq.0) then
            if (lzero) then
               if (j.eq.0) write(iecho,1000) nlv
               if (j.eq.1) write(iecho,2000)
            else
               if (j.eq.0) call printf(f,ndof,numnp,nlv)
               if (j.eq.1) then
                  rotulo=" n o d a l  b o d y  f o r c e s"
                  call printd (rotulo, f,ndof,numnp,iecho)
               end if
            endif
         endif
      100 continue
      return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
   end subroutine leituraValoresCondContornoDS
   
   !**********************************************************************

   subroutine leituraValoresCondContornoDS_PressaoFratura(keyword_name, f,ndof,numnp,j, nlvect,iprtin)
      use mleituraEscrita, only: iecho, printf, printd
      use mMalha, only: nenFratura, numelBordas, numelFratura, conecNodaisFratura, conecNodaisBordas
      ! use mInputReader,    only: findKeyword, genflDS

      implicit none

      integer*4 :: ndof, numnp, j, nlvect, iprtin, keyword_line
      character(len=50) :: keyword_name
      real*8 f(ndof,numnp,nlvect)
      integer :: contador(numnp)

      logical lzero
      integer*4 nlv, i, n, nel, nelBorda
      character(len=35) :: rotulo

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.-1) return

      f = 0d0

      do 100 nlv=1,nlvect
         ! call genflDS(f(1,1,nlv),ndof,keyword_line)
         call genflDS_ptj(f(1,1,nlv),ndof,keyword_line)

         call ztest(f(1,1,nlv),ndof*numnp,lzero)
                        
         if (iprtin.eq.0) then
            if (lzero) then
               if (j.eq.0) write(iecho,1000) nlv
               if (j.eq.1) write(iecho,2000)
            else
               if (j.eq.0) call printf(f,ndof,numnp,nlv)
               if (j.eq.1) then
                  rotulo=" n o d a l  b o d y  f o r c e s"
                  call printd (rotulo, f,ndof,numnp,iecho)
               end if
            endif
         endif
      100 continue

      do i=2,ndof
         f(i,:,:)=f(1,:,:)
      end do

      !........numero de graus de liberdade        
      ! f(2,:,:)=f(1,:,:)
      ! f(3,:,:)=f(1,:,:)
      
      ! Versão alternativa
      ! do nel=1,numelFratura
      !    do i=1,nenFratura
      !       n = conecNodaisFratura(i,nel)
      !       ! Verificando se o no atual e uma extremidade
      !       if (contador(n) == 1) then
      !          ! Verificando se a extremidade esta na borda
      !          do nelBorda = 1, numelBordas
      !             ! Todos os nos de borda aparecem duas vezes em conecNodaisBordas,
      !             ! por isso so precisamos olhar o primeiro no de cada elemento de borda
      !             if (n == conecNodaisBordas(1,nelBorda)) then
      !                ! Graus adicionais possuem o mesmo id que seus equivalentes
      !                f(2,:,:)=f(1,:,:)
      !                f(3,:,:)=f(1,:,:)
      !                ! print*, "aqui1", nel, n, id(1:3,n)
      !             end if
      !          end do
      !       else
      !          ! Graus adicionais possuem o mesmo id que seus equivalentes
      !          f(2,:,:)=f(1,:,:)
      !          f(3,:,:)=f(1,:,:)
      !          ! print*, "aqui2", nel, n, id(1:3,n)
      !       end if
      !    end do
      ! end do

      return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
   end subroutine leituraValoresCondContornoDS_PressaoFratura
   
   !**********************************************************************
   subroutine leituraValoresIniciaisDS(keyword_name, campo, ndof, numnp)
      !> Leitura e geracao de valores iniciais
    
      !! @param keyword_name  Palavra-chave que indica o tipo de valores iniciais.
      !! @param campo         Campo que recebera os valores iniciais.
      !! @param ndof          Numero de graus de liberdade do campo
      !! @param numnp         Numero de nos

      implicit none

      integer*4,intent(in) :: ndof, numnp
      character(len=50), intent(in) :: keyword_name
      real*8, intent(out) :: campo(ndof,numnp)

      integer*4 :: keyword_line

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.0) return

      call genflDS_ptj(campo(:,:),ndof,keyword_line)

   end subroutine leituraValoresIniciaisDS

   !**********************************************************************

   subroutine leituraCodigosCondContornoFluxoNormalDS(id, ndof, idP, ndofP, numnp, neq)
      !> Leitura e geracao de codigos de condicoes de contorno.
      !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
      !! @param id            id.
      !! @param ndof          ndof
      !! @param idP           idP
      !! @param ndofP         ndofP
      !! @param numnp         numnp
      !! @param neq           neq

      integer*4, intent(in) :: ndof, ndofP, numnp, idP(ndofP,numnp)
      integer*4 ::  neq
      integer*4, intent(inout) :: id(ndof,numnp)

      integer*4:: n, i

      do i=1, numnp
         if(idP(1,i).ne.0) then
            id(1,i)=1
         else
            id(1,i)=0
         endif
      end do

      !.... establish equation numbers
      neq = 0
      do 400 n=1,numnp
         do 300 i=1,ndof
            if (id(i,n).eq.0) then
               neq = neq + 1
               id(i,n) = neq
            else
               id(i,n) = 1 - id(i,n)
            endif
         300 continue
      400 continue
      return
   end subroutine leituraCodigosCondContornoFluxoNormalDS 
      
   !**********************************************************************

   subroutine leituraCodigosCondContornoDS(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
      !> Leitura e geracao de codigos de condicoes de contorno.
    
      !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
      !! @param id            id.
      !! @param ndof          numnp
      !! @param numnp         numnp
      !! @param neq           neq
      !! @param iecho         iecho
      !! @param iprtin        iprtin

      integer*4, intent(in) :: ndof, numnp, iecho, iprtin
      integer*4 ::  neq, keyword_line
      integer*4, intent(inout) :: id(ndof,numnp)
      character(len=50) keyword_name

      integer*4:: nn, n, i
      logical pflag

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.-1) return

      id = 0
      call igenDS(id,ndof, keyword_line)

      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
            pflag = .false.

            do 100 i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
            100    continue

            if (pflag) then
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
               write(iecho,2000) n,(id(i,n),i=1,ndof)
            endif
         200    continue
      endif

      !.... establish equation numbers

      neq = 0

      do 400 n=1,numnp
         do 300 i=1,ndof
            if (id(i,n).eq.0) then
               neq = neq + 1
               id(i,n) = neq
            else
               id(i,n) = 1 - id(i,n)
            endif
         300 continue
      400 continue

      return

 1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n &
     &         c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
   end subroutine leituraCodigosCondContornoDS
   
   !********************************************************************

   subroutine leituraCodigosCondContornoDS_PressaoFratura(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
      !> Leitura e geracao de codigos de condicoes de contorno para deslocamento com fraturas
      
      !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
      !! @param id            id.
      !! @param ndof          numnp
      !! @param numnp         numnp
      !! @param neq           neq
      !! @param iecho         iecho
      !! @param iprtin        iprtin
      
      use mMalha, only: nsd,numel,conecNodaisFratura, conecNodaisBordas, numelFratura, nenFratura, numelBordas

      implicit none

      integer*4, intent(in) :: ndof, numnp, iecho, iprtin
      integer*4 ::  neq, keyword_line
      integer*4, intent(inout) :: id(ndof,numnp)
      character(len=50) keyword_name

      integer*4:: nn, n, i, j, nel, nelBorda
      integer*4 :: contador(numnp)
      logical :: pflag
      character(*), parameter :: fmtTitulo = "('1',' n o d a l   b o u n d a r y   c o n d i t i o n          c o  d e s' /// &
               & 5x,' node no.',3x,6(6x,'dof',i1:)//)"

      keyword_line = findKeyword(keyword_name)
      if (keyword_line == 0) return

      id = 0
      call igenDS(id, ndof, keyword_line)

      ! Verifica quais nos sao extremidades das fraturas:
      ! esses nos aparecem apenas uma vez nas conectividades de fratura
      contador(:) = 0
      do nel=1,numelFratura
         do i=1,nenFratura
            n = conecNodaisFratura(i,nel)
            contador(n) = contador(n) + 1
         end do
      end do

      ! prende os graus de liberdade extras
      do i=2,ndof
         id(i,:) = 1
      end do
         
      ! ! Todos os nos de fratura possuem o mesmo id que seus equivalentes...
      ! do nel=1, numelFratura
      !    do i=1, nenFratura
      !       n = conecNodaisFratura(i,nel)
      !       id(2,n) = id(1,n)
      !       id(3,n) = id(1,n)
      !    end do
      ! end do
         
      ! Versão alternativa
      do nel=1,numelFratura
         do i=1,nenFratura
            n = conecNodaisFratura(i,nel)
            ! Verificando se o no atual e uma extremidade
            if (contador(n) == 1) then
               ! Verificando se a extremidade esta na borda
               do nelBorda = 1, numelBordas
                  ! Todos os nos de borda aparecem duas vezes em conecNodaisBordas,
                  ! por isso so precisamos olhar o primeiro no de cada elemento de borda
                  if (n == conecNodaisBordas(1,nelBorda)) then
                     ! Graus adicionais possuem o mesmo id que seus equivalentes
                     do j=2,ndof
                        id(j,n) = id(1,n)
                     end do
                     ! id(3,n) = id(1,n)
                     ! print*, "aqui1", nel, n, id(1:3,n)
                  end if
               end do
            else
               ! Graus adicionais possuem o mesmo id que seus equivalentes
               do j=2,ndof
                  id(j,n) = id(1,n)
               end do
               ! id(2,n) = id(1,n)
               ! id(3,n) = id(1,n)
               ! print*, "aqui2", nel, n, id(1:3,n)
            end if
         end do
      end do
         
      ! do nel=1,numelFratura
      !    do i=1,nenFratura
      !       n = conecNodaisFratura(i,nel)
      !       print*, n, id(:, n)
      !    end do
      ! end do 

      if (iprtin.eq.0) then
         nn=0
         do n=1,numnp
            pflag = .false.
               
            do i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
            end do
               
            if (pflag) then
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho, fmtTitulo) (i,i=1,ndof)
               write(iecho, "(6x,i10,5x,6(5x,i10))") n,(id(i,n),i=1,ndof)
            endif
         end do
      endif
         
      !.... establish equation numbers
      neq = 0
      do n=1,numnp
         do i=1,ndof
            if (id(i,n).eq.0) then
               neq = neq + 1
               id(i,n) = neq
            else
               id(i,n) = 1 - id(i,n)
            endif
         end do
      end do

   end subroutine leituraCodigosCondContornoDS_PressaoFratura
   
   !********************************************************************

   subroutine leituraCodigosCondContornoDS_PressaoFratura_3D(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
      use mMalha, only: conecNodaisFratura, conecNodaisBordas, conectividadesInterseccao
      use mMalha, only: nsd, numel, numelFratura, nenFratura, numelBordas

      implicit none

      integer*4, intent(in) :: ndof, numnp, iecho, iprtin
      integer*4 ::  neq, keyword_line
      integer*4, intent(inout) :: id(ndof,numnp)
      character(len=50) keyword_name

      integer*4:: nn, n, i, j, nel, nelBorda
      integer*4 :: contador(numnp),noInterseccao
      logical :: pflag
      character(*), parameter :: fmtTitulo = "('1',' n o d a l   b o u n d a r y   c o n d i t i o n          c o  d e s' /// &
               & 5x,' node no.',3x,6(6x,'dof',i1:)//)"

      keyword_line = findKeyword(keyword_name)
      if (keyword_line == 0) return

      id = 0
      call igenDS(id, ndof, keyword_line)

      id(2,:) = 1
      id(3,:) = 1
         
      ! Todos os nos de fratura possuem o mesmo id que seus equivalentes...
      do nel=1, numelFratura
         do i=1, nenFratura
            n = conecNodaisFratura(i,nel)
            id(2,n) = id(1,n)
            id(3,n) = id(1,n)
         end do
      end do

      !...exceto as bordas que devem estar presas
      ! este trecho de código conta o numero de contribuições do nó. se for menor que 4, sabemos que 
      ! está em uma das bordas.
      ! localiza as bordas das fraturas... se contador(n)==1 ou contador(n)==2 estou na borda
      ! contador(:) = 0
      ! do nel=1,numelFratura
      !    do i=1,nenFratura
      !       n = conecNodaisFratura(i,nel)
      !       contador(n) = contador(n) + 1
      !    end do
      ! end do
         
      ! do nel=1,numelFratura
      !    do i=1,nenFratura
      !       n = conecNodaisFratura(i,nel)
      !       if(contador(n)==1 .or. contador(n)==2) then
      !          id(2,n) = 1
      !          id(3,n) = 1
      !       endif
      !    end do
      ! end do
         
      ! porem se a fratura encostar na borda, ela pode ser aberta.
      ! do nelBorda = 1, numelBordas
      !    do i=1, 2
      !       if(conectividadesInterseccao(i,nelBorda)>0)then
      !          noInterseccao=conectividadesInterseccao(i,nelBorda)
      !          id(2,noInterseccao)=id(1,noInterseccao)
      !          id(3,noInterseccao)=id(1,noInterseccao)
      !       endif
      !    enddo
      ! enddo

      ! do nel=1,numelFratura
      !    do i=1,nenFratura
      !       n = conecNodaisFratura(i,nel)
      !       print*, n, id(:, n)
      !    end do
      ! end do

      if (iprtin.eq.0) then
         nn=0
         do n=1,numnp
            pflag = .false.

            do i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
            end do

            if (pflag) then
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho, fmtTitulo) (i,i=1,ndof)
               write(iecho, "(6x,i10,5x,6(5x,i10))") n,(id(i,n),i=1,ndof)
            endif
         end do
      endif

      !.... establish equation numbers
      neq = 0
      do n=1,numnp
         do i=1,ndof
            if (id(i,n).eq.0) then
               neq = neq + 1
               id(i,n) = neq
            else
               id(i,n) = 1 - id(i,n)
            endif
         end do
      end do

   end subroutine leituraCodigosCondContornoDS_PressaoFratura_3D
   
   !********************************************************************

   subroutine leituraCodigosCondContornoDS_PressaoDescontinuidade(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
      !> Leitura e geracao de codigos de condicoes de contorno para deslocamento com descontinuidades
      !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
      !! @param id            id.
      !! @param ndof          numnp
      !! @param numnp         numnp
      !! @param neq           neq
      !! @param iecho         iecho
      !! @param iprtin        iprtin

      use mMalha,    only: numelFratura, nenFratura, numelBordas
      use mMalha,    only: conecNodaisFratura, conecNodaisBordas
      use mMLayer,   only: LayerMaterial
      use mGlobaisArranjos,  only: matFratura

      implicit none

      integer*4, intent(in) :: ndof, numnp, iecho, iprtin
      integer*4 ::  neq, keyword_line
      integer*4, intent(inout) :: id(ndof,numnp)
      character(len=50) keyword_name

      integer*4:: nn, n, i, j, k, nel, nelBorda, nlayers
      integer*4 :: contador(numnp)
      logical :: pflag
      character(*), parameter :: fmtTitulo = "('1',' n o d a l   b o u n d a r y   c o n d i t i o n          c o  d e s' /// &
               & 5x,' node no.',3x,6(6x,'dof',i1:)//)"

      keyword_line = findKeyword(keyword_name)
      if (keyword_line == 0) return

      id = 0
      call igenDS(id, ndof, keyword_line)

      ! Verifica quais nos sao extremidades das fraturas:
      ! esses nos aparecem apenas uma vez nas conectividades de fratura
      contador(:) = 0
      do nel=1,numelFratura
         do i=1,nenFratura
            n = conecNodaisFratura(i,nel)
            contador(n) = contador(n) + 1
         end do
      end do

      ! prende os graus de liberdade extras
      do i=2,ndof
         id(i,:) = 1
      end do

      
      ! Versão alternativa
      do nel=1,numelFratura
         nlayers = LayerMaterial(matFratura(nel),1) + 2
         do i=1,nenFratura
            n = conecNodaisFratura(i,nel)
            ! Verificando se o no atual e uma extremidade
            if (contador(n) == 1) then
               ! Verificando se a extremidade esta na borda
               if ((n.eq.3).or.(n.eq.4))then
                  print*,"id(:,n)", id(:,n)
               endif
               do nelBorda = 1, numelBordas
                  ! Todos os nos de borda aparecem duas vezes em conecNodaisBordas,
                  ! por isso so precisamos olhar o primeiro no de cada elemento de borda
                  do k=1,nenFratura
                     if (n == conecNodaisBordas(k,nelBorda)) then
                        ! Graus adicionais possuem o mesmo id que seus equivalentes
                        do j=2,nlayers
                           id(j,n) = id(1,n)
                        end do
                        exit
                     end if
                  enddo
               end do
            else
               ! Graus adicionais possuem o mesmo id que seus equivalentes
               do j=2,nlayers
                  id(j,n) = id(1,n)
               end do
            end if
            if ((n.eq.3).or.(n.eq.4))then
               print*,"id(:,n)", id(:,n)
            endif
         end do
      end do

      if (iprtin.eq.0) then
         nn=0
         do n=1,numnp
            pflag = .false.
            do i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
            end do

            if (pflag) then
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho, fmtTitulo) (i,i=1,ndof)
               write(iecho, "(6x,i10,5x,6(5x,i10))") n,(id(i,n),i=1,ndof)
            endif
         end do
      endif

      !.... establish equation numbers
      neq = 0
      do n=1,numnp
         do i=1,ndof
            if (id(i,n).eq.0) then
               neq = neq + 1
               id(i,n) = neq
            else
               id(i,n) = 1 - id(i,n)
            endif
         end do
      end do

   end subroutine leituraCodigosCondContornoDS_PressaoDescontinuidade
   
   !********************************************************************

   subroutine leituraCodigosCondContornoDS_PressaoDescontinuidade_3D(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
      use mMalha, only: conecNodaisFratura, conecNodaisBordas, conectividadesInterseccao
      use mMalha, only: nsd, numel, numelFratura, nenFratura, numelBordas

      implicit none

      integer*4, intent(in) :: ndof, numnp, iecho, iprtin
      integer*4 ::  neq, keyword_line
      integer*4, intent(inout) :: id(ndof,numnp)
      character(len=50) keyword_name

      integer*4:: nn, n, i, j, nel, nelBorda
      integer*4 :: contador(numnp),noInterseccao
      logical :: pflag
      character(*), parameter :: fmtTitulo = "('1',' n o d a l   b o u n d a r y   c o n d i t i o n          c o  d e s' /// &
               & 5x,' node no.',3x,6(6x,'dof',i1:)//)"

      keyword_line = findKeyword(keyword_name)
      if (keyword_line == 0) return

      id = 0
      call igenDS(id, ndof, keyword_line)

      id(2,:) = 1
      id(3,:) = 1
         
      ! Todos os nos de fratura possuem o mesmo id que seus equivalentes...
      do nel=1, numelFratura
         do i=1, nenFratura
            n = conecNodaisFratura(i,nel)
            id(2,n) = id(1,n)
            id(3,n) = id(1,n)
         end do
      end do

      if (iprtin.eq.0) then
         nn=0
         do n=1,numnp
            pflag = .false.

            do i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
            end do

            if (pflag) then
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho, fmtTitulo) (i,i=1,ndof)
               write(iecho, "(6x,i10,5x,6(5x,i10))") n,(id(i,n),i=1,ndof)
            endif
         end do
      endif

      !.... establish equation numbers
      neq = 0
      do n=1,numnp
         do i=1,ndof
            if (id(i,n).eq.0) then
               neq = neq + 1
               id(i,n) = neq
            else
               id(i,n) = 1 - id(i,n)
            endif
         end do
      end do

   end subroutine leituraCodigosCondContornoDS_PressaoDescontinuidade_3D

   !********************************************************************

   subroutine leituraCodigosCondContornoDS_PressaoFratura_old(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
      use mMalha, only: nsd,numel,conecNodaisFratura, conecNodaisBordas, numelFratura, nenFratura, numelBordas

      implicit none

      integer*4, intent(in) :: ndof, numnp, iecho, iprtin
      integer*4 ::  neq, keyword_line
      integer*4, intent(inout) :: id(ndof,numnp)
      character(len=50) keyword_name

      integer*4:: nn, n, i, j, nel, nelBorda
      integer*4 :: contador(numnp)
      logical :: pflag
      character(*), parameter :: fmtTitulo = "('1',' n o d a l   b o u n d a r y   c o n d i t i o n          c o  d e s' /// &
               & 5x,' node no.',3x,6(6x,'dof',i1:)//)"

      keyword_line = findKeyword(keyword_name)
      if (keyword_line == 0) return

      id = 0
      call igenDS(id, ndof, keyword_line)

      ! Verifica quais nos sao extremidades das fraturas:
      ! esses nos aparecem apenas uma vez nas conectividades de fratura
      contador(:) = 0
      do nel=1,numelFratura
         do i=1,nenFratura
            n = conecNodaisFratura(i,nel)
            contador(n) = contador(n) + 1
         end do
      end do

      id(2,:) = 1
      id(3,:) = 1

      ! Versão alternativa
      do nel=1,numelFratura
         do i=1,nenFratura
            n = conecNodaisFratura(i,nel)
            print*, n
            ! Verificando se o no atual e uma extremidade
            if (contador(n) == 1) then
               ! Verificando se a extremidade esta na borda
               do nelBorda = 1, numelBordas
                  ! Todos os nos de borda aparecem duas vezes em conecNodaisBordas,
                  ! por isso so precisamos olhar o primeiro no de cada elemento de borda
                  if (n == conecNodaisBordas(1,nelBorda)) then
                     ! Graus adicionais possuem o mesmo id que seus equivalentes
                     id(2,n) = id(1,n) ! quando está na borda e é extremidadade da fratura
                     id(3,n) = id(1,n) 
                  end if
                  if (nsd==3 .and. n == conecNodaisBordas(3,nelBorda)) then  !Gambiarra da tuane
                     ! Graus adicionais possuem o mesmo id que seus equivalentes
                     id(2,n) = id(1,n) ! quando está na borda e é extremidadade da fratura
                     id(3,n) = id(1,n) 
                  end if
               end do
            else
               ! Graus adicionais possuem o mesmo id que seus equivalentes
               id(2,n) = id(1,n) !Quando o nó está no meio da fratura
               id(3,n) = id(1,n)
            end if
         end do
      end do

      if (iprtin.eq.0) then
         nn=0
         do n=1,numnp
            pflag = .false.

            do i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
            end do

            if (pflag) then
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho, fmtTitulo) (i,i=1,ndof)
               write(iecho, "(6x,i10,5x,6(5x,i10))") n,(id(i,n),i=1,ndof)
            endif
         end do
      endif

      !.... establish equation numbers

      neq = 0
      do n=1,numnp
         do i=1,ndof
            if (id(i,n).eq.0) then
               neq = neq + 1
               id(i,n) = neq
            else
               id(i,n) = 1 - id(i,n)
            endif
         end do
      end do

   end subroutine leituraCodigosCondContornoDS_PressaoFratura_old

   !********************************************************************

   subroutine leituraGeracaoCoordenadasDS(x, nsd, numnp, icoords, iprtin)
      !> Leitura e geracao de coordenadas
    
      !! @param x          Matriz onde serao armazeadas as coordenadas.
      !! @param nsd        Corresponde ao nao de linhas da matriz x.
      !! @param numnp      Numero de xxx
      !! @param icoords    Handle para o arquivo de coordenadas.
      !! @param iprtin     iprtin

      implicit none

      integer*4, intent(in)   :: nsd, numnp, icoords, iprtin
      real*8, intent(inout) ::  x(nsd,*)
      integer*4:: i, n, keyword_line
      character(len=50) keyword_name

      keyword_name = 'coordenadas_nodais'
      keyword_line = findKeyword(keyword_name)

      call genflDS(x,nsd,keyword_line)

      if (iprtin.eq.1) return
      write(icoords,*) "# Coordenadas ", nsd
      do n=1,numnp
        write(icoords,2000) n,(x(i,n),i=1,nsd)
      end do

      return
      2000 format(6x,i12,10x,3(1pe15.8,2x))
   end subroutine leituraGeracaoCoordenadasDS
   
   !********************************************************************

   subroutine readInputFileDS()
      !> Le arquivo de input e armazena seu conteudo em um array.
      !! @param file_name Nome do arquivo a ser lido.
      use mLeituraEscrita,   only: iin

      implicit none

      integer*4 success, lines_count
      character(len=200) file_line

      integer*4 :: main_number_of_lines, main_number_of_includes, i, include_index, inc_nlines, inc_inc, merge_lines
      character(len=200), allocatable :: main_file_lines(:)
      character(len=200), allocatable :: include_files(:)
      character(len=200) include_file
      integer*4, allocatable :: include_indexes(:), include_number_of_lines(:)

      main_number_of_lines = 0
      main_number_of_includes = 0

      call analyzeFileInput(main_number_of_lines, main_number_of_includes)

      if (main_number_of_includes.eq.0) then
         call createSimpleInputFile()
         return
      end if

      allocate(main_file_lines(main_number_of_lines))
      allocate(include_indexes(main_number_of_includes))
      allocate(include_number_of_lines(main_number_of_includes))
      allocate(include_files(main_number_of_includes))

      lines_count = 1
      do
         read(iin, "(A)", iostat=success) file_line
         if (success.ne.0) exit
         main_file_lines(lines_count) = file_line
         lines_count = lines_count + 1
      end do
      rewind(iin)

      number_of_lines = main_number_of_lines

      ! Number of lines
      do i=1, main_number_of_includes
         include_index = findInclude(i, main_file_lines, main_number_of_lines)
         read(main_file_lines(include_index), '(A)') include_file

         include_file = adJustl(include_file)
         call analyzeFile(include_file, inc_nlines, inc_inc)
         number_of_lines = number_of_lines + inc_nlines
         include_indexes(i) = include_index
         include_number_of_lines(i) = inc_nlines
         include_files(i) = include_file
      end do

      ! Prepare final struct.
      call prepareFileLines(include_indexes, include_number_of_lines, main_number_of_includes, main_file_lines)

      ! Merge contensts.
      merge_lines = 0
      do i=1, main_number_of_includes
         call mergeIncludeContents(include_files(i), include_indexes(i) + merge_lines)
         merge_lines = merge_lines + include_number_of_lines(i)
      end do

      deallocate(main_file_lines)
      deallocate(include_indexes)
      deallocate(include_number_of_lines)
      return
   end subroutine readInputFileDS
   
   !********************************************************************

   subroutine createSimpleInputFile()
      !> Cria a estretura de input usando um arquivo de entrada sem includes
      !! @param file_name Nome do arquivo a ser lido.
    
      use mLeituraEscrita,   only: iin

      implicit none

      integer*4 success, lines_count
      character(len=200) file_line

      number_of_lines = 0

      do
         read(iin, "(A)", iostat=success) file_line
         if (success.ne.0) exit
         number_of_lines = number_of_lines + 1
      end do
      rewind(iin)

      allocate(file_lines(number_of_lines))
      ! TO-DO avoid two-times read
      lines_count = 1
      do
         read(iin, "(A)", iostat=success) file_line
         if (success.ne.0) exit
         file_lines(lines_count) = file_line
         lines_count = lines_count + 1
      end do
      rewind(iin)
   end subroutine createSimpleInputFile
   
   !********************************************************************

   subroutine mergeIncludeContents(include_file, include_line)
      !> Le o conteudo do arquivo de include e armazena no array principal.
      !! @param   include_index   O index do include.
      !! @param   include_files   Array com includes.
      !! @param   include_line    A linha do include.

      implicit none

      integer*4 include_line
      character(len=200) include_file

      character(len=200) file_line
      integer*4 file_channel, success, current_index

      file_channel = 1

      current_index = include_line

      open(unit=file_channel, file=include_file)
      do
         read(file_channel, "(A)", iostat=success) file_line
         if (success.ne.0) exit
         file_lines(current_index) = file_line
         current_index = current_index + 1
      end do
      close(file_channel)
   end subroutine mergeIncludeContents
   
   !********************************************************************

   subroutine prepareFileLines(include_indexes, include_number_of_lines, number_of_includes, original_file_lines)
      !> Efetua a aloca��o da estrutura definitiva, preparando a linha dos arquivos originais para receber os includes
      !! @param   include_indexes             Array os indices de ocorrencias dos includes.
      !! @param   include_number_of_lines     Array com o numero de linhas de cada include
      !! @param   number_of_includes          Numero de includes.
      !! @param   original_file_lines         Linhas do arquivo de entrada original.

      integer*4 number_of_includes, number_of_original_lines, line_index, shift_lines
      integer*4 include_indexes(:), include_number_of_lines(:)
      character(len=200) original_file_lines(:)

      integer*4 current_include_index, original_index

      allocate(file_lines(number_of_lines))

      current_include_index = 1
      original_index = 1
      line_index = 1
      shift_lines = 0
      do while ( line_index <= number_of_lines)
         if (original_index.eq.(include_indexes(current_include_index))) then
            line_index = line_index + include_number_of_lines(current_include_index)
            current_include_index = current_include_index + 1
         end if
         file_lines(line_index) = original_file_lines(original_index)
         line_index = line_index + 1
         original_index = original_index + 1
      end do

   end subroutine prepareFileLines
   
   !********************************************************************

   subroutine analyzeFileInput(number_of_lines, number_of_includes)
      !> Efetua algumas análises no arquivo recebido.
      !! @param   number_of_lines     N�mero de linhas.
      !! @param   number_of_include   N�mero de ocorr�ncias da palavra include.
    
      use mLeituraEscrita,   only: iin

      character(len=200) file_line
      integer*4 number_of_lines, number_of_includes

      character(len=50) include_keyword, formated_keyword
      integer*4 keyword_len, success

      include_keyword = "include"
      keyword_len = len(trim(include_keyword)) + 2
      formated_keyword = trim('*' // trim(include_keyword) // '{')

      number_of_lines = 0
      number_of_includes = 0

      ! lunitInicial = 15
      do
         read(iin, "(A)", iostat=success) file_line
         if (success.ne.0) exit
         number_of_lines = number_of_lines + 1
         if (formated_keyword.eq.file_line(1:keyword_len)) then
            number_of_includes = number_of_includes + 1
         end if
      end do
      rewind(iin)

   end subroutine analyzeFileInput
   
   !********************************************************************

   subroutine analyzeFile(file_name, number_of_lines, number_of_includes)
      !> Efetua algumas an�lises no arquivo recebido.
      !! @param   file_name           O nome do arquivo.
      !! @param   number_of_lines     N�mero de linhas.
      !! @param   number_of_include   N�mero de ocorr�ncias da palavra include.
    
      character(len=200) file_name, file_line
      integer*4 number_of_lines, number_of_includes

      character(len=50) include_keyword, formated_keyword
      integer*4 keyword_len, file_channel, success

      include_keyword = "include"
      keyword_len = len(trim(include_keyword)) + 2
      formated_keyword = trim('*' // trim(include_keyword) // '{')

      number_of_lines = 0
      number_of_includes = 0

      file_channel = 2
      lunitInicial = 15
      file_channel = lunitInicial

      open(unit=file_channel, file=file_name)
      do
         read(file_channel, "(A)", iostat=success) file_line
         if (success.ne.0) exit
         number_of_lines = number_of_lines + 1
         if (formated_keyword.eq.file_line(1:keyword_len)) then
            number_of_includes = number_of_includes + 1
         end if
      end do
      close(file_channel)

   end subroutine analyzeFile
   
   !********************************************************************

   integer*4 function findInclude(position, file_lines, number_of_lines)
      !> Procura a n-esima palavra-chave include.
      !! @param  position         Corresponde a posicao desejada.
      !! @param  file_lines       Linhas do arquivo.
      !! @param  number_of_lines  Numero de linhas atuais.
      !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.

      implicit none
      integer*4 position, number_of_lines, current_position
      character(len=200) file_lines(:)
      character(50) keyword, formated_keyword
      character(len=120) file_line
      integer*4 i, keyword_len

      keyword = "include"
      keyword_len = len(trim(keyword)) + 2
      formated_keyword = trim('*' // trim(keyword) // '{')
      current_position = 0

      do i=1, number_of_lines
         file_line = file_lines(i)
         if (formated_keyword.eq.file_line(1:keyword_len)) then
            current_position = current_position + 1
            if (current_position.eq.position) then
               findInclude = i + 1
               return
            end if
         end if
      end do
      findInclude = 0
      return
   end function findInclude 
   
   !********************************************************************

   integer*4 function findKeyword(keyword)
      !> Procura uma palavra-chave.
      !! @param  keyword A palavra-chave.
      !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.

      implicit none
      character(50) keyword, formated_keyword
      character(len=120) file_line
      integer*4 i, keyword_len
      do i=1, number_of_lines, 1
         file_line = file_lines(i)
         keyword_len = len(trim(keyword)) + 2
         formated_keyword = trim('*' // trim(keyword) // '{')
         if (formated_keyword.eq.file_line(1:keyword_len)) then
            findKeyword = i + 1
            return
         end if
      end do
      findKeyword = 0
      return
   end function findKeyword
   
   !********************************************************************

   subroutine readIntegerKeywordValue(keyword, target, default_value)
      !> Efetua a leitura de uma palavra-chave to tipo inteiro. Se nao encontrado, associa o valor defualt fornecido.
      !! @param keyword       A palavra-chave a ser encontrada.
      !! @param target        Variavel onde o valor inteiro sera atribuido.
      !! @param default_value Valor default.
    
      implicit none
      character(50) keyword
      character(120) file_line
      integer*4 target, default_value, keyword_line
      keyword_line = findKeyword(keyword)
      if (keyword_line.eq.0) then
         target = default_value
         return
      end if
      file_line = adjustL(trim(file_lines(keyword_line)))
      read(file_line, *) target
      print*, keyword, target
      return
   end subroutine readIntegerKeywordValue
   
   !********************************************************************

   subroutine readStringKeywordValue(keyword, value, value_default)
      !> Efetua a leitura de uma palavra-chave to tipo string. Se nao encontrado, associa o valor defualt fornecido.
      !! @param keyword       A palavra-chave a ser encontrada.
      !! @param target        Variavel onde a string sera atribuido.
      !! @param default_value Valor default.
    
      implicit none
      character(50) keyword
      character(120) file_line
      character(len=*) :: value, value_default
      integer*4 keyword_line

      keyword_line = findKeyword(keyword)
      if (keyword_line.eq.0) then
         value=value_default
         return
      end if
      read(file_lines(keyword_line), '(a)') value
      return
   end subroutine readStringKeywordValue
   
   !********************************************************************

   subroutine readRealKeywordValue(keyword, target, default_value)
      !> Efetua a leitura de uma palavra-chave to tipo real. Se nao encontrado, associa o valor defualt fornecido.
      !! @param keyword       A palavra-chave a ser encontrada.
      !! @param target        Variavel onde o real sera atribuido.
      !! @param default_value Valor default.
    
      implicit none
      character(50) keyword
      character(120) file_line
      real(8) target,default_value
      integer*4 keyword_line
      keyword_line = findKeyword(keyword)
      if (keyword_line.eq.0) then
         target = default_value
         return
      end if
      file_line = adjustL(trim(file_lines(keyword_line)))
      read(file_line, *) target
      return
   end subroutine readRealKeywordValue
   
   !********************************************************************

   subroutine readLogicalKeywordValue(keyword, target, default_value)
      !! Efetua a leitura de uma palavra-chave to tipo logico. 
      !! Se nao encontrado, associa o valor defualt fornecido.
      !! @param keyword       A palavra-chave a ser encontrada.
      !! @param target        Variavel onde o real sera atribuido.
      !! @param default_value Valor default.
   
      implicit none
      character(50) keyword
      character(120) file_line
      logical target, default_value
      integer*4 keyword_line
      keyword_line = findKeyword(keyword)
      if (keyword_line.eq.0) then
         target = default_value
         return
      end if
      file_line = adjustL(trim(file_lines(keyword_line)))
      read(file_line, *) target
      return
   end subroutine readLogicalKeywordValue

   !********************************************************************

   subroutine readOutFlagKeyword(keyword, target, var_out, var_n)
      !> Efetua a leitura de valores definidos para propriedades de sa�da. Na prática são lidas 3 variáveis.
      !! @param keyword       A palavra-chave a ser encontrada.
      !! @param target        A variável destino.
      !! @param var_out       O valor da variável out.
      !! @param var_n         Valor n.
    
      implicit none
      character(50) keyword
      integer*4 target, var_n
      character(120) var_out
      character(120) file_line
      integer*4 keyword_line
      keyword_line = findKeyword(keyword)
      if (keyword_line.eq.0) then
         return
      end if
      read(file_lines(keyword_line), *) target
      keyword_line = keyword_line + 1
      read(file_lines(keyword_line), "(a)") var_out
      keyword_line = keyword_line + 1
      read(file_lines(keyword_line), *) var_n
   end subroutine readOutFlagKeyword

   !********************************************************************

   subroutine genflDS(a,nra, nLinhaArqInput )
      !> Efetua a geracao de coordeadas, de acordo com parametros.
    
      !! @param a      Matriz onde serao armazenados os dados.
      !! @param nra    Inteiro indicando nra
      !! @param nLinhaArqInput    Indice da linha onde as coordenadas estao posicionadas no array linhas no arquivo de entrada.
    
      use mMalha, only: genfl1   
      use mGlobaisEscalares

      implicit none

      integer*4:: nra, nLinhaArqInput
      real*8  :: a(nra,*)
      real*8  :: temp(6,20)
      integer*4:: n,numgp,ninc(3),inc(3)
      integer*4:: i, j, m, mgen

      100 continue
      read(file_lines(nLinhaArqInput),1000) n,numgp,(temp(i,1),i=1,nra)
      nLinhaArqInput = nLinhaArqInput + 1

      if (n.eq.0) return
      ! call move(a(1,n),temp,nra)
      a(1:nra,n) = temp(1:nra,1)

      if (numgp.ne.0) then
         do 200 j=2,numgp
            read(file_lines(nLinhaArqInput),1000) m,mgen,(temp(i,j),i=1,nra)
            nLinhaArqInput = nLinhaArqInput + 1

            if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m)         !B
         200    continue
         read(file_lines(nLinhaArqInput),2000) (ninc(i),inc(i),i=1,3)
         nLinhaArqInput = nLinhaArqInput + 1

         call genfl1(a,nra, temp, n, numgp, ninc, inc)
      endif
      go to 100

      1000 format(2i10,6f10.0)
      2000 format(16i10)

   end subroutine
   
   !********************************************************************

   subroutine genflDS_ptj(a,nra, nLinhaArqInput )
      use mMalha, only: genfl1   
      use mGlobaisEscalares

      implicit none

      integer*4:: nra, nLinhaArqInput
      real*8  :: a(nra,*)
      real*8  :: temp(6,20)
      integer*4:: n,numgp,ninc(3),inc(3)
      integer*4:: i, j, m, mgen

      100 continue
      read(file_lines(nLinhaArqInput),1000) n,numgp,(temp(i,1),i=1,nra)
      nLinhaArqInput = nLinhaArqInput + 1

      if (n.eq.0) return
      ! call move(a(1,n),temp,nra)
      a(1:nra,n) = temp(1:nra,1)

      if (numgp.ne.0) then
         do 200 j=2,numgp
            read(file_lines(nLinhaArqInput),1000) m,mgen,(temp(i,j),i=1,nra)
            nLinhaArqInput = nLinhaArqInput + 1

            if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m)         !B
         200    continue
         read(file_lines(nLinhaArqInput),2000) (ninc(i),inc(i),i=1,3)
         nLinhaArqInput = nLinhaArqInput + 1

         call genfl1(a,nra, temp, n, numgp, ninc, inc)
      endif
      go to 100

      1000 format(2i10,7f15.0) ! Limita a quantidade de dof para nlayer
      2000 format(16i10)

   end subroutine

   !********************************************************************

   subroutine igenDS(ia, m, nLinhaArqInput)
      !> Subrotina para ler e gerar dados nodais inteiros.
      !! @param ia              Array de entrada.
      !! @param m               Numero de linhas na matriz de entrada.
      !! @param nLinhaArqInput  Indice da linha no array de linhas do arquivo de entrada.
    
      use mGlobaisEscalares

      integer*4:: m, ia(m,*), nLinhaArqInput
      integer*4:: ib(m)
      integer*4:: n, ne, ng
      integer*4:: i
        
      100 continue
         read(file_lines(nLinhaArqInput),1000) n,ne,ng,(ib(i),i=1,m)
         nLinhaArqInput = nLinhaArqInput + 1

         if (n.eq.0) return

         if (ng.eq.0) then
            ne = n
            ng = 1
         else
            ne = ne - mod(ne-n,ng)
         endif

         do 200 i=n,ne,ng
            ia(:,i)=ib
         200 continue

      go to 100

      1000 format(16i10)
   end subroutine igenDS
   
   !********************************************************************

   subroutine leituraGeracaoConectividadesDs(keyword_name, conecElem, mat, nen)
      !> Subrotina respons�vel por ler e gerar conectividades nodais e ladais.
      !> @param keyword_name  O nome da keyword associada.
      !> @param conecElem     Código do elemento
      !> @param mat           Código do material
      !> @param nen           Número de elementos.
    
      use mLeituraEscrita , only: genel1
      integer*4:: n,nel(3),incel(3),inc(3)

      integer*4:: nen
      integer*4:: conecElem(nen,*),mat(*) 
      integer*4:: itemp(27) !B
      character(len=50) keyword_name

      integer*4:: m,ng,i, keyword_line

      ! write(*,'(a)') " em subroutine leituraGeracaoConectividadesDS(keyword_name, conecElem, mat, nen)"
      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.0) return

      100 continue
         read(file_lines(keyword_line),1000) n,m,(itemp(i),i=1,nen),ng
         ! write(*,1000) n,m,(itemp(i),i=1,nen),ng
         keyword_line = keyword_line + 1

         if (n.eq.0) return
         ! call imove(conecElem(1,n),itemp,nen)
         conecElem(1:nen,n)=itemp(1:nen)
         mat(n)=m
         if (ng.ne.0) then
            !Generate data
            read(file_lines(keyword_line),1000) (nel(i),incel(i),inc(i),i=1,3)
            keyword_line = keyword_line + 1
            !write(*,1000) (nel(i),incel(i),inc(i),i=1,3)
            call genel1(conecElem,mat,nen,n,nel,incel,inc)
         endif
      go to 100
      1000 format(16i10,10x,14i10)
   end subroutine

   !********************************************************************

   subroutine leituraMatrizInteiros(keyword_name, matrx, nlin, ncol)
      !> Subrotina responsável por ler uma matriz de inteiros
      !> @param keyword_name  O nome da keyword associada.
      !> @param matrx         Matriz para armazenar os elementos
      !> @param nlin          Número de linhas da matriz
      !> @param ncol          Número de colunas da matriz
      INTEGER, INTENT(IN)  :: nlin, ncol
      INTEGER, INTENT(OUT) :: matrx(nlin, ncol)
      character(len=50) keyword_name

      INTEGER :: i, j, m

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.0) return

      do i=1,nlin
         read (file_lines(keyword_line), *) m, (matrx(i,j),j=1,(ncol))
         keyword_line = keyword_line + 1
      enddo
   end subroutine

   !********************************************************************

   subroutine leituraMatrizInteiros2(keyword_name, matrx, nlin, ncol)
      !> Subrotina responsável por ler uma matriz de inteiros sem indice
      !> @param keyword_name  O nome da keyword associada.
      !> @param matrx         Matriz para armazenar os elementos
      !> @param nlin          Número de linhas da matriz
      !> @param ncol          Número de colunas da matriz
      INTEGER, INTENT(IN)  :: nlin, ncol
      INTEGER, INTENT(OUT) :: matrx(nlin, ncol)
      character(len=50) keyword_name

      INTEGER :: i, j

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.0) return

      do i=1,nlin
         read (file_lines(keyword_line), *) (matrx(i,j),j=1,(ncol))
         keyword_line = keyword_line + 1
      enddo
   end subroutine

   !********************************************************************

   subroutine genelFacesDS(keyword_name, conecElem, nen, nelx, nely, nelz)
      !> Subrotina responsável por ler e gerar elementos de face.
      !> @param keyword_name     A plavra-chave associada.
      !> @param conecElem        Código dos nós dos elementos.
      !> @param nen              N�mero de lementos.
      !> @param nelx             N�mero de elmentos em x.
      !> @param nely             N�mero de elmentos em y.
      !> @param nelz             N�mero de elmentos em z.

      use mGlobaisEscalares
      use mMalha, only: numel

      implicit none
      integer*4:: nen, nelx, nely, nelz
      integer*4:: conecElem(nen,*)
      character(len=50) keyword_name

      integer*4:: ng, n, m, nel, i
      integer*4:: condicao, condicao2, keyword_line
        
      print*, "nelx", nelx
      print*, "nely", nely

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.0) return

      read(file_lines(keyword_line),1000) n,m,(conecElem(i,1),i=1,nen),ng
      keyword_line = keyword_line + 1
      ! write(*,1000) n,m,(conecElem(i,1),i=1,nen),ng

      condicao=0
      condicao2=0

      do nel=2, numel
         if(condicao==0.and.condicao2==0) then
            do i=1, nen
               conecElem(i,nel)=conecElem(i,nel-1)+1
            end do
         else
            if(condicao==1.and.condicao2==0) then
               do i=1, nen
                  if(i<=4) then
                     conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+1
                  else
                     conecElem(i,nel)=conecElem(i,nel-1)+1
                  endif
               end do
            else
               if(condicao==1.and.condicao2==1) then
                  do i=1, nen
                     if(i<=4) then
                        conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+nelx+(nelx*nely)+1
                     else
                        conecElem(i,nel)=conecElem(i,nel-1)+(nelx*(nely+1))+(nely*(nelx+1))+1
                     endif
                  enddo
               end if
            end if
         end if
         if(mod(nel, nelx)==0) then
            condicao=1
         else
            condicao=0
         end if
         if(mod(nel, nelx*nely)==0) then
            condicao2=1
         else
            condicao2=0
         end if
      end do
      1000 format(16i10,10x,14i10)
   end subroutine
    
   !********************************************************************

   subroutine readNodeElementsDS
      use mGlobaisEscalares
      use mMalha,          only: nen
      implicit none
      character(len=50) keyword_name
      ! integer*4 :: ntype, numat

      keyword_name = "ntype"
      call readIntegerKeywordValue(keyword_name, ntype, ntype)
      keyword_name = "numat"
      call readIntegerKeywordValue(keyword_name, numat, numat)
      keyword_name = "nen"
      call readIntegerKeywordValue(keyword_name, nen, nen)
      keyword_name = "nicode"
      call readIntegerKeywordValue(keyword_name, nicode, nicode)
   end subroutine
    
   !********************************************************************

   subroutine readMaterialPropertiesDS(keyword_name)
      !> Efetua a leitura de propriedades de materiais.
      !> @param keyword_name  Keyword especifica das  propriedades de materiais.
      use mGlobaisEscalares
      use mGlobaisArranjos
      use mleituraEscrita, only: iecho

      implicit none
      character(len=50) keyword_name
      integer*4 n, m, i, keyword_line

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.0) return

      do 400 n=1,numat
         if (mod(n,50).eq.1) write(iecho,4000) numat
         read (file_lines(keyword_line),  5000) m,(c(i,m),i=1,3)
         keyword_line = keyword_line + 1
         write(iecho,6000) m,(c(i,m),i=1,3)
      400 continue
      5000  format(i10,5x,5f10.0)
      6000  format(2x,i3,1x,5(1x,1pe11.4))
      4000  format(///,&
               ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
               ' number of material sets . . . . . . . . . . (numat ) = ',i10///,2x,'set',4x,'Kx ',&
               10x,'Ky',10x,'Kz')
   end subroutine readMaterialPropertiesDS
   
   !********************************************************************

   subroutine readConstantBodyForcesDS(keyword_name)
      !> Faz a leitura de constant body forces.
      !> @param keyword_name  Keyword especifica para constant body forces.
    
      use mGlobaisArranjos, only: grav
      use mleituraEscrita, only: iecho

      implicit none
      character(len=50) keyword_name
      integer*4 i, keyword_line

      keyword_line = findKeyword(keyword_name)
      if (keyword_line.eq.0) return

      read  (file_lines(keyword_line),  7000) (grav(i),i=1,3)
      write (iecho,8000) (grav(i),i=1,3)

      7000 format(8 f10.0)
      8000 format(///,&
               ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
               ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
               ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
               ' exemplo 3............................ = ',      1pe15.8,//)
   end subroutine readConstantBodyForcesDS
    
   !********************************************************************

end module mInputReader


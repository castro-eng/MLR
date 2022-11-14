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
!
module mLeituraEscrita
      
   use mGlobaisEscalares, only: saltoPressao

   implicit none

   integer*4:: iin, iecho, icoords, iconects, iconectsL
   integer*4:: ignuplotPotencial, ignuplotFluxo
   integer*4 :: iparaviewPotencial, iparaviewFluxo, iparaviewPerm, iparaviewMat
   integer*4 :: iparaview
   integer*4:: tipo_arq_saida, passosPorImpressao
   logical :: plotProfundidade

   integer*4:: iconectIntersecao, iconectIntersecaoOut, iaberturasBorda, iaberturasBordaOut
      

   ! funcoes e subrotinas
   public :: echo, leituraValoresCondContorno, leituraGeracaoCoordenadas
   public :: printf, printd, printp, prntel
   public :: printResultado, prtgnup, prtvB
   public :: escreverArqParaview, escreverPontosNodais, escreverConectividades
   public :: escreverTiposElementos, escreverEscalaresNodais
   public :: leituraGeracaoConectividades
   public :: genel1, escreverArquivosSaida_Fluxo, escreverArquivosSaida_Potencial, escreverValoresReferencia

   contains

   !=================================================================================

   subroutine echo
      implicit none

      !.... program to echo input data
      character*4 ia(50)
      integer*4:: iech, i

      ! cabeçalho
      write(iecho,500)

      read(iin,1000) iech
      if (iech.eq.0) return

      write(iecho,2000) iech
      backspace iin

      do 100 i=1,100000
         read(iin,3000,end=200) ia
         if (mod(i,50).eq.1) write(iecho,4000)
         write(iecho,5000) ia
      100 continue

      200 continue
      rewind iin
      read(iin,1000) iech

      return

      500  format('programa de elementos finitos em fortran 90 baseado em:',// &
               'The Finite Element Method, Hughes, T. J. R., (2003)'//)
      1000 format(16i10)
      2000 format('1',' i n p u t   d a t a   f i l e               ',  //5x,&
               ' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x,&
               '    eq. 0, no echo of input data                        ',   /5x,&
               '    eq. 1, echo input data                              ',   ///)
      3000 format(20a4)
      4000 format(' ',8('123456789*'),//)
      5000 format(' ',20a4)
   end subroutine echo

   !=================================================================================

   subroutine leituraGeracaoConectividades(conecElem,mat,nen, iin)   

      use mGlobaisEscalares 

      implicit none

      !.... program to read and generate element node and material numbers    

      !         conecElem(nen,numel) = element node numbers                         
      !         mat(numel)     = element material numbers                     
      !         nen            = number of element nodes (le.27)              
      !         n              = element number                               
      !         ng             = generation parameter                         
      !         nel(i)         = number of elements in direction i            
      !         incel(i)       = element number increment for direction i     
      !         inc(i)         = node number increment for direction i        

      integer*4:: nen, iin
      integer*4:: conecElem(nen,*),mat(*)

      integer*4:: m,ng, i, itemp(27)
      integer*4:: n,nel(3),incel(3),inc(3)        
                       
      100 continue                                                          
         read(iin,1000) n,m,(itemp(i),i=1,nen),ng     
                
         if (n.eq.0) return                                                
         conecElem(1:nen,n)=itemp(1:nen)                                    
         mat(n)=m                                                          
         if (ng.ne.0) then
            !....... generate data                                                     
            read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)                     
            call genel1(conecElem,mat,nen,n,nel,incel,inc)                     
         endif
      go to 100                                                         
                                                                      
      1000 format(16i10,10x,14i10)                                             
                                                                      
   end subroutine                                                               

   !******************************************************************************

   subroutine genel1(conecElem,mat,nen, n,nel,incel,inc) 
      !.... program to generate element node and material numbers             

      integer*4 ::  nen, conecElem(nen,*),mat(*)                                       
      integer*4:: n,nel(3),incel(3),inc(3)                          

      integer*4:: i,j,k,ii,jj,kk,ie,le,je,ke
                                                                      
      !.... set defaults                                                      
                                                                      
      call geneld                                                 
                                                                      
      !.... generation algorithm                                              
                                                                      
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
                                                                     
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       

      do 300 k=1,kk
         do 200 j=1,jj
            do 100 i=1,ii

               if (i.ne.ii) then
                  le = ie
                  ie = le + incel(1)
                  call geneli(conecElem(1,ie),conecElem(1,le),inc(1),nen)
                  mat(ie) = mat(le)
               endif
            100 continue

            if (j.ne.jj) then
               le = je
               je = le + incel(2)
               call geneli(conecElem(1,je),conecElem(1,le),inc(2),nen)
               mat(je) = mat(le)
               ie = je
            endif
         200 continue                                                          

         if (k.ne.kk) then
            le = ke
            ke = le + incel(3)
            call geneli(conecElem(1,ke),conecElem(1,le),inc(3),nen)
            mat(ke) = mat(le)
            ie = ke
            je=ke
         endif
      300 continue
                                                                     
      return
      contains

      !******************************************************************************
      subroutine geneld
         !.... program to set defaults for element node
         !        and material number generation

         if (nel(1).eq.0) nel(1) = 1
         if (nel(2).eq.0) nel(2) = 1
         if (nel(3).eq.0) nel(3) = 1
                                                                      
         if (incel(1).eq.0) incel(1) = 1
         if (incel(2).eq.0) incel(2) = nel(1)
         if (incel(3).eq.0) incel(3) = nel(1)*nel(2)
                                                                      
         if (inc(1).eq.0) inc(1) = 1
         if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)
         if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)

         return
      end subroutine

      !******************************************************************************

      subroutine geneli(conecElem2,conecElem1,inc,nen)
         !.... program to increment element node numbers

         integer*4 ::  conecElem1(*),conecElem2(*)
         integer*4:: inc, nen

         integer*4:: i

         do 100 i=1,nen
            if (conecElem1(i).eq.0) then
               conecElem2(i) = 0
            else
               conecElem2(i) = conecElem1(i) + inc
            endif
         100 continue

         return                                                            
      end  subroutine
   end subroutine genel1                                                              

   !**** new **********************************************************************

   subroutine leituraGeracaoCoordenadas(x, nsd, numnp, iin, icoords, iprtin)
      !.... program to read, generate and write coordinate data
      use mMalha, only: genfl

      implicit none

      integer*4, intent(in)   :: nsd, numnp, iin, icoords, iprtin
      real*8, intent(inout) ::  x(nsd,*)

      integer*4:: i, n
     
      call genfl(x,nsd,iin)

      if (iprtin.eq.1) return

      write(icoords,*) "# Coordenadas "
      do n=1,numnp
         write(icoords,2000) n,(x(i,n),i=1,nsd)   
      end do

      return

      ! 1000 format('1',' n o d a l   c o o r d i n a t e   d a t a '///5x,&
      !          ' node no.',3(13x,' x',i1,' ',:)//)
      2000 format(6x,i12,10x,3(1pe15.8,2x))
   end subroutine

   !**** new **********************************************************************

   subroutine leituraCodigosCondContorno(id, ndof, numnp, neq, iin, iecho,iprtin)
      !.... program to read, generate and write boundary condition data
      !        and establish equation numbers

      use mMalha, only: igen

      integer*4, intent(in) :: ndof, numnp, iin, iecho, iprtin
      integer*4, intent(inout) :: neq
      integer*4, intent(inout) :: id(ndof,numnp)

      integer*4:: nn, n, i
      logical pflag

      id(1:ndof,1:numnp) = 0
      call igen(id,ndof, iin)

      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
            pflag = .false.

            do 100 i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
               ! print*, i, 'id =',  id(i,n)
            100 continue

            if (pflag) then      
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
               write(iecho,2000) n,(id(i,n),i=1,ndof)
            endif
         200 continue
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
               c o  d e s'/// &
               5x,' node no.',3x,6(6x,'dof',i1:)//)
      2000 format(6x,i10,5x,6(5x,i10))

   end subroutine

   !**** new **********************************************************************

   subroutine leituraCodigosCondContornoBidu(id, ndof, numnp, neq, iin, iecho,iprtin)
      !.... program to read, generate and write boundary condition data
      !        and establish equation numbers

      use mMalha, only: igen

      integer*4, intent(in) :: ndof, numnp, iin, iecho, iprtin
      integer*4, intent(inout) :: neq
      integer*4, intent(inout) :: id(ndof,numnp)

      integer*4:: nn, n, i
      logical pflag

      write(*,*) "subroutine leituraCodigosCondContornoB(id, ndof, numnp, ..."
      id(1:ndof,1:numnp) = 0
      call igen(id,ndof, iin)

      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
            pflag = .false.

            do 100 i=1,ndof
               if (id(i,n).ne.0) pflag = .true.
               ! print*, i, 'id =',  id(i,n)
            100 continue

            if (pflag) then      
               nn = nn + 1
               if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
               write(iecho,2000) n,(id(i,n),i=1,ndof)
            endif
         200 continue
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

      ! do i =1, numnp    !BD
      !    print*,'i=', i, 'id =',  id(1,i)
      ! end do
      i = 3
      i = i+1; id(1,i) = 4 + 1
      i = i+1; id(1,i) = 6 + 1
      i = i+1; id(1,i) = 11 + 1
      i = i+1; id(1,i) = 1 + 1
      i = i+1; id(1,i) = 5 + 1
      i = i+1; id(1,i) = 9 + 1
      i = i+1; id(1,i) = 10 + 1
      i = i+1; id(1,i) = 2 + 1
      i = i+1; id(1,i) = 3 + 1
      i = i+1; id(1,i) = 8 + 1
      i = i+1; id(1,i) = 7 + 1
      i = i+1; id(1,i) = 0 + 1

      ! do i =1, numnp    !BD
      !    print*,'i=', i, 'id =',  id(1,i)
      ! end do

      return

      1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n & 
               c o  d e s'/// &
               5x,' node no.',3x,6(6x,'dof',i1:)//)
      2000 format(6x,i10,5x,6(5x,i10))

   end subroutine

   !**** new **********************************************************************
   subroutine leituraValoresCondContorno(f,ndof,numnp,j,nlvect,iprtin)
      !.... program to read, generate and write nodal input data

      !        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
      !                             = nodal body forces(j=1)

      use mMalha, only : genfl
      implicit none

      integer*4:: ndof, numnp, j, nlvect, iprtin
      real*8 f(ndof,numnp,nlvect)

      logical lzero
      integer*4 nlv
      character(len=30) :: rotulo

      ! call clear(f,nlvect*numnp*ndof)
      f=0.0      

      do 100 nlv=1,nlvect
         call genfl(f(1,1,nlv),ndof,iin)
         call ztest(f(1,1,nlv),ndof*numnp,lzero)

         if (iprtin.eq.0) then

            if (lzero) then
               if (j.eq.0) write(iecho,1000) nlv
               if (j.eq.1) write(iecho,2000)
            else
               if (j.eq.0) call printf(f,ndof,numnp,nlv)
               if (j.eq.1) then
                  rotulo=" n o d a l  b o d y  f o r c e s  "
                  call printd (rotulo, f,ndof,numnp,iecho)
               end if
            endif
         endif
      100 continue

      return
      1000 format('1'//,' there are no nonzero prescribed forces and ',&
               'kinematic boundary conditions for load vector number ',i10)
      2000 format('1'//,' there are no nonzero nodal body forces')
   end subroutine

   !**** new **********************************************************************

   subroutine printf(f,ndof,numnp,nlv)
      !.... program to print prescribed force and boundary condition data

      implicit none

      integer*4 :: ndof, numnp, nlv
      real*8 :: f(ndof,numnp,*)

      logical lzero
      integer*4:: nn, n, i

      nn = 0

      do 100 n=1,numnp
         call ztest(f(1,n,nlv),ndof,lzero)
         if (.not.lzero) then
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) nlv,(i,i=1,ndof)
            write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
         endif
      100 continue

      return

      1000 format('1',&
               ' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',&
               '  b o u n d a r y   c o n d i t i o n s'//5x,&
               ' load vector number = ',i10///5x,&
               ' node no.',6(13x,'dof',i1,:)/)
      2000 format(6x,i10,10x,6(1pe15.8,2x))
   end subroutine

   !**** new **********************************************************************

   subroutine printd(name,dva,ndof,numnp,icode)
      !.... program to print kinematic data

      implicit none

      integer*4:: ndof,numnp, icode
      character (LEN=*) ::  name
      real*8 dva(ndof,*)

      logical lzero
      integer*4 nn, n, i

      nn = 0

      do 100 n=1,numnp
         call ztest(dva(1,n),ndof,lzero)
         if (.not.lzero) then
            nn = nn + 1
            if (mod(nn,50).eq.1) &
            write(icode,1000) name,(i,i=1,ndof)
            write(icode,2000) n,(dva(i,n),i=1,ndof)
         endif
      100 continue

      return

      1000 format('1',11a4//1x,'node',6(11x,'dof',i1)/)
      2000 format(1x,i10,2x,6(1pe13.6,2x))
   end subroutine

   !**** new **********************************************************************

   subroutine printp(a,idiag,neq,nsq,*)
      !.... program to print array d after Crout factorization 
      !        a = u(transpose) * d * u

      use mGlobaisEscalares

      implicit none

      integer*4:: neq, nsq
      real*8 :: a(*)
      integer*4:: idiag(*)

      integer*4:: n, i

      do 100 n=1,neq
         if (mod(n,50).eq.1) write(iecho,1000) nsq
         write(iecho,1000)
         i = idiag(n)
         write(iecho,2000) n,a(i)
      100 continue

      return 1

      1000 format('1',' array d of factorization',/&
               ' a = u(transpose) * d * u ',                                //5x,&
               ' time sequence number   . . . . . . . . . . . . (nsq) = ',i10//5x)
      2000 format(1x,i10,4x,1pe20.8)
   end subroutine

   !**** new **********************************************************************

   subroutine prntel(mat,conectElem,nen,numel,tipo)
      !.... program to print data for element with "nen" nodes
      !        note: presently the label formats are limited to
      !              elements with one to nine nodes
      implicit none

      integer*4:: nen, numel
      integer*4:: mat(*),conectElem(nen,*)
      integer*4:: tipo

      integer*4 n, i

      if(tipo==1) then
         write(iconects,*) "# Conectividades nodais"
         do n=1,numel
            write(iconects,2000) n,mat(n),(conectElem(i,n),i=1,nen)
         end do
      end if

      if(tipo==2) then
         write(iconectsL,*) "# Conectividades ladais"
         do  n=1,numel
            write(iconectsL,3000) n,mat(n),(conectElem(i,n),i=1,nen)
         end do
      end if

      return

      1000 format(///10x,&
               ' e l e m e n t   d a t a',//1x,&
               ' element   material ',9('node ',i1,1x),/1x,&
               '  number    number',//)
      2000 format(1x,i10,9(2x,i10))
      3000 format(1x,i10,7(2x,i10))
   end subroutine

   !**** new **********************************************************************

   subroutine printResultado(dva, ndof, numnp, inicio, fim, icode)
      !.... program to print kinematic data

      implicit none

      integer*4:: ndof, numnp, inicio, fim, icode
      real*8  :: dva(ndof,numnp)

      integer*4:: n, i

      write(icode,*) "# Solucao"
      do 100 n=inicio,fim
         write(icode,2000) n,(dva(i,n),i=1,ndof)
         ! write(*,*) n,(dva(i,n),i=1,ndof)
      100 continue

      return
      2000 format(1x,i10,2x,6(1pe13.6,2x))
   end subroutine

   !**** new **********************************************************************

   subroutine prtgnup(name,x,dva,nsd,ndof,numnp,icode)
      !.... program to print kinematic data

      implicit none

      integer*4:: nsd, ndof, numnp, icode
      character(len=*) :: name
      real*8 :: x(nsd,*),dva(ndof,*)

      integer*4:: n
      real*8  :: tmp, t1, t2

      call timing(t1)
      write(icode,'(a)') name
      tmp=0.0
      do 100 n=1,numnp
         write(icode,2000) x(1:nsd,n), dva(1:ndof,n)
         if(x(nsd,n+1)>x(nsd,n).and.n<numnp) write(icode,*) ""
      100 continue
      ! do 100 n=1,numnp
      !    write(icode,2000) (x(j,n),j=1,nsd), (dva(i,n),i=1,ndof)
      ! 100 continue

      call timing(t2)
      write(*,*) 'tempo de escrita em arquivo = ', t2 - t1
      return

      2000 format(6(1pe13.6,2x))
   end subroutine

   !=================================================================================

   subroutine escreverArquivosSaida_Potencial(estrutSistEqP, passoTempo, tempo)
      use mEstruturasDadosSistEq, only: estruturasArmazenamentoSistemaEq
      use mMalha, only: nsd, numelFratura, numnpFratura, numnp, nen
      use mMalha, only: conecNodaisElem, conecNodaisFratura,nosFratura_impressao_global

      type(estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqP
      integer, intent(in) :: passoTempo
      real*8, intent(in) :: tempo
         
      real*8 :: pressao(2,numnpFratura), pressaoF(numnpFratura), saltoP(numnpFratura)
      integer :: cont, no, i

      character(21) :: labelTempo

      if (tipo_arq_saida == 1) then
         call escreverArquivoVTU_Pressao(estrutSistEqP%u, estrutSistEqP%ndof, numnp, passoTempo)
         if (numelFratura > 0)  then
            do i=1,numnpFratura
               no = nosFratura_impressao_global(i)
               pressaoF(i) = estrutSistEqP%u(1, no)
            enddo              
            call escreverArquivoVTU_Fraturas_Pressao(pressao, saltoP, pressaoF, estrutSistEqP%ndof, numnpFratura, passoTempo)
         endif   

      else if (tipo_arq_saida == 2) then
         call escreverArqUnicoVTU(estrutSistEqP%u(1,:), passoTempo)
         ! call escreverArqUnicoVTU(estrutSistEqP%uIterAnt(1,:), passoTempo)

      else if (tipo_arq_saida == 0) then
         call gerarLabel(labelTempo,tempo)
         if(passoTempo==0) then
            open(unit=iparaviewPotencial ,file= './out/resultadoPotencial.vtk')
            if((saltoPressao.eqv..true.).and.(numelFratura>0)) then
               call escreverArqParaviewVector(estrutSistEqP%u, estrutSistEqP%ndof, numnp, nen, &
                        conecNodaisElem, 2, trim(labelTempo), len(trim(labelTempo)), iparaviewPotencial)                     
            else
               call escreverArqParaview(estrutSistEqP%u, estrutSistEqP%ndof, numnp, nen, &
                        conecNodaisElem, 2, trim(labelTempo), len(trim(labelTempo)), iparaviewPotencial)
            endif

         else
            call escreverArqParaviewIntermed(estrutSistEqP%u, estrutSistEqP%ndof, numnp, 'potencial', &
                     trim(labelTempo), len(trim(labelTempo)), iparaviewPotencial)
         end if
      end if
   end subroutine escreverArquivosSaida_Potencial

   !=================================================================================
  
   subroutine escreverArquivosSaida_Fluxo(estrutSistEqF, passoTempo, tempo)
      use mEstruturasDadosSistEq, only: estruturasArmazenamentoSistemaEq
      use mMalha, only: nsd, numnp, nen, nosFratura_impressao_global
      use mMalha, only: conecNodaisElem, conecNodaisFratura, numelFratura, numnpFratura

      type(estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqF
      integer, intent(in) :: passoTempo
      real*8, intent(in) :: tempo

      character(21) :: labelTempo
         
      real*8 :: velocidadeFratura(nsd,numnpFratura)
      integer :: cont, no, i
         
      if (tipo_arq_saida == 1) then
         call escreverArquivoVTU_Velocidade(estrutSistEqF%u, nsd, numnp, passoTempo)
         if (numelFratura > 0)  then
            do i=1,numnpFratura
               no = nosFratura_impressao_global(i)
               velocidadeFratura(:,i) = estrutSistEqF%u(:, no)
            enddo
            call escreverArquivoVTU_Fraturas_Velocidade(velocidadeFratura, nsd, numnpFratura, passoTempo)
         endif   

      elseif (tipo_arq_saida == 0) then
         call gerarLabel(labelTempo,tempo) 
         if(passoTempo==1) then
            open(unit=iparaviewFluxo     ,file= './out/resultadoFluxo.vtk')
            call escreverArqParaviewVector(estrutSistEqF%u, estrutSistEqF%ndof, numnp, nen, &
                     conecNodaisElem, 2, trim(labelTempo), len(trim(labelTempo)), iparaviewFluxo)

         else                     
            call escreverArqParaviewIntermed(estrutSistEqF%u, estrutSistEqF%ndof, numnp, 'velocidade',   &
                     trim(labelTempo), len(trim(labelTempo)), iparaviewFluxo)
         endif
      end if

   end subroutine escreverArquivosSaida_Fluxo

   !=================================================================================

   subroutine escreverArqParaviewMateriais(campo, nen, conectElem, rotulo, tamRot, iparaview)
      use mMalha, only: x, nsd, numel, numnp

      implicit none
      integer*4, intent(in) :: iparaview
      integer, intent(in) :: campo(numel)
      integer*4:: nen
      integer*4:: conectElem(nen,numel)
      integer*4:: tamRot, j

      character(len=tamRot) :: rotulo
    
      write(iparaview,'(a)')'# vtk DataFile Version 3.0'
      write(iparaview,'(a)')'vtk output'
      write(iparaview,'(a)')'ASCII'
      write(iparaview,'(a)')'DATASET UNSTRUCTURED_GRID'
      write(iparaview,'(a,i10,a)')'POINTS', numnp,' float '

      call escreverPontosNodais  (x, numnp, nsd, iparaview)

      write(iparaview,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
      call escreverConectividades(conectElem, numel, nen, nsd, iparaview)

      write(iparaview,'(a,i10)')'CELL_TYPES ', numel
      call escreverTiposElementos(numel, nen, iparaview)

      write(iparaview,'(a,i10)')'CELL_DATA ', numel

      write(iparaview,'(3a)')'SCALARS ', trim(rotulo), ' float '
      write(iparaview,'(a)')'LOOKUP_TABLE default'

      ! call escreverEscalaresNodais(campo, dim1, dim2,rotulo,tamRot,iparaview)
      do j=1, numel
         write(iparaview,*) campo(j)
      end do
     
   end subroutine escreverArqParaviewMateriais

   !=================================================================================

   subroutine escreverArqParaview(campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot, iparaview)
      use mMalha, only: x, nsd, numel

      implicit none
      integer*4, intent(in) :: dim1, dim2, iparaview
      double precision, intent(in) :: campo(dim1, dim2)
      integer*4:: nen
      integer*4:: conectElem(nen,numel)
      integer*4:: tipo  !tipo=1 para elemento, e tipo=2 para no
      integer*4:: tamRot

      character(len=tamRot) :: rotulo
    
      write(iparaview,'(a)')'# vtk DataFile Version 3.0'
      write(iparaview,'(a)')'vtk output'
      write(iparaview,'(a)')'ASCII'
      write(iparaview,'(a)')'DATASET UNSTRUCTURED_GRID'
      write(iparaview,'(a,i10,a)')'POINTS', dim2,' float '

      call escreverPontosNodais  (x, dim2, nsd, iparaview)

      write(iparaview,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
      call escreverConectividades(conectElem, numel, nen, nsd, iparaview)

      write(iparaview,'(a,i10)')'CELL_TYPES ', numel
      call escreverTiposElementos(numel, nen, iparaview)
      ! call escreverTiposElementos_old(numel, nsd, iparaview)

      if(tipo==1) write(iparaview,'(a,i10)')'CELL_DATA ', numel

      if(tipo==2) write(iparaview,'(a,i10)')'POINT_DATA',  dim2

      write(iparaview,'(3a)')'SCALARS ', trim(rotulo), ' float '
      write(iparaview,'(a)')'LOOKUP_TABLE default'

      call escreverEscalaresNodais(campo, dim1, dim2,rotulo,tamRot,iparaview)

   end subroutine escreverArqParaview

   !=================================================================================

   subroutine escreverArqParaviewVector(campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot, iparaview)
      use mMalha, only: x, nsd, numel

      implicit none
      integer*4, intent(in) :: dim1, dim2, iparaview
      double precision, intent(in) :: campo(dim1, dim2)
      integer*4:: nen
      integer*4:: conectElem(nen,numel)
      integer*4:: tipo  !tipo=1 para elemento, e tipo=2 para no
      integer*4:: tamRot

      character(len=tamRot) :: rotulo
    
      write(iparaview,'(a)')'# vtk DataFile Version 3.0'
      write(iparaview,'(a)')'vtk output'
      write(iparaview,'(a)')'ASCII'
      write(iparaview,'(a)')'DATASET UNSTRUCTURED_GRID'
      write(iparaview,'(a,i10,a)')'POINTS', dim2,' float '

      call escreverPontosNodais  (x, dim2, nsd, iparaview)

      write(iparaview,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
      call escreverConectividades(conectElem, numel, nen, nsd, iparaview)

      write(iparaview,'(a,i10)')'CELL_TYPES ', numel
      call escreverTiposElementos(numel, nen, iparaview)
      ! call escreverTiposElementos_old(numel, nsd, iparaview)

      if(tipo==1) write(iparaview,'(a,i10)')'CELL_DATA ', numel

      if(tipo==2) write(iparaview,'(a,i10)')'POINT_DATA',  dim2

      write(iparaview,'(3a)')'VECTORS ', trim(rotulo), ' float '

      if(tipo==1) call escreverVetoresNodais(campo, dim1, numel,rotulo,tamRot,iparaview)
      if(tipo==2) call escreverVetoresNodais(campo, dim1, dim2, rotulo,tamRot,iparaview)

   end subroutine escreverArqParaviewVector

   !=================================================================================

   subroutine escreverPontosNodais  (coords, numnp, nsd, iparaview)
      implicit none
      integer*4, intent(in) :: numnp, nsd, iparaview
      real*8,  intent(in) :: coords(nsd,numnp)

      real*8  :: coordZ = 0.0 
      integer*4:: d, i

      if(nsd==2) then
         do i=1,numnp
            write(iparaview,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd), coordZ 
         end do
      end if

      if(nsd==3) then
         do i=1,numnp
            write(iparaview,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd)
         end do
      end if
   end subroutine escreverPontosNodais

   !=================================================================================
   subroutine escreverConectividades(conectElem, numel, nen, nsd, iparaview)
      implicit none
      integer*4, intent(in)  :: numel, nen, nsd, iparaview
      integer*4, intent(in)  :: conectElem(nen,numel)

      integer*4 n, i

      if(nsd==2) then
         do  n=1,numel
            write(iparaview,'(i10,9(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
         end do
      end if

      if(nsd==3) then
         do  n=1,numel
            write(iparaview,'(i10,18(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
         end do
      end if

   end subroutine escreverConectividades

   !=================================================================================

   subroutine escreverTiposElementos(numel, nen, iparaview)
      
      use mMalha, only: nsd
      
      implicit none
      integer*4, intent(in)   :: numel, nen, iparaview

      integer*4:: i
      character*2 :: tipo

      if (nen==3) then !triangulo
         tipo = "5"
      else if (nen==4.and.nsd==2) then !quadrilátero
         tipo = "9"
      else if (nen==4.and.nsd==3) then !tetraedro
         tipo = "10"
      else if(nen==8) then !hexaedro
         tipo = "12"
      else
         write(*,*) "Erro: subrotina 'escreverTiposElementos' nao esta implementada para esse tipo de elemento. Valor de nen: ", nen
         stop 1
      end if

      do i = 1, numel
         write(iparaview,'(a)') trim(adjustl(tipo))
      end do

   end subroutine escreverTiposElementos

   !=================================================================================

   subroutine escreverTiposElementos_old(numel, nsd, iparaview)
      implicit none
      integer*4, intent(in)   :: numel, nsd, iparaview

      integer*4:: i

      if(nsd==2) then
         do  i =1,numel
            write(iparaview,'(a)') '9'!trim(adjustl(tipo))
         end do
      end if 

      if(nsd==3) then
         do  i =1,numel
            write(iparaview,'(a)') '12'!trim(adjustl(tipo))
         end do
      end if 
   end subroutine escreverTiposElementos_old

   !=================================================================================

   subroutine escreverEscalaresNodais(v, tam1, tam2, rotulo, tamRot, iparaview)
      implicit none
      integer*4, intent(in)  :: tam1,tam2,iparaview
      real*8, intent(in)   :: v(tam1,tam2)
      integer*4:: tamRot
      character(len=tamRot) :: rotulo

      character(len=tamRot+5) ::  rotuloN
      integer*4:: i,j
      character(len=5):: eixo
      real*8 :: limite

      limite=1.e-20
      do i=1,tam1
         if(i>1) then
            write(eixo,'(i0)') i
            if(rotulo.ne.'potencial') then
               rotuloN=trim(rotulo)//'Dir'//trim(eixo)
               write(iparaview,'(3a)')'SCALARS ', trim(rotuloN), ' float '
               write(iparaview,'(a)')'LOOKUP_TABLE default'
            endif
         endif

         do j=1, tam2
            write(iparaview,*) v(i,j)
         end do
      end do

   end subroutine escreverEscalaresNodais

   !=================================================================================

   subroutine escreverVetoresNodais(v, tam1, tam2, rotulo, tamRot, iparaview)
      implicit none
      integer*4, intent(in)  :: tam1,tam2,iparaview
      real*8, intent(in)   :: v(tam1,tam2)
      integer*4:: tamRot
      character(len=tamRot) :: rotulo

      character(len=tamRot+5) ::  rotuloN
      integer*4:: i,j
      character(len=5):: eixo
      real*8 :: limite
      
      real*8   :: vPrint(tam1,tam2)

      vPrint=v
      
      limite=1.e-20

      write(eixo,'(i0)') i
      if(rotulo.ne.'potencial') then
         rotuloN=trim(rotulo)//'Dir'//trim(eixo)
      endif

      do j=1, tam2
         if(tam1==2) then
            if(abs(vPrint(1,j))<limite) vPrint(1,j)=0.d0
            if(abs(vPrint(2,j))<limite) vPrint(2,j)=0.d0
              
            write(iparaview,*) vPrint(1:tam1,j), '   0.0'
         else
        
            if(abs(vPrint(1,j))<limite) vPrint(1,j)=0.d0
            if(abs(vPrint(2,j))<limite) vPrint(2,j)=0.d0
            if(abs(vPrint(3,j))<limite) vPrint(3,j)=0.d0

            write(iparaview,*) vPrint(1:tam1,j)
         endif
      end do

   end subroutine escreverVetoresNodais

   !=================================================================================

   subroutine escreverEscalaresPorElemento(v, tam1, tam2, rotulo, tamRot, iparaview)
      implicit none
      integer*4, intent(in)  :: tam1,tam2,iparaview
      real*8, intent(in)   :: v(tam1,tam2)
      integer*4:: tamRot
      character(len=tamRot) :: rotulo

      character(len=tamRot+5) ::  rotuloN
      integer*4:: i,j
      character(len=5):: eixo
      real*8 :: limite

      limite=1.e-20
      do i=1,tam1
         if(i>1) then
            write(eixo,'(i0)') i
            if(rotulo.ne.'potencial') then
               rotuloN=trim(rotulo)//'Dir'//trim(eixo)
               write(iparaview,'(3a)')'SCALARS ', trim(rotuloN), ' float '
               write(iparaview,'(a)')'LOOKUP_TABLE default'
            endif
         endif

         do j=1, tam2
            write(iparaview,*) v(i,j)
         end do
      end do

   end subroutine escreverEscalaresPorElemento

   !=================================================================================

   subroutine escreverArqParaviewIntermed(campo, dim1, dim2, nomeCampo, rotulo, tamRot,iparaview)
      use mMalha, only: x, nsd, numel, numnp, numelFratura

      implicit none
      integer*4, intent(in) :: dim1, dim2, iparaview
      double precision, intent(in) :: campo(dim1, dim2)

      integer*4:: tamRot
      character(len=tamRot) :: rotulo
      character(len=*) :: nomeCampo

      if(nomeCampo.eq.'potencial')then
         if((saltoPressao.eqv..true.) .and. (numelFratura>0)) then
            write(iparaview,'(3a)')'VECTORS ', trim(rotulo), ' float '
            call escreverVetoresNodais(campo, dim1, dim2,rotulo,tamRot,iparaview)
         else
            write(iparaview,'(3a)')'SCALARS ', trim(rotulo), ' float '
            write(iparaview,'(a)')'LOOKUP_TABLE default'
            call escreverEscalaresNodais(campo, dim1, dim2, rotulo, tamRot, iparaview)
         endif
      else    
         write(iparaview,'(3a)')'VECTORS ', trim(rotulo), ' float '
         call escreverVetoresNodais(campo, dim1, dim2,rotulo,tamRot,iparaview)
      endif

   end subroutine escreverArqParaviewIntermed

   !=================================================================================
   ! *** BEGIN -- ROTINA PARA ESCREVER MALHA PARA O PARAVIEW (AUTOR: KADU/PROMEC 06/2011)

   subroutine escreverArqParaview_geo()
      use mMalha, only: x, conecNodaisElem, nsd, numel, numnp

      implicit none
      integer*4:: i , iparaview_geo  

      iparaview_geo  = 31
      open(unit=iparaview_geo , file= 'modelo.geo')

      write(iparaview_geo,'(a)')'Title1'
      write(iparaview_geo,'(a)')'Title2'
      write(iparaview_geo,'(a)')'node id given'
      write(iparaview_geo,'(a)')'element id given'
      write(iparaview_geo,'(a)')'coordinates'
      write(iparaview_geo,'(i8)')  numnp

      do i = 1, numnp
         WRITE (iparaview_geo,'(I8,3E12.5)') I,x(1,i),x(2,i),0.d0
      enddo

      WRITE (iparaview_geo,'(A,/,A,/,A,/,I8)') 'part 1', 'malha', 'quad4', numel

      WRITE (iparaview_geo,'(5I8)')  ( I, conecNodaisElem(1,i), conecNodaisElem(2,i), &
               conecNodaisElem(3,i), conecNodaisElem(4,i),i=1, numel ) 

   end subroutine escreverArqParaview_geo

   !=================================================================================


   subroutine escreverArqParaview_res(campo, ndim , npasso, var)
      ! *** BEGIN -- ROTINA PARA ESCREVER ESCALAR POR ELEMENTO PARA O PARAVIEW (AUTOR: KADU/PROMEC 06/2011)

      implicit none
      integer*4, intent(in) :: ndim, npasso
      double precision, intent(in) :: campo(ndim)
      character*1, intent(in) :: var

      integer*4:: i
      character*30 :: filename
      integer*4, parameter :: iparaview_res=41

      write (filename(1:3),'(i3.3)') npasso
      if (var=='p') write (filename,'(a)')   'potencial'//'.'//filename(1:3)
      if (var=='s') write (filename,'(a)') 'saturation'//'.'//filename(1:3)

      open  (unit=iparaview_res,file=filename,form='formatted')

      write (iparaview_res,'(a,i5,1x,a)') 'Ensight Escalar passo ',npasso
      write (iparaview_res,'(a/a)') 'part 1','quad4' 
      write (iparaview_res,'(1p,6e12.5)') (campo(i),i=1,ndim)

      close  (unit=iparaview_res)

   end subroutine escreverArqParaview_res

   !=======================================================================

   subroutine prt(nen,nsd,numel,conectElem,t0,u,iunit)
      ! imprime campos escalares para o gnuplot ou para o matlab

      use mMalha, only: x
      use mMalha, only: local
     
      implicit none

      integer*4                  :: nen,nsd,numel
      integer*4                  :: conectElem(nen,numel)
      real(8), dimension(*)     :: u
      real(8)                   :: t0
    
      integer*4:: nel
      real(8) :: xx,yy,uu
    
      integer*4:: iunit

      real*8 :: xg, yg
      real*8 :: xl(nsd,nen)
    
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
      write(iunit,*)
    
      do nel=1,numel
         call local(conectElem(1,nel),x,xl,nen,nsd,nsd)
         xg = sum(xl(1,1:nen))/nen
         yg = sum(xl(2,1:nen))/nen
         xx=xg !xc(1,nel)
         yy=yg !xc(2,nel)

         uu=u(nel)

         write(iunit,"(5(f25.15,2x))") xx,yy,uu 
      end do
    
      write(iunit,*)    
   end subroutine
    
   !=======================================================================

   subroutine prtvB(nen,nsd,numel,conecNodaisElem,t0,velocLadal,ndofV, numLados, conecLadaisElem, numLadosElem, iunit)
      ! imprime campos vetoriais para o gnuplot ou para o matlab

      use mMalha, only: x
      use mMalha, only: local

      implicit none

      integer*4                  :: nen,numel,nsd, ndofV ,numLados,numLadosElem
      real(8), dimension(ndofV,numLados) :: velocLadal
      real(8)                   :: t0
      integer*4                  :: conecNodaisElem(nen,numel),conecLadaisElem(numLadosElem,numel)

      integer*4:: nel
      real*8 :: xl(nsd,nen)
      real(8) :: xg,yg, vcx,vcy
      integer*4:: lado1, lado2, lado3, lado4

      integer*4:: iunit

      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
      write(iunit,*)

      do nel=1,numel
        call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
        xg = sum(xl(1,1:nen))/nen
        yg = sum(xl(2,1:nen))/nen

        lado1 = conecLadaisElem(1,nel);  lado2 = conecLadaisElem(2,nel);
        lado3 = conecLadaisElem(3,nel);  lado4 = conecLadaisElem(4,nel);
        vcx = (velocLadal(1,lado2)+velocLadal(1,lado4))/2.0
        vcy = (velocLadal(1,lado1)+velocLadal(1,lado3))/2.0
        write(iunit,"(5(f25.15,2x))") xg,yg,vcx,vcy ! velocLadal(1,nel),velocLadal(2,nel)
      end do

      write(iunit,*)
   end subroutine

   !=======================================================================

   subroutine prtv(nen,nsd,numel,ndofV,numLados,conecElem,t0,u,iunit)
      ! imprime campos vetoriais para o gnuplot ou para o matlab

      use mMalha, only: x
      use mMalha, only: local

      implicit none

      integer*4                  :: nen,numel,nsd, ndofV, numLados
      real(8), dimension(ndofV,numLados) :: u
      real(8)                   :: t0
      integer*4                  :: conecElem(nen,numel)

      integer*4:: nel
      real*8 :: xl(nsd,nen)
    
      real(8) :: xx,yy
    
      integer*4:: iunit
    
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
      write(iunit,*)
    
      do nel=1,numel    
         call local(conecElem(1,nel),x,xl,nen,nsd,nsd)
         xx = sum(xl(1,1:nen))/nen
         yy = sum(xl(2,1:nen))/nen

         write(iunit,"(5(f25.15,2x))") xx,yy,u(1,nel),u(2,nel)
         ! write(*,"(5(f25.15,2x))") xx,yy,u(1,nel),u(2,nel)
      end do
    
      write(iunit,*)    
   end subroutine

   !=================================================================================

   subroutine abrirArquivos()
      ! iin    = input unit number
      ! iecho  = output unit of input data
      ! iouter  = output unit of error norms

      character(len=20) :: nomeIn, nomeEcho

      iin        = 15
      iecho      = 16 
      icoords    = 18
      iconects   = 19

      iparaviewPotencial = 30
      iparaviewFluxo= 31
      iparaviewPerm = 32

      iparaview = 50

      nomeIn='inputDS.dat'
      write(*,*) " +++", nomeIn
      nomeEcho='./out/echo.dat'

      open(unit=iin,    file=nomeIn, status='old', err=100) 
      open(unit=iecho , file=nomeEcho)

      open(unit=icoords   ,file= './out/coordenadas.dat')
      open(unit=iconects  ,file= './out/conectsNodais.dat')

      return 

      100 continue
      write(*,*) ' arquivo: ', nomeIn, ' NAO encontrado'
      stop 1
       
   end subroutine abrirArquivos

   !=================================================================================

   subroutine fecharArquivos()
      close(iin      )
      close(iecho    )
      close(icoords  )
      close(iconects )
      close(ignuplotPotencial )
      close(ignuplotFluxo   )
      close(iparaviewPotencial)
      close(iparaviewFluxo)
   end subroutine fecharArquivos

   !=================================================================================

   subroutine escreverArquivoVTU_Velocidade(velocidade, nsd, numnp, tempo)
      use mGlobaisEscalares, only: zero
      use mGlobaisArranjos,  only: mat, matFratura
      use mMalha,            only: numel, numelFratura, nen, nenFratura, x, conecNodaisElem, conecNodaisFratura

      real*8,  intent(in) :: velocidade(nsd, numnp)
      integer, intent(in) :: numnp, nsd, tempo

      character(50) :: fileName, tempFormat
      integer :: node, nel, idFile, idCellType
      integer :: sinalProfundidade

      idFile = 500000000

      write(fileName, "(a, i8.8, a)") "./out/campoVel-", tempo, ".vtu"
      ! open(unit=idFile, file=trim(dir_saidaLocal)//fileName)
      open(unit=idFile, file=fileName)

      write(idFile,"(a)") '<?xml version="1.0"?>'
      write(idFile,"(a)") '<VTKFile type="UnstructuredGrid" version="0.1">'
      write(idFile,"(a)") '<UnstructuredGrid>'
      write(idFile, "(a23, i10, a17, i10, a2)") '<Piece NumberOfPoints="', numnp, '" NumberOfCells="', numel, '">'

      ! write(idFile,"(a)") '<PointData>'
      write(idFile,"(a)") '<PointData Vectors="v">'
      write(tempFormat, "(a1, i10, a12)") "(", numnp, "(es15.7,1x))"

      write(idFile,"(a, i1,a)") '<DataArray type="Float64" Name="v" NumberOfComponents="3" format="ascii">'
      do node=1,numnp
         write(idFile, tempFormat, advance="no") velocidade(1:nsd,node)!, sinalProfundidade*velocidade(nsd,node)
         if (nsd == 2) write(idFile, tempFormat, advance="no") 0d0
         write(idFile, *)
      end do
      write(idFile,"(a)") '</DataArray>'       
      ! write(idFile,"(a, i1,a)") '<DataArray type="Float32" Name="v" NumberOfComponents="', nsd,'" format="ascii">'
      ! write(idFile, tempFormat) velocidade(:,:)
      ! write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</PointData>'

      write(idFile,"(a)") '<CellData>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="materiais" format="ascii">'
      write(tempFormat, "(a1, i10, a12)") "(", numel, "(i3,1x))"
      write(idFile, tempFormat) mat(:)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</CellData>'

      write(idFile,"(a)") '<Points>'
      write(idFile,"(a)") '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      write(tempFormat, "(a1, i10, a12, a1)") "(", numnp*3, "(es15.7, 1x)", ")"
      if (plotProfundidade) then
         sinalProfundidade = -1
      else
         sinalProfundidade = 1
      end if
      select case(nsd)
         case(2)
            write(idFile, tempFormat) (x(1,node), sinalProfundidade*x(2,node), zero, node=1, numnp)
         case(3)
            write(idFile, tempFormat) (x(1:2,node), sinalProfundidade*x(3,node), node=1, numnp)
         case default
            write(*,*) "Erro na subrotina 'escreverArquivoVTU': nao implementado para a dimensao igual a ", nsd
      end select
      write(idFile,"(a)") '</DataArray>'
      write(idFile,"(a)") '</Points>'

      write(idFile,"(a)") '<Cells>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="connectivity" format="ascii">'
      ! Conectividade dos elementos da rocha intacta
      write(tempFormat, "(a1, i2, a9, a1)") "(", nen, "(i10, 1x)", ")"
      write(idFile, tempFormat) (conecNodaisElem(:, nel)-1, nel=1, numel)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="offsets" format="ascii">'
      ! Offsets dos elementos da rocha intacta
      write(tempFormat, "(a1, i10, a9, a1)") "(", numel, "(i10, 1x)", ")"
      write(idFile, tempFormat) (nen*nel, nel=1, numel)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="types" format="ascii">'
      ! Tipos dos elementos da rocha intacta
      write(tempFormat, "(a1, i10, a9, a1)") "(", numel, "(i10, 1x)", ")"
      idCellType = cellTypeVTK(nen)
      write(idFile, tempFormat) (idCellType, nel=1, numel)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</Cells>'

      write(idFile,"(a)") '</Piece>'
      write(idFile,"(a)") '</UnstructuredGrid>'
      write(idFile,"(a)") '</VTKFile>'

      close(idFile)
   end subroutine escreverArquivoVTU_Velocidade

   !=================================================================================

   subroutine escreverArquivoVTU_Pressao(pressao, ndof, numnp, tempo)
      use mGlobaisEscalares, only: zero
      use mGlobaisArranjos,  only: mat, matFratura
      use mMalha,            only: nsd, numel, numelFratura, nen, nenFratura, x, conecNodaisElem, conecNodaisFratura

      real*8,  intent(in) :: pressao(ndof, numnp)
      integer, intent(in) :: ndof, numnp, tempo

      character(50) :: fileName, tempFormat
      integer :: node, nel, idFile, idCellType
      integer :: sinalProfundidade

      idFile = 500000000

      write(fileName, "(a, i8.8, a)") "./out/campoPressao-", tempo, ".vtu"
      ! open(unit=idFile, file=trim(dir_saidaLocal)//fileName)
      open(unit=idFile, file=fileName)

      write(idFile,"(a)") '<?xml version="1.0"?>'
      write(idFile,"(a)") '<VTKFile type="UnstructuredGrid" version="0.1">'
      write(idFile,"(a)") '<UnstructuredGrid>'
      write(idFile, "(a23, i10, a17, i10, a2)") '<Piece NumberOfPoints="', numnp, '" NumberOfCells="', numel, '">'

      write(idFile,"(a)") '<PointData>'
      write(tempFormat, "(a1, i10, a12)") "(", ndof*numnp, "(es15.7,1x))"

      write(idFile,"(a, i1,a)") '<DataArray type="Float32" Name="p" NumberOfComponents="', ndof,'" format="ascii">'
      write(idFile, tempFormat) pressao(:,:)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</PointData>'

      write(idFile,"(a)") '<CellData>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="materiais" format="ascii">'
      write(tempFormat, "(a1, i10, a12)") "(", numel, "(i3,1x))"
      write(idFile, tempFormat) mat(:)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</CellData>'

      write(idFile,"(a)") '<Points>'
      write(idFile,"(a)") '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      write(tempFormat, "(a1, i10, a12, a1)") "(", numnp*3, "(es15.7, 1x)", ")"
      if (plotProfundidade) then
         sinalProfundidade = -1
      else
         sinalProfundidade = 1
      end if
      select case(nsd)
         case(2)
            write(idFile, tempFormat) (x(1,node), sinalProfundidade*x(2,node), zero, node=1, numnp)
         case(3)
            write(idFile, tempFormat) (x(1:2,node), sinalProfundidade*x(3,node), node=1, numnp)
         case default
            write(*,*) "Erro na subrotina 'escreverArquivoVTU': nao implementado para a dimensao igual a ", nsd
      end select
      write(idFile,"(a)") '</DataArray>'
      write(idFile,"(a)") '</Points>'

      write(idFile,"(a)") '<Cells>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="connectivity" format="ascii">'
      ! Conectividade dos elementos da rocha intacta
      write(tempFormat, "(a1, i2, a9, a1)") "(", nen, "(i10, 1x)", ")"
      write(idFile, tempFormat) (conecNodaisElem(:, nel)-1, nel=1, numel)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="offsets" format="ascii">'
      ! Offsets dos elementos da rocha intacta
      write(tempFormat, "(a1, i10, a9, a1)") "(", numel, "(i10, 1x)", ")"
      write(idFile, tempFormat) (nen*nel, nel=1, numel)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="types" format="ascii">'
      ! Tipos dos elementos da rocha intacta
      write(tempFormat, "(a1, i10, a9, a1)") "(", numel, "(i10, 1x)", ")"
      idCellType = cellTypeVTK(nen)
      write(idFile, tempFormat) (idCellType, nel=1, numel)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</Cells>'

      write(idFile,"(a)") '</Piece>'
      write(idFile,"(a)") '</UnstructuredGrid>'
      write(idFile,"(a)") '</VTKFile>'

      close(idFile)
   end subroutine escreverArquivoVTU_Pressao

   !=================================================================================

   subroutine escreverArquivoVTU_Fraturas_Velocidade(velocidade, nsd, numnpFratura, tempo)
      use mGlobaisEscalares, only: zero
      use mGlobaisArranjos, only: matFratura
      use mMalha, only: numelFratura, nenFratura, x, conecNodaisFratura_impressao, nosFratura_impressao_global
      use mMalha, only: normal

      real*8, intent(in)  :: velocidade(nsd, numnpFratura)
      integer, intent(in) :: nsd, numnpFratura, tempo

      character(50) :: fileName, tempFormat
      integer :: node, nel, idFile, idCellType
      integer :: sinalProfundidade

      idFile = 500000000

      write(fileName, "(a, i8.8, a)") "./out/campos_fratura_Velocidade-", tempo, ".vtu"
      ! open(unit=idFile, file=trim(dir_saidaLocal)//fileName)
      open(unit=idFile, file=fileName)

      write(idFile,"(a)") '<?xml version="1.0"?>'
      write(idFile,"(a)") '<VTKFile type="UnstructuredGrid" version="0.1">'
      write(idFile,"(a)") '<UnstructuredGrid>'
      write(idFile, "(a23, i10, a17, i10, a2)") '<Piece NumberOfPoints="', numnpFratura, '" NumberOfCells="', numelFratura, '">'

      write(idFile,"(a)") '<PointData Vectors="normais">'
      write(tempFormat, "(a1, i10, a12)") "(", numnpFratura, "(es15.7,1x))"

      write(tempFormat, "(a1, i10, a12)") "(", nsd*numnpFratura, "(es15.7,1x))"

      write(idFile,"(a, i1,a)") '<DataArray type="Float32" Name="v" NumberOfComponents="', nsd,'" format="ascii">'
      write(idFile, tempFormat) velocidade(:,:)
      write(idFile,"(a)") '</DataArray>'


      write(idFile,"(a,i1,a)") '<DataArray type="Float32" NumberOfComponents="3" Name="normais" format="ascii">'
      write(tempFormat, "(a)") "(3(es15.7,1x))"
      select case(nsd)
         case(2)
            write(idFile, tempFormat) (normal(:, nosFratura_impressao_global(node)), zero, node=1, numnpFratura)
         case(3)
            write(idFile, tempFormat) (normal(:, nosFratura_impressao_global(node)), node=1, numnpFratura)
         case default
            write(*,*) "Erro na subrotina 'escreverArquivoVTU_fraturas': nao implementado para a dimensao igual a ", nsd
      end select
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</PointData>'

      write(idFile,"(a)") '<CellData>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="materiais" format="ascii">'
      write(tempFormat, "(a1, i10, a12)") "(", numelFratura, "(i3,1x))"
      write(idFile, tempFormat) matFratura(:)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</CellData>'

      write(idFile,"(a)") '<Points>'
      write(idFile,"(a)") '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      write(tempFormat, "(a1, i10, a12, a1)") "(", numnpFratura*3, "(es15.7, 1x)", ")"
      if (plotProfundidade) then
         sinalProfundidade = -1
      else
         sinalProfundidade = 1
      end if
      select case(nsd)
         case(2)
            write(idFile, tempFormat) (x(1, nosFratura_impressao_global(node)), sinalProfundidade*x(2, nosFratura_impressao_global(node)), zero, node=1, numnpFratura)
         case(3)
            write(idFile, tempFormat) (x(1:2, nosFratura_impressao_global(node)), sinalProfundidade*x(3, nosFratura_impressao_global(node)), node=1, numnpFratura)
         case default
            write(*,*) "Erro na subrotina 'escreverArquivoVTU_fraturas': nao implementado para a dimensao igual a ", nsd
      end select
      write(idFile,"(a)") '</DataArray>'
      write(idFile,"(a)") '</Points>'

      write(idFile,"(a)") '<Cells>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="connectivity" format="ascii">'
      ! Conectividade
      write(tempFormat, "(a1, i2, a9, a1)") "(", nenFratura, "(i10, 1x)", ")"
      write(idFile, tempFormat) (conecNodaisFratura_impressao(:, nel)-1, nel=1, numelFratura)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="offsets" format="ascii">'
      ! Offsets dos elementos
      write(tempFormat, "(a1, i10, a9, a1)") "(", numelFratura, "(i10, 1x)", ")"
      write(idFile, tempFormat) (nenFratura*nel, nel=1, numelFratura)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="types" format="ascii">'
      ! Tipos dos elementos da rocha intacta
      write(tempFormat, "(a1, i10, a9, a1)") "(", numelFratura, "(i10, 1x)", ")"
      idCellType = cellTypeVTK(nenFratura)
      write(idFile, tempFormat) (idCellType, nel=1, numelFratura)

      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</Cells>'

      write(idFile,"(a)") '</Piece>'
      write(idFile,"(a)") '</UnstructuredGrid>'
      write(idFile,"(a)") '</VTKFile>'

      close(idFile)
   end subroutine escreverArquivoVTU_Fraturas_Velocidade

   !=================================================================================    

   subroutine escreverArquivoVTU_Fraturas_Pressao(pressao, salto, pressaoF, ndof, numnpFratura, tempo)
      use mGlobaisEscalares, only: zero
      use mGlobaisArranjos, only: matFratura
      use mMalha, only: nsd, numelFratura, nenFratura, x, conecNodaisFratura_impressao, nosFratura_impressao_global
      use mMalha, only: normal

      real*8, intent(in)  :: pressao(nsd,numnpFratura), salto(numnpFratura), pressaoF(numnpFratura)
      integer, intent(in) :: ndof, numnpFratura, tempo

      character(50) :: fileName, tempFormat
      integer :: node, nel, idFile, idCellType
      integer :: sinalProfundidade

      idFile = 500000000

      write(fileName, "(a, i8.8, a)") "./out/campos_fratura_pressao-", tempo, ".vtu"
      ! open(unit=idFile, file=trim(dir_saidaLocal)//fileName)
      open(unit=idFile, file=fileName)

      write(idFile,"(a)") '<?xml version="1.0"?>'
      write(idFile,"(a)") '<VTKFile type="UnstructuredGrid" version="0.1">'
      write(idFile,"(a)") '<UnstructuredGrid>'
      write(idFile, "(a23, i10, a17, i10, a2)") '<Piece NumberOfPoints="', numnpFratura, '" NumberOfCells="', numelFratura, '">'

      write(idFile,"(a)") '<PointData Vectors="normais">'
      write(tempFormat, "(a1, i10, a12)") "(", numnpFratura, "(es15.7,1x))"

      ! write(tempFormat, "(a1, i10, a12)") "(", nsd*numnpFratura, "(es15.7,1x))"

      ! write(idFile,"(a, i1,a)") '<DataArray type="Float32" Name="p1 e p2" NumberOfComponents="', nsd,'" format="ascii">'
      ! write(idFile, tempFormat) pressao(:,:)
      ! write(idFile,"(a)") '</DataArray>'
       
      write(idFile,"(a, i1,a)") '<DataArray type="Float32" Name="saltoPressao" format="ascii">'
      write(idFile, tempFormat) salto(:)
      write(idFile,"(a)") '</DataArray>'
       
      write(idFile,"(a, i1,a)") '<DataArray type="Float32" Name="pressaoF" format="ascii">'
      write(idFile, tempFormat) pressaoF(:)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a,i1,a)") '<DataArray type="Float32" NumberOfComponents="3" Name="normais" format="ascii">'
      write(tempFormat, "(a)") "(3(es15.7,1x))"
      select case(nsd)
         case(2)
            write(idFile, tempFormat) (normal(:, nosFratura_impressao_global(node)), zero, node=1, numnpFratura)
         case(3)
            write(idFile, tempFormat) (normal(:, nosFratura_impressao_global(node)), node=1, numnpFratura)
         case default
            write(*,*) "Erro na subrotina 'escreverArquivoVTU_fraturas': nao implementado para a dimensao igual a ", nsd
      end select
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</PointData>'

      write(idFile,"(a)") '<CellData>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="materiais" format="ascii">'
      write(tempFormat, "(a1, i10, a12)") "(", numelFratura, "(i3,1x))"
      write(idFile, tempFormat) matFratura(:)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</CellData>'

      write(idFile,"(a)") '<Points>'
      write(idFile,"(a)") '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      write(tempFormat, "(a1, i10, a12, a1)") "(", numnpFratura*3, "(es15.7, 1x)", ")"
      if (plotProfundidade) then
         sinalProfundidade = -1
      else
         sinalProfundidade = 1
      end if
      select case(nsd)
         case(2)
            write(idFile, tempFormat) (x(1, nosFratura_impressao_global(node)), sinalProfundidade*x(2, nosFratura_impressao_global(node)), zero, node=1, numnpFratura)
         case(3)
            write(idFile, tempFormat) (x(1:2, nosFratura_impressao_global(node)), sinalProfundidade*x(3, nosFratura_impressao_global(node)), node=1, numnpFratura)
         case default
            write(*,*) "Erro na subrotina 'escreverArquivoVTU_fraturas': nao implementado para a dimensao igual a ", nsd
      end select
      write(idFile,"(a)") '</DataArray>'
      write(idFile,"(a)") '</Points>'

      write(idFile,"(a)") '<Cells>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="connectivity" format="ascii">'
      ! Conectividade
      write(tempFormat, "(a1, i2, a9, a1)") "(", nenFratura, "(i10, 1x)", ")"
      write(idFile, tempFormat) (conecNodaisFratura_impressao(:, nel)-1, nel=1, numelFratura)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="offsets" format="ascii">'
      ! Offsets dos elementos
      write(tempFormat, "(a1, i10, a9, a1)") "(", numelFratura, "(i10, 1x)", ")"
      write(idFile, tempFormat) (nenFratura*nel, nel=1, numelFratura)
      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '<DataArray type="Int32" Name="types" format="ascii">'
      ! Tipos dos elementos da rocha intacta
      write(tempFormat, "(a1, i10, a9, a1)") "(", numelFratura, "(i10, 1x)", ")"
      idCellType = cellTypeVTK(nenFratura)
      write(idFile, tempFormat) (idCellType, nel=1, numelFratura)

      write(idFile,"(a)") '</DataArray>'

      write(idFile,"(a)") '</Cells>'

      write(idFile,"(a)") '</Piece>'
      write(idFile,"(a)") '</UnstructuredGrid>'
      write(idFile,"(a)") '</VTKFile>'

      close(idFile)
   end subroutine escreverArquivoVTU_Fraturas_Pressao

   !=================================================================================
   subroutine escreverArqUnicoVTU(pressao, iteracao_tempo)
      use mGlobaisEscalares, only: zero, numPassosTempo, deltaT
      use mMalha, only: numnp, nsd, numel, nen, x, conecNodaisElem

      real*8, intent(in) :: pressao(numnp)
      integer, intent(in) :: iteracao_tempo

      character(50) :: fileName, tempFormat
      integer :: node, nel, idFile, idCellType, i

      idFile = 500000000

      write(fileName, "(a, i6.6, a)") "campos.vtu"
      open(unit=idFile, file=fileName, position='append')

      if (iteracao_tempo == 0) then
         write(idFile, "(a)") '<?xml version="1.0"?>'
         write(idFile, "(a)") '<VTKFile type="UnstructuredGrid" version="0.1">'
         write(idFile, "(a)", advance="no") '<UnstructuredGrid TimeValues="'

         do i=0, numPassosTempo, passosPorImpressao
            write(idFile, "(es15.7,1x)", advance="no") i*deltaT
         end do
         write(idFile, "(a)") '">'

         write(idFile, "(a23, i10, a17, i10, a2)") '<Piece NumberOfPoints="', numnp, '" NumberOfCells="', numel, '">'
         write(idFile,"(a)") '<PointData>'
      end if

      write(tempFormat, "(a1, i10, a12)") "(", numnp, "(es15.7,1x))"

      write(idFile,"(a, i5.5, a)") '<DataArray type="Float32" Name="p" TimeStep="',&
               iteracao_tempo/passosPorImpressao, '" format="ascii">'
      write(idFile, tempFormat) pressao(:)
      write(idFile,"(a)") '</DataArray>'

      if ((iteracao_tempo == numPassosTempo) .or. (iteracao_tempo + passosPorImpressao > numPassosTempo)) then
         write(idFile,"(a)") '</PointData>'
         write(idFile,"(a)") '<Points>'
         write(idFile,"(a)") '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
         write(tempFormat, "(a1, i10, a12, a1)") "(", numnp*3, "(es15.7, 1x)", ")"
         select case(nsd)
            case(2)
               write(idFile, tempFormat) (x(1:2,node), zero, node=1, numnp)
            case(3)
               write(idFile, tempFormat) (x(1:3,node), node=1, numnp)
            case default
               write(*,*) "Erro na subrotina 'escreverArqUnicoVTU': nao implementado para a dimensao igual a ", nsd
         end select
         write(idFile,"(a)") '</DataArray>'
         write(idFile,"(a)") '</Points>'

         write(idFile,"(a)") '<Cells>'
         write(idFile,"(a)") '<DataArray type="Int32" Name="connectivity" format="ascii">'
         write(tempFormat, "(a1, i2, a9, a1)") "(", nen, "(i10, 1x)", ")"
         write(idFile, tempFormat) (conecNodaisElem(:, nel)-1, nel=1, numel)
         write(idFile,"(a)") '</DataArray>'
         write(idFile,"(a)") '<DataArray type="Int32" Name="offsets" format="ascii">'
         write(tempFormat, "(a1, i10, a9, a1)") "(", numel, "(i10, 1x)", ")"
         write(idFile, tempFormat) (nen*nel, nel=1, numel)
         write(idFile,"(a)") '</DataArray>'
         write(idFile,"(a)") '<DataArray type="Int32" Name="types" format="ascii">'
         write(tempFormat, "(a1, i10, a9, a1)") "(", numel, "(i10, 1x)", ")"
         idCellType = cellTypeVTK(nen)
         write(idFile, tempFormat) (idCellType, nel=1, numel)
         write(idFile,"(a)") '</DataArray>'
         write(idFile,"(a)") '</Cells>'
         write(idFile,"(a)") '</Piece>'
         write(idFile,"(a)") '</UnstructuredGrid>'
         write(idFile,"(a)") '</VTKFile>'
      end if

      close(idFile)
   end subroutine escreverArqUnicoVTU

   !=================================================================================

   function cellTypeVTK(nen)    
      use mMalha, only:nsd
       
      integer, intent(in) :: nen
      integer :: cellTypeVTK

      select case(nen)
         case(2)
            ! VTK_LINE
            cellTypeVTK = 3
         case(3)
            ! VTK_TRIANGLE
            cellTypeVTK = 5
         case(4)
            ! VTK_QUAD
            if(nsd==2) cellTypeVTK = 9
            if(nsd==3) cellTypeVTK = 10
         case(6)
            ! VTK_QUADRATIC_TRIANGLE
            cellTypeVTK = 22
         case(8)
            ! VTK_BIQUADRATIC_QUAD
            cellTypeVTK = 12
         case(9)
            ! VTK_BIQUADRATIC_QUAD
            cellTypeVTK = 28
         case default
            write(*, "(a23, i2, a24)") "Error on function cellTypeVTK: Cell type of ", nen," nodes, not implemented!"
            stop
      end select
   end function cellTypeVTK

   !=================================================================================

   subroutine escreverValoresReferencia(estrutSistEqP)
      use mEstruturasDadosSistEq, only: estruturasArmazenamentoSistemaEq
      use mMalha, only: numnp

      type(estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqP
      integer :: no, idArq

      idArq=50000000
      open(unit=idArq, file="p_ref.inc")
      write(idArq, "(a)") "*valores_potencial_ref{"
      do no=1, numnp
         write(idArq, "(2i10,es15.7)") no, 0, estrutSistEqP%u(1,no)
      end do
      write(idArq,*)
      write(idArq, "(a)") "}"
      close(idArq)
   end subroutine escreverValoresReferencia

   !=================================================================================   

   subroutine escreverArquivoParaview(estrutSistEqF, passoTempo, tempo)
      use mEstruturasDadosSistEq, only: estruturasArmazenamentoSistemaEq
      use mGlobaisEscalares, only: zero, numPassosTempo, deltaT
      use mMalha, only: nsd, numnp, numel, nen, nosFratura_impressao_global
      use mMalha, only: conecNodaisElem, conecNodaisFratura, numelFratura, numnpFratura
      use mMalha, only: numelFratura, nenFratura, x

      type(estruturasArmazenamentoSistemaEq), intent(in) :: estrutSistEqF
      integer, intent(in) :: passoTempo
      real*8, intent(in) :: tempo

      ! Lista de variaveis locais
      character(15) :: fINT1,fINT2,fFLOAT,fASCII,fBINARY,fCONC,fOFFS,fTYPES
      character(21) :: labelTempo
      character(50) :: tempFormat
         
      real*8 :: velocidadeFratura(nsd,numnpFratura)
      integer :: cont, no, i, nel, node,idCellType

      ! Variaveis de texto para criação do arquivo de vizualização
      fINT1    = 'UInt8'              ! Variavel do tipo inteira de 8 bits
      fINT2    = 'Int32'              ! Variavel do tipo inteira de 32 bits
      fFLOAT   = 'Float32'            ! Variavel do tipo float de 32 bits
      fASCII   = 'ascii'              ! Variavel alocada da forma ASCII
      fBINARY  = 'binary'             ! Variavel alocada da forma binária
      fCONC    = 'connectivity'       ! Conectividades dos nós
      fOFFS    = 'offsets'            ! Offsets das celulas
      fTYPES   = 'types'              ! Tipos das celulas

      ! Arquivo de saida do tipo VTU
      write(*,*)tipo_arq_saida
      if (tipo_arq_saida == 1) then
         call gerarLabel(labelTempo,tempo)
         open(unit=iparaview,file= './out/resultadoFluxo.vtk')
         ! cabeçalho
         write(iparaview,2000)
         write(iparaview,2010)
         write(iparaview,2020)
         write(iparaview,2030) numnp, numel

         ! Coordenadas nodais
         write(iparaview,2040)
         write(iparaview,2130) trim(adjustl(fFLOAT)),3,trim(adjustl(fASCII))
         write(iparaview, '(10x, 3es15.7)') (x(1,node), x(2,node), zero, node=1, numnp)
         write(iparaview,2101)
         write(iparaview,2041)

         ! Conectividade e tipo de celula
         write(iparaview,2050)
         write(iparaview,2110) trim(adjustl(fINT2)),trim(adjustl(fCONC)),trim(adjustl(fASCII))
         write(tempFormat, "(a1, i10, a9, a1)") "(", nen, "(i10, 1x)", ")"
         write(iparaview,tempFormat) (conecNodaisElem(:, nel)-1, nel=1, numel)
         write(iparaview,2101)
         write(iparaview,2110) trim(adjustl(fINT2)),trim(adjustl(fOFFS)),trim(adjustl(fASCII))
         write(iparaview, '(10x, i10)') (nen*nel, nel=1, numel)
         write(iparaview,2101)
         write(iparaview,2110) trim(adjustl(fINT2)),trim(adjustl(fTYPES)),trim(adjustl(fASCII))
         idCellType = cellTypeVTK(nen)
         write(iparaview, '(10x, i10)') (idCellType, nel=1, numel)
         write(iparaview,2101)
         write(iparaview,2051)

         ! Propriedades avaliadas
         ! Informações nos nós
         ! write(iparaview,2060)
         ! write(iparaview,2061)

         ! Informações nos elementos
         ! write(iparaview,2070)
         ! write(iparaview,2071)

         ! fechamento do arquivo
         write(iparaview,2031)
         write(iparaview,2021)
         write(iparaview,2011)
         close(iparaview)

         call escreverArquivoVTU_Velocidade(estrutSistEqF%u, nsd, numnp, passoTempo)
             
         if (numelFratura > 0)  then
            do i=1,numnpFratura
               no = nosFratura_impressao_global(i)
               velocidadeFratura(:,i) = estrutSistEqF%u(:, no)
            enddo
            call escreverArquivoVTU_Fraturas_Velocidade(velocidadeFratura, nsd, numnpFratura, passoTempo)
         endif   
      ! Arquivo de saida do tipo VTK
      elseif (tipo_arq_saida == 0) then

         ! if(passoTempo==1) then
         !    call escreverArqParaviewVector(estrutSistEqF%u, estrutSistEqF%ndof, numnp, nen, &
         !             conecNodaisElem, 2, trim(labelTempo), len(trim(labelTempo)), iparaviewFluxo)  
         ! else
         !    call escreverArqParaviewIntermed(estrutSistEqF%u, estrutSistEqF%ndof, numnp, 'velocidade',   &
         !             trim(labelTempo), len(trim(labelTempo)), iparaviewFluxo)
         ! endif
      end if
      ! Lista de formato de leitura  
      2000 format('<?xml version="1.0"?>') 
      2010 format('<VTKFile type="UnstructuredGrid" version="0.1">')
      2011 format('</VTKFile>') 
      2020 format(2x,'<UnstructuredGrid>')
      2021 format(2x,'</UnstructuredGrid>') 
      2030 format(4x,'<Piece NumberOfPoints="',i0,'" NumberOfCells="',i0,'">')
      2031 format(4x,'</Piece>')

      2040 format(6x,'<Points>')
      2041 format(6x,'</Points>') 
      2050 format(6x,'<Cells>')
      2051 format(6x,'</Cells>')

      2060 format(6x,'<PointData>')
      2061 format(6x,'</PointData>') 
      2070 format(6x,'<CellData>')
      2071 format(6x,'</CellData>') 
      2101 format(8x,'</DataArray>')
      ! Para utilização (de um valor escalar/propriedades da celula
      2110 format(8x,'<DataArray type="',A,'" Name="',A,'" format="',A,'">')
      2111 format(10x,i5)
      ! Para utilização de um vetor
      2120 format(8x,'<DataArray type="',A,'" Name="',A,'" NumberOfComponents="',i0,'" format="',A,'">')
      2121 format(10x,2i5)
      ! Para utilização das coordenadas nodais
      2130 format(8x,'<DataArray type="',A,'" NumberOfComponents="',i0,'" format="',A,'">')
      2131 format(10x,i0)

      3000 format(10x,e15.8,2e18.8)

   end subroutine

   !=================================================================================
   ! subroutine escreverVTUcabecalho(pressao, salto, pressaoF, ndof, numnpFratura, tempo)
   !    use mGlobaisEscalares, only: zero
   !    use mGlobaisArranjos, only: matFratura
   !    use mMalha, only: nsd, numelFratura, nenFratura, x, conecNodaisFratura_impressao, nosFratura_impressao_global
   !    use mMalha, only: normal

   !    real*8, intent(in)  :: pressao(nsd,numnpFratura), salto(numnpFratura), pressaoF(numnpFratura)
   !    integer, intent(in) :: ndof, numnpFratura, tempo
   ! end subroutine
end module

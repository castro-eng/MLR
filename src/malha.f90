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
module mMalha
   implicit none

   integer*4 :: nsd, numnp
   real*8, allocatable :: x(:,:)

   ! Variaveis relativas ao maciço
   integer*4 :: numel, nen
   integer*4, allocatable :: conecNodaisElem(:,:)

   ! Variaveis relativas as descontinuidades
   integer*4 :: numelFratura, numnpFratura, nenFratura
   integer*4, allocatable :: conecNodaisFratura(:,:), conecNodaisReduzido(:,:)
   ! A variável conecNodaisFratura contém as conectividades dos elementos reduzidos em relação ao domínio integro
   ! A variável conecNodaisReduzido contém as conectividades dos elementos reduzido em relação ao domínio reduzido

   ! Variaveis relativas ao bordo
   integer*4 :: numelBordas, nenFraturasNaBorda
   integer*4, allocatable :: conecNodaisBordas(:,:)

   ! Variaveis mistas
   integer*4, allocatable :: conectividadesInterseccao(:,:), conectividadesInterseccaoFLuxCons(:,:,:)
   integer*4, allocatable :: listaDosElemsPorNo(:,:)
   integer*4, pointer  :: listaDosElemsPorNoCSR(:,:) => null()
   
   integer*4 :: numConexoesPorElem 
   real*8, allocatable    :: normal(:,:), normalElem(:,:), normalElemBorda(:,:)
   integer*4, allocatable :: nosFratura_impressao_global(:), conecNodaisFratura_impressao(:,:)

   ! Variaveis relativas ao fluxo consistente
   integer*4 :: numelDescontNaBorda, numelBordasDirichlet
   integer*4, allocatable :: conecNodaisDesconNaBorda(:,:), conecNodaisBordasDirichlet(:,:)

   ! integer*4 :: nenFraturasNaBorda, nenBorda
   ! integer*4 :: numelFraturasNaBorda, numelBordasDirichlet
   ! integer*4 :: numelVizDirFrat, numelVizEsqFrat, numnpVizDirFrat, numnpVizEsqFrat
   ! integer*4, allocatable :: conecNodaisFraturasNaBorda(:,:), conecNodaisBordasDirichlet(:,:)
   ! integer*4, allocatable :: posicaoElemEmRelacaoA_Fratura(:)
   ! integer*4, allocatable :: nosVizDirFrat_impressao_global(:), conecNodaisVizDirFrat_impressao(:,:)
   ! integer*4, allocatable :: nosVizEsqFrat_impressao_global(:), conecNodaisVizEsqFrat_impressao(:,:)

   ! funcoes e subrotinas
   public :: genfl, genfl1
   public :: gensh, gensh1, gensh2, gensh3, igen
   public :: conectIntegroParaReduzido

   contains

   !=======================================================================
   subroutine conectIntegroParaReduzido()
      ! A rotina serve para transcrever as conectividades dos elementos reduzidos
      ! descritos nos nós da malha integra, para uma conectividade para elementos em
      ! uma malha reduzida

      implicit none

      logical :: flag
      integer :: i,j,k,pos_conc, tam_array
      integer, allocatable :: array_mapeamento(:)

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

      ! Montagem e escrita das conectividades
      conecNodaisReduzido(:,:)=0
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
                  conecNodaisReduzido(j,i)=k
               end if
               k=k+1
            end do
         end do
      end do
      deallocate(array_mapeamento)

      RETURN
   end subroutine

   !=======================================================================

   subroutine genfl(a,nra,iin)
      !.... program to read and generate floating-point nodal data                                                                                        
      !         a       = input array                                         
      !         nra     = number of rows in a (le.6)                          
      !         n       = node number                                         
      !         numgp   = number of generation points                         
      !         ninc(i) = number of increments for direction i                
      !         inc(i)  = increment for direction i                           
                                                                            
      use mGlobaisEscalares  

      implicit none
                                                                      
      integer*4:: nra, iin                                               
      real*8  :: a(nra,*)

      integer*4:: n,numgp,ninc(3),inc(3)   
      integer*4:: i, j, m, mgen
      integer*4, parameter :: dim1 = 6, dim2 = 20
      real*8  :: temp(dim1,dim2)
      !real*8  :: temp(6,20)
                                                                      
      100 continue                                                          
         read(iin,1000) n,numgp,(temp(i,1),i=1,nra)      
         ! write(*,1000) n,numgp,(temp(i,1),i=1,nra)      

         if (n.eq.0) return                                                
         ! call move(a(1,n),temp,nra)           
         a(1:nra,n) = temp(1:nra,1)

         if (numgp.ne.0) then
            do 200 j=2,numgp
               read(iin,1000) m,mgen,(temp(i,j),i=1,nra)    

               ! if (mgen.ne.0) call move(temp(1,j),a(1,m),nra) 
               if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m) 

            200 continue
            read(iin,2000) (ninc(i),inc(i),i=1,3)
            call genfl1(a,nra, temp, n, numgp, ninc, inc)                                                
         endif
      go to 100                                                         
                                                                      
      ! 1000 format(2i10,6f10.0)                                                
      1000 format(2i10,6f10.0)                                                
      2000 format(16i10)                                                      
   end  subroutine

   !******************************************************************************

   subroutine genfl1(a,nra, temp, n, numgp, ninc, inc)     
      !.... program to generate floating-point nodal data 
      !        via isoparametric interpolation         
      !         iopt = 1, generation along a line                             
      !              = 2, generation over a surface                           
      !              = 3, generation within a volume                            
      use mGlobaisEscalares

      implicit none                                          

      integer*4:: nra
      real*8  :: a(nra,*)
      ! real*8  :: temp(6,20)
      integer*4, parameter :: dim1 = 6, dim2 = 20, um = 1
      real*8  :: temp(dim1,dim2)
      integer*4:: n,numgp,ninc(3),inc(3)                 

      real*8  :: sh(20), dr, ds, dt, r, s, t
      integer*4:: iopt, ni, nj, nk, ii, jj, kk, i, j, k

      iopt = 3                                                          
      if (ninc(3).eq.0) iopt = 2                                        
      if (ninc(2).eq.0) iopt = 1                                        
                                                                      
      dr = zero                                                           
      ds = zero                                                           
      dt = zero                                                           
                                                                      
      if (ninc(1).ne.0) dr = two/ninc(1)                                
      if (ninc(2).ne.0) ds = two/ninc(2)                                
      if (ninc(3).ne.0) dt = two/ninc(3)                                
                                                                      
      ii = ninc(1)+1                                                    
      jj = ninc(2)+1                                                    
      kk = ninc(3)+1                                                    
                                                                      
      ni = n                                                            
      nj = n                                                            
      nk = n                                                            
                                                                      
      t = -one
      do 300 k=1,kk

         s = -one
         do 200 j=1,jj

            r = -one
            do 100 i=1,ii
               call gensh(r,s,t,sh,numgp,iopt)
               call multab(temp,sh,a(1,ni),dim1,dim2,nra,numgp,nra,um,um)
               ni = ni + inc(1)
               r = r + dr
            100 continue
                                                                      
            nj = nj + inc(2)
            ni = nj
            s = s + ds
         200 continue
                                                                      
         nk = nk + inc(3)
         ni = nk
         t = t + dt
      300 continue
                                                                      
      return                                                            
   end  subroutine

   !****************************************************************************** 
  
   subroutine gensh(r,s,t,sh,numgp,iopt)                             
      !.... program to call shape function routines         
      !        for isoparametric generation         
                                                                      
      implicit none                                      

      real*8  :: r, s, t, sh(*)                                                        
      integer*4:: numgp, iopt
                                                                      
      go to (100,200,300),iopt                                                
                                                                      
         100 call gensh1(r,sh,numgp)                                           
         return                                                            
                                                                      
         200 call gensh2(r,s,sh,numgp)                                         
         return                                                            
                                                                      
         300 call gensh3(r,s,t,sh,numgp)                                       
         return                                                            
                                                                      
   end subroutine

   !******************************************************************************

   subroutine gensh1(r,sh,n)  
      !.... program to compute 1d shape functions           
      !        for isoparametric generation                     

      use mGlobaisEscalares
                                                                      
      implicit none

      real*8  :: r, sh(*)                                                   
      integer*4:: n
                                                                      
      sh(2) = pt5*r                                                       
      sh(1) = pt5 - sh(2)                                                   
      sh(2) = pt5 + sh(2)                                                   
      if (n.eq.3) then
         sh(3) = one - r*r                                                     
         sh(1) = sh(1) - pt5*sh(3)                                             
         sh(2) = sh(2) - pt5*sh(3)                                             
      endif
                                                                      
      return                                                            
   end subroutine

   !******************************************************************************

   subroutine gensh2(r,s,sh,n)     
      !.... program to compute 2d shape functions 
      !        for isoparametric generation    
                                                                      
      use mGlobaisEscalares
      implicit none                                         

      real*8  :: r, s, sh(*)                                                   
      integer*4:: n    

      real*8  :: r1, r2, r3, s1, s2, s3

      r2 = pt5*r                                                          
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s                                                          
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      sh(1) = r1*s1                                                       
      sh(2) = r2*s1                                                       
      sh(3) = r2*s2                                                       
      sh(4) = r1*s2                                                       
      if (n.eq.4) return                                                
                                                                      
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      sh(5) = r3*s1                                                       
      sh(6) = s3*r2                                                       
      sh(7) = r3*s2                                                       
      sh(8) = s3*r1                                                       
      sh(1) = sh(1) - pt5*(sh(5) + sh(8))
      sh(2) = sh(2) - pt5*(sh(6) + sh(5))
      sh(3) = sh(3) - pt5*(sh(7) + sh(6))
      sh(4) = sh(4) - pt5*(sh(8) + sh(7))
                                                                      
      return
   end subroutine

   !******************************************************************************

   subroutine gensh3(r,s,t,sh,n) 
      !.... program to compute 3d shape functions            
      !        for isoparametric generation   
                                                                      
      use mGlobaisEscalares
      implicit none                                         

      real*8  :: r, s, t, sh(*)                                                   
      integer*4:: n    

      real*8  :: r1, r2, r3, rs1, rs2, rs3, rs4
      real*8  :: s1, s2, s3, t1, t2, t3                                              
                                                                      
      r2 = pt5*r
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      t2 = pt5*t
      t1 = pt5 - t2                                                         
      t2 = pt5 + t2                                                         
                                                                      
      rs1 = r1*s1                                                         
      rs2 = r2*s1                                                         
      rs3 = r2*s2                                                         
      rs4 = r1*s2                                                         
      sh(1) = rs1*t1                                                      
      sh(2) = rs2*t1                                                      
      sh(3) = rs3*t1                                                      
      sh(4) = rs4*t1                                                      
      sh(5) = rs1*t2                                                      
      sh(6) = rs2*t2                                                      
      sh(7) = rs3*t2                                                      
      sh(8) = rs4*t2                                                      
      if (n.eq.8) return                                                 
                                                                      
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      t3 = one - t*t                                                        
      sh(17) = t3*rs1                                                     
      sh(18) = t3*rs2                                                     
      sh(19) = t3*rs3                                                     
      sh(20) = t3*rs4                                                     
      rs1 = r3*s1                                                         
      rs2 = s3*r2                                                         
      rs3 = r3*s2                                                         
      rs4 = s3*r1                                                         
      sh( 9) = rs1*t1                                                     
      sh(10) = rs2*t1                                                     
      sh(11) = rs3*t1                                                     
      sh(12) = rs4*t1                                                     
      sh(13) = rs1*t2                                                     
      sh(14) = rs2*t2                                                     
      sh(15) = rs3*t2                                                     
      sh(16) = rs4*t2                                                     
                                                                      
      sh(1) = sh(1) - pt5*(sh( 9) + sh(12) + sh(17))
      sh(2) = sh(2) - pt5*(sh( 9) + sh(10) + sh(18))
      sh(3) = sh(3) - pt5*(sh(10) + sh(11) + sh(19))
      sh(4) = sh(4) - pt5*(sh(11) + sh(12) + sh(20))
      sh(5) = sh(5) - pt5*(sh(13) + sh(16) + sh(17))
      sh(6) = sh(6) - pt5*(sh(13) + sh(14) + sh(18))
      sh(7) = sh(7) - pt5*(sh(14) + sh(15) + sh(19))
      sh(8) = sh(8) - pt5*(sh(15) + sh(16) + sh(20))
                                                                      
      return                                                            
   end subroutine                                                              

   !**** new **********************************************************************

   subroutine igen(ia,m,iin)
      !.... program to read and generate integer*4, nodal data
      !        ia = input array
      !         m = number of rows in ia
      !         n = node number
      !        ne = end node in generation sequence
      !        ng = generation increment

      use mGlobaisEscalares
      integer*4:: m, ia(m,*), iin

      integer*4:: ib(13)
      integer*4:: n, ne, ng
      integer*4:: i

      100 continue
         read(iin,1000) n,ne,ng,(ib(i),i=1,m)
         ! write(*,1000) n,ne,ng,(ib(i),i=1,m)
         if (n.eq.0) return

         if (ng.eq.0) then
            ne = n
            ng = 1
         else
            ne = ne - mod(ne-n,ng)
         endif

         do 200 i=n,ne,ng
            ! call imove(ia(1,i),ib,m)
            ia(1:m,i)=ib(1:m)
         200 continue

      go to 100

      1000 format(16i10)
   end subroutine

   !**** new **********************************************************************

   subroutine local(conectElem,x,xl,nen,nrowx,nrowxl)
      !.... program to localize a global array
      !        note: it is assumed nrowxl.le.nrowx

      implicit none

      integer*4:: conectElem(*)
      integer*4:: nrowx, nrowxl, nen
      double precision :: x(nrowx,*),xl(nrowxl,*)

      integer*4:: i, j, node

      do 200 j=1,nen
         node = conectElem(j)

         do 100 i=1,nrowxl
            xl(i,j)= x(i,node)
         100 continue

      200 continue

      return
   end subroutine

   !**** new **********************************************************************

   subroutine multab(a,b,c,ma,mb,mc,l,m,n,iopt)
      !.... program to multiply two matrices
      !        l = range of dot-product index
      !        m = number of active rows in c
      !        n = number of active columns in c

      implicit none

      integer*4:: ma,mb,mc,l,m,n,iopt
      real*8  :: a(ma,*),b(mb,*),c(mc,*)

      integer*4:: i,j

      go to (1000,2000,3000,4000),iopt

         !.... iopt = 1, c(i,j) = a(i,k)*b(k,j) , (c = a * b)
         1000 do 1200 i=1,m

            do 1100 j=1,n
               c(i,j) = rcdot(a(i,1),b(1,j),ma,l)
            1100 continue

         1200 continue
         return

         !.... iopt = 2, c(i,j) = a(k,i)*b(k,j) (c = a  * b)
         2000 do 2200 i=1,m
            do 2100 j=1,n
               c(i,j) = dot_product (a(:,i),b(:,j)) !,l)
               ! c(i,j) = coldot(a(1,i),b(1,j),l)
            2100 continue
         2200 continue
         return

         !.... iopt = 3, c(i,j) = a(i,k)*b(j,k) (c = a * b )
         3000 do 3200 i=1,m
            do 3100 j=1,n
               c(i,j) = rowdot(a(i,1),b(j,1),ma,mb,l)
            3100 continue
         3200 continue
         return

         !.... iopt = 4, c(i,j) = a(k,i)*b(j,k) (c = a  * b )
         4000 do 4200 i=1,m
            do 4100 j=1,n
               c(i,j) = rcdot(b(j,1),a(1,i),mb,l)
            4100 continue
         4200 continue
         return

   end subroutine

   !**** new **********************************************************************
   function rcdot(a,b,ma,n)
      !.... program to compute the dot product of a vector stored row-wise
      !        with a vector stored column-wise

      implicit none

      integer*4:: ma, n
      real*8  :: a(ma,*),b(*)

      real*8  :: rcdot
      integer*4:: i

      rcdot = 0.0d0

      do 100 i=1,n
         rcdot = rcdot + a(1,i)*b(i)
      100 continue

      return
   end function

   !**** new **********************************************************************

   function rowdot(a,b,ma,mb,n)
      !.... program to compute the dot product of vectors stored row-wise

      implicit none

      integer*4:: ma, mb, n
      real*8  :: a(ma,*),b(mb,*)

      real*8  :: rowdot
      integer*4:: i

      rowdot = 0.0d0

      do 100 i=1,n
         rowdot = rowdot + a(1,i)*b(1,i)
      100 continue

     return
   end function

   !=======================================================================
    
   ! subroutine criarListaVizinhos(nen,numnp,numel,conecElem,listaDosElemsPorNo)
   !    ! Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no

   !    implicit none

   !    integer*4 :: nen,numnp,numel
   !    integer*4, dimension(nen,numel)  :: conecElem
   !    integer*4, dimension(nen,numnp)  :: listaDosElemsPorNo
   !    integer*4:: no,nel,l
      
   !    do nel=1, numel
   !       do l=1, nen
   !          no=conecElem(l,nel)
   !          listaDosElemsPorNo(l,no) = nel
   !       end do
   !    end do
    
   ! end subroutine
      
   !=======================================================================
   subroutine criarListaVizinhosCRS(nen,numnp,numel,conecElem,nVizinMax)
      ! Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no

      implicit none

      integer nen,numnp,numel,nVizinMax
      integer, dimension(nen,numel)  :: conecElem
      integer, allocatable  :: listaDosElemsPorNo(:,:)!(nVizinMax,numnp)
      integer :: no,nel,l,numVizAtual, i, j, cont, noElemento
      integer :: listaNumVizPorNo(numnp)

      listaNumVizPorNo=0
      nVizinMax=0
      do nel=1, numel
         do l=1, nen
            cont=0
            no=conecElem(l,nel)
            listaNumVizPorNo(no)=listaNumVizPorNo(no)+1
          end do
      end do

      nVizinMax=maxval(listaNumVizPorNo)
      
      print*, "nVizinMax=", nVizinMax
      
      allocate(listaDosElemsPorNo(nVizinMax+1, numnp))
      allocate(listaDosElemsPorNoCSR(nVizinMax, numnp))
     
      listaDosElemsPorNo = 0 
      do nel=1, numel
         do l=1, nen
            no=conecElem(l,nel)
            listaDosElemsPorNo(nVizinMax+1,no) = listaDosElemsPorNo(nVizinMax+1,no) + 1 
            numVizAtual                        = listaDosElemsPorNo(nVizinMax+1,no)
            listaDosElemsPorNo(numVizAtual,no) = nel
         end do
      end do
          
      do no=1, numnp
            listaDosElemsPorNoCSR(:,no)=listaDosElemsPorNo(1:nVizinMax,no)
      end do
      deallocate(listaDosElemsPorNo)

      ! do no=1, numnp
      !    print*,no, ",", listaDosElemsPorNoCSR(:,no)
      ! end do
      ! stop

      ! do no=1, numnp
      !    numVizAtual = listaDosElemsPorNo(nVizinMax,no)  
      !    write(*,'(8(1X,I3))') no, numVizAtual , listaDosElemsPorNo(1:numVizAtual,no)
      ! end do
      ! write(*,*) '---------------'
      ! stop

   end subroutine

   !=======================================================================

   subroutine calcularNenFratura
      select case(nen)
         case(4)
            if(nsd==2)nenFratura = 2 !elemento quadrilatero, fratura linha
            if(nsd==3)nenFratura = 3 !elemento tetraedrico, fratura triangulo
         case(8)
            nenFratura = 4
         case(9)
            nenFratura = 3
         case(16)
            nenFratura = 4
         case(25)
            nenFratura = 5
         case(36)
            nenFratura = 6
         case(49)
            nenFratura = 7
         case(64)
            nenFratura = 8
         case(3) ! elemento triangular, fratura linha
            nenFratura = 2 
         case(6)
            nenFratura = 3
      end select
   end subroutine

   !======================================================================= 

   subroutine setNormalElem(normal, xl)
      ! Set the element normal vector (only linear and bilinear elements).
      use mGlobaisEscalares, only: zero

      real*8, intent(in) :: xl(nsd, nenFratura)
      real*8, intent(out) :: normal(nsd)
      real*8 :: length, diff(2), tangent(2)
      real*8 :: v1(3), v2(3)

      if (nsd == 2) then
         diff = xl(:,2) - xl(:,1)
         length = sqrt(diff(1)*diff(1) + diff(2)*diff(2))
         tangent(:) = diff(:)/length

         ! Tangent must have positive y-coordinate. It's an arbitrary
         ! choice to assure that setNormal(normal, pointA, pointB) and
         ! setNormal(normal, pointB, pointA) yield the same result
         if (tangent(2) < 0d0) then
            tangent(:) = -tangent(:)
         end if

         normal(1) = tangent(2)
         normal(2) = -tangent(1)

         ! When the line is approximately parallel to the x-axis, the normal vector must point upward.
         if ((abs(normal(1)) < 1d-10) .and. (normal(2) < 0d0)) then
            normal(:) = -normal(:)
         end if
      else if (nsd == 3) then
         ! v1 = xl(:,1) - xl(:,2)
         v1 = xl(:,2) - xl(:,1)
         v2 = xl(:,3) - xl(:,2)

         ! Cross product
         normal(1) = v1(2)*v2(3) - v1(3)*v2(2)
         normal(2) = v1(3)*v2(1) - v1(1)*v2(3)
         normal(3) = v1(1)*v2(2) - v1(2)*v2(1)

         length = sqrt(dot_product(normal(:), normal(:)))
         normal(:) = normal(:)/length

         ! Normal must have non-negative z-value (or a negative near-zero value)
         ! This choice is arbitrary and intend to let the normal vector independent of the nodes' order
         if (normal(3) < -1d-10) then
            normal(:) = -normal(:)
         end if
      else
         print*, "Erro na subrotina setNormalElem: so e aceito nsd = 2 ou nsd = 3"
         stop
      end if
   end subroutine setNormalElem

   !=======================================================================

   subroutine gerarNormais
      integer :: nel, nosLocais(nenFratura), no, i
      real*8 :: xl(nsd, nenFratura), normalTemp(nsd)

      if (.not. allocated(normal)) then
         allocate(normal(nsd,numnp))
      end if

      do nel=1,numelFratura
         ! call local(conecNodaisFratura, x, xl, nenFratura, nsd, nsd)
         call local(conecNodaisFratura(:,nel), x, xl, nenFratura, nsd, nsd)
         call setNormalElem(normalTemp, xl)
         do i=1,nenFratura
            no = conecNodaisFratura(i,nel)
            normal(:,no) = normalTemp(:)
            ! write(*,'()'), nel, no, normal(:,no)
         end do
      end do
   end subroutine gerarNormais

   !=======================================================================
 
   subroutine gerarNormaisElemFratura
      integer :: nel
      real*8 :: xl(nsd,nenFratura)

      if (.not. allocated(normalElem)) then
         allocate(normalElem(nsd,numelFratura))
      end if

      do nel=1,numelFratura
         call local(conecNodaisFratura(:,nel), x, xl, nenFratura, nsd, nsd)
         ! call local(conecNodaisFratura, x, xl, nenFratura, nsd, nsd)
         call setNormalElem(normalElem(:,nel), xl)
      end do
   end subroutine gerarNormaisElemFratura

   !=======================================================================

   subroutine gerarNormaisElemBordas

      use mGlobaisArranjos, only: tiposElemBordas

      integer :: nel, tipoBorda
      real*8 :: xl(nsd,nenFratura)
      real*8 :: vetorI(nsd), vetorJ(nsd), vetorK(nsd), vetor(nsd)

      if (.not. allocated(normalElemBorda)) then
         allocate(normalElemBorda(nsd,numelBordas))
      end if
   
      ! vetorI(:) = 0d0
      ! vetorJ(:) = 0d0
      ! vetorK(:) = 0d0
   
      ! vetorI(1) = 1d0
      ! vetorJ(2) = 1d0
      ! vetorK(nsd) = 1d0
   
      do nel=1,numelBordas
   
         tipoBorda = tiposElemBordas(nel)
         if(tipoBorda==1) then
            normalElemBorda(1,nel)= 0.d0
            normalElemBorda(2,nel)=-1.d0
            if(nsd==3) normalElemBorda(3,nel)=0.d0
         elseif(tipoBorda==2) then
            normalElemBorda(1,nel)= 1.d0
            normalElemBorda(2,nel)= 0.d0
            if(nsd==3) normalElemBorda(3,nel)=0.d0
         elseif(tipoBorda==3) then
            normalElemBorda(1,nel)= 0.d0
            normalElemBorda(2,nel)= 1.d0
            if(nsd==3) normalElemBorda(3,nel)=0.d0
         elseif(tipoBorda==4) then
            normalElemBorda(1,nel)=-1.d0
            normalElemBorda(2,nel)= 0.d0
            if(nsd==3) normalElemBorda(3,nel)=0.d0
         elseif(tipoBorda==5) then
            normalElemBorda(1,nel)= 0.d0
            normalElemBorda(2,nel)= 0.d0
            if(nsd==3) normalElemBorda(3,nel)=-1.d0
         else
            normalElemBorda(1,nel)= 0.d0
            normalElemBorda(2,nel)= 0.d0
            if(nsd==3) normalElemBorda(3,nel)= 1.d0
         endif

         ! call local(conecNodaisBordas, x, xl, nenFratura, nsd, nsd)
         ! call setNormalElem(normalElemBorda(:,nel), xl)
      
         ! tipoBorda = tiposElemBordas(nel)
         ! if (tipoBorda==1) then
         !    vetor(:) = -vetorJ(:)
         ! else if (tipoBorda==2) then
         !    vetor(:) = vetorI(:)
         ! else if (tipoBorda==3) then
         !    vetor(:) = vetorJ(:)
         ! else
         !    vetor(:) = -vetorI(:)
         ! end if
      
         ! if (dot_product(normalElemBorda(:,nel), vetor(:)) < 0d0) then
         !    normalElemBorda(:,nel) = -normalElemBorda(:,nel)
         ! end if
      end do
   end subroutine gerarNormaisElemBordas

   !=======================================================================
 
   subroutine gerarNormaisElemBordas_old
      integer :: nel
      real*8 :: xl(nsd,nenFratura)

      if (.not. allocated(normalElemBorda)) then
         allocate(normalElemBorda(nsd,numelBordas))
      end if

      do nel=1,numelBordas
         call local(conecNodaisBordas, x, xl, nenFratura, nsd, nsd)
         call setNormalElem(normalElemBorda(:,nel), xl)
      end do
   end subroutine gerarNormaisElemBordas_old  

   !=======================================================================
 
   subroutine gerarNumnpFratura
      ! nosFratura_global_impressao dado um no global o vetor retorna o no de impressao correspondente
      integer :: noLocal, noGlobal, nel, nosFratura_global_impressao(numnp)

      allocate(conecNodaisFratura_impressao(nenFratura, numelFratura))

      nosFratura_global_impressao(:) = 0
      numnpFratura = 0
      do nel = 1, numelFratura
         do noLocal=1,nenFratura
            noGlobal = conecNodaisFratura(noLocal, nel)
            if (nosFratura_global_impressao(noGlobal) == 0) then
               numnpFratura = numnpFratura + 1
               nosFratura_global_impressao(noGlobal) = numnpFratura
            end if
            conecNodaisFratura_impressao(noLocal, nel) = nosFratura_global_impressao(noGlobal)
         end do
      end do

      allocate(nosFratura_impressao_global(numnpFratura))

      do noGlobal=1,numnp
         if (nosFratura_global_impressao(noGlobal) > 0) then
            nosFratura_impressao_global(nosFratura_global_impressao(noGlobal)) = noGlobal
         end if
      end do
   end subroutine gerarNumnpFratura

end module

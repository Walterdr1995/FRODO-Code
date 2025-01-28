Program FRODO

use omp_lib ! import the library for parallel computing

implicit none
integer, parameter :: dp = 8

real(kind = dp), parameter :: alpha = 2.0 ! fractional parameter

real(kind = dp), parameter :: x_min = -40.0 !lower x-boundary of the domain
real(kind = dp), parameter :: x_max = 40.0  !upper x-boundary of the domain

real(kind = dp), parameter :: s_min = -3.0  !lower s=ln(p/p_0) of the domain
real(kind = dp), parameter :: s_max = 7.0   !upper s=ln(p/p_0) of the domain

real(kind = dp), parameter :: t_min = 0.0   !starting time of the computation
real(kind = dp), parameter :: t_max = 400.0 !ending time of the computation

integer :: NX,NS, NT  !points in x,s and t

write(*,*) 'x-axis = [',x_min,',',x_max,']'
write(*,*) 'Discretisation of the x-axis:'
!read(*,*) NX
NX = 2000  !see line 19

write(*,*) 's-axis = [',s_min,',',s_max,']'
write(*,*) 'Discretisation of the s-axis:'
NS = 600 !see line 19

write(*,*) ' '

write(*,*) 't-axis = [',t_min,',',t_max,']'
write(*,*) 'Discretisation of the t-axis:'
!read(*,*) NT
NT = 100001 !see line 19

call Hauptprogramm(NX,NT,NS,x_min,x_max,t_min,t_max,alpha,s_min,s_max)  ! Here, the main Programm is called and the relevant parameters are delivered
End Program

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine Hauptprogramm(NX,NT,NS,x_min,x_max,t_min,t_max,alpha,s_min,s_max)  ! main routine of the programm
implicit none
integer, parameter :: dp = 8
! Here, the main parameters of the system are defined again, for an explanation see the first paragraph of the code
real(kind = dp):: alpha

real(kind = dp):: x_min 
real(kind = dp):: x_max

real(kind = dp):: s_min 
real(kind = dp):: s_max

real(kind = dp):: t_min 
real(kind = dp):: t_max

integer :: NX, NT, NS

!---------- Passed on values ------------------------

!functions
real(kind=dp)   :: d_plus
real(kind=dp)   :: d_minus
real(kind=dp)   :: q
real(kind=dp)   :: x_bl, x_bh
real(kind=dp)   :: s_bl, s_bh
real(kind=dp)   :: f_initial
real(kind=dp)   :: konvektion
real(kind=dp)   :: konvektionderv

!Here the necessary arrays for the main computation takes place. Explanations
! for the individual arrays below:

real(kind=dp), dimension (0:NX,0:NS) :: u ! main array, denoting our solution f

real(kind=dp), dimension (0:NX,0:NS) :: rhs !A*u = u_rhs
real(kind=dp), dimension (0:NS,0:NX) :: rhs_sad !B*u = u_rhs_sad
real(kind=dp), dimension (0:NX) :: rhs_pr 
real(kind=dp), dimension (0:NS) :: rhs_sad_pr 

real(kind=dp), dimension (0:NX,0:NX) :: A_i !Matrix for the x advancements
real(kind=dp), dimension (0:NS,0:NS) :: B! Matrix for the s advancements

real(kind=dp), dimension (0:NX,0:NX) ::Proxy_NX !Those Matrices are simple placeholders for some methods
real(kind=dp), dimension (0:NS,0:NS) ::Proxy_NS

real(kind=dp), dimension (0:NX,0:NX) :: A_invers ! Inverse matrices to A_i and B
real(kind=dp), dimension (0:NS,0:NS,0:NX) :: B_invers


!x and t arrays
real(kind=dp), dimension (0:NX) :: x
real(kind=dp), dimension (0:NS) :: s
real(kind=dp), dimension (0:NT) :: t

!Grünwald-factors
real(kind=dp), dimension (0:NX) :: g_alpha

!Matrix entries
!real(kind=dp), dimension (0:NX) :: spe, sme !spe=sigma_plus_explizit , m = minus
real(kind=dp), dimension (0:NX) :: spi, smi 
real(kind=dp), dimension (0:NX) :: k,k0,k1
real(kind=dp), dimension (0:NX) :: adi


!step-widths
real(kind=dp) :: dx, dy, dt ,ds

!loop variables
integer :: ZeitSchritt
integer :: i,j,l,m

!fractional-factor
real(kind=dp) :: fracFac ! this factor will later be defined to take care of the pre-factors of the 
                         ! Riesz derivative in comparison with the Liouville derivative

!--------------------------------------------------------------


dx = (x_max-x_min)/NX ! assign a value to the step-widths
dt = (t_max-t_min)/NT
ds = (s_max-s_min)/NS

fracFac = -0.5/dcos(alpha/2.*3.14159265) ! assign a value to the fractional factor

!assign the value to the x-array
do i = 0, NX
   x(i) = x_min + i*dx
end do

!assign the value to the s-array
do l = 0, NS
   s(l) = s_min + l*ds
end do

!assign the value to the t-array
do i = 0, NT
   t(i) = t_min + i*dt
end do

! assign values to the grünwald coefficients
g_alpha(0) = 1.0
do i = 0, NX-1
   g_alpha(i+1) = - g_alpha(i)*(alpha-i)/real(i+1)
end do



!boundary conditions
do l = 0,NS ! in these two loops we assign the values of u on the four boundary lines
u(0,l)  = x_bl(t(0))
u(NX,l) = x_bh(t(0))
end do

do i = 0,NX
u(i,0)  = s_bl(t(0))
u(i,NS) = s_bh(t(0))
end do

!initial value assignment
do i = 1, NX-1 ! here all the u-values of the grid except for the boundaries will be assigned
   do l= 1, NS-1
      u(i,l) = f_initial(x(i),s(l))
   end do   
end do

 

call Dateiausgabe(NX,NS,0,t(0),x,s,u) !(Dateiausgabe = data-output) Calling this method will produce an output-data file 
                                      ! with all given values of u at time t   
! Here, the Matrices A and B will be filled and inverted. As long as none
! of the coefficients are time-Dependent, we we only need to do this once and 
! be done with it for the rest of code.
!.................................................................................
! First A will be filled
!do l = 1, NS-1  ! This loop has to be included, in case of an s-dependent diffusion coefficient
   do i = 0, NX   ! fill the coefficient arrays, that will be needed for the A-Matrix
      spi(i) = dt/(dx**alpha)*d_plus(x(i))*fracFac  ! coefficient of the upstream part of the Riesz-derivative
      smi(i) = dt/(dx**alpha)*d_minus(x(i))*fracFac ! coefficient of the downstream part of the Riesz-derivative
      k(i)   = dt*konvektion(x(i))/(2*dx)   !coefficient for the convection terms
      k0(i)  = dt*konvektion(x(i))/(2*dx)
      k1(i)  = dt*konvektionderv(x(i))      !coefficient for the convection-derivative terms
   end do!i
   

   
   !Fill the A matrix
   do j = 1,NX-1 ! loop over the inner points of the A-matrix to fill it with the values resulting from the fractional diffusion
      do i = 0, NX
         if (i == j) then
            A_i(i,j) = 1 - (spi(j) + smi(j))*g_alpha(1)
         else if (i == j+1) then
            A_i(i,j) = - (spi(j)*g_alpha(0) + smi(j)*g_alpha(2))
         else if (i == j - 1) then
            A_i(i,j) = - (spi(j)*g_alpha(2) + smi(j)*g_alpha(0))
         else if (i > j + 1) then
            A_i(i,j) = - smi(j)*g_alpha(i - j + 1)
         else if (i < j - 1) then
            A_i(i,j) = - spi(j)*g_alpha(j - i + 1)
         end if                  
      end do!i
     
      A_i(j-1,j) = A_i(j-1,j) - 0.5*k(j)   ! add values of the convection term
      A_i(j+1,j) = A_i(j+1,j) + 0.5*k(j) 
   end do!j    
   !Here we define the very first and the very last line of the matrix
   A_i(0,0)   = 1.0 
   A_i(NX,NX) = 1.0
   
   do i = 1, NX
      A_i(i,0)    = 0.0
      A_i(i-1,NX) = 0.0
   end do!i
   ! Now that we filled the matrix completely we can calculate the inverse
   call Mat_inv(A_i,A_invers,NX+1)! this is the method to inverse a matrix, we end up with A_invers   
   
   
   

! Fill and inverse all B-Matrices
do i=1,NX  ! We need this loop, because every b-matrix is x-dependent over the term dV/dx.
           ! THerefore we need Nx-diffferent B -Matrices
   do l=1,NS-1  ! Now fill all the inner rows of all the B-Matrices
      do m=0,NS
         if(l==m) then
            B(m,l) = 1
         elseif (m == l+1) then
            B(m,l) = -0.5*dt/2*konvektionderv(x(i))/3/(2*ds)
         elseif (m == l-1) then
            B(m,l) =  0.5*dt/2*konvektionderv(x(i))/3/(2*ds)
         else
          B(m,l) = 0
         endif         
      enddo!m
   enddo!l
   ! Here we fill the very first and the very last row of B
   B(0,0) = 1.
   B(NS,NS) = 1.
   
   do l = 1, NS
      B(l,0)    = 0.0
      B(l-1,NS) = 0.0
   end do!l
   
   
   call Mat_inv(B,B_invers(:,:,i),NS+1)! call the method to invert the matrix again, the output is B_invers(:,:,i)
   
   
end do!i  
! Now we still need to calculate the inverse matrix of i=0 and i=nX, which should just be the identity
do l = 0,NS
   do m =0,NS
   
      if(l==m) then
         B_invers(l,m,0) = 1
         B_invers(l,m,NX) = 1
      else
         B_invers(l,m,0) = 0
         B_invers(l,m,NX) = 0
      end if       
   
   end do!m
end do!l  
!--------------------------------------------------------------------------------------
do ZeitSchritt = 0, NT-1  !(Zeitschritt = time-step) Now we start with the core calculation
! Recipe: First do half a time step in the adiabatic 'direction' , then a full one in the spatial 'direction' and 
! then again half a step in the adiabatic 'direction'

   write(*,*) ZeitSchritt+1 ! write the current timestep to the screen (optional)

! Let us begin with the halfstep in the adiabatic changes
!------------------------------------------------------------------------
! begin to use parallel computing for a faster performance
  !$OMP Parallel DO    
   do i=1,NX-1 ! here we will compute the RightHandSide term of the SAdiabatic-step (therefore called rhs_sad)
      
      do l=1,NS-1
         rhs_sad(l,i) = u(i,l) + 0.5 * (u(i,l+1) - u(i,l-1))*dt/2*konvektionderv(x(i))/3/(2*ds)  ! see manual
      enddo !l
      rhs_sad(0,i) = s_bl(t(ZeitSchritt+1)) !  the left and right side values in s of u are given by the boundary conditions 
      rhs_sad(NS,i) = s_bh(t(ZeitSchritt+1)) 
   enddo !i
  !$OMP End Parallel DO   

  !$OMP Parallel DO    
   do  i=1,NX-1  ! Here we just do the matrix - vector multiplication of inverse-B and rhs_sad

      call Matrix_mal_Vektor(B_invers(:,:,i),rhs_sad(:,i),u(i,:),NS+1) !(Matrix_mal_Vektor = matrix_times_vector) multiplies a matrix and a vector. The result is u(i:)
      

   enddo!i
  !$OMP End Parallel DO
!-------------------------------------------------------------------------   
! We now do a full time step with the spatial terms
! Again, we do parallel computing
   !$OMP Parallel DO 
   do l = 1,NS-1
      rhs(0,l)  = x_bl(t(ZeitSchritt+1)) ! boundary conditions for the left and right boundaries in x for the right hand side term
      rhs(NX,l) = x_bh(t(ZeitSchritt+1))
      do i = 1 , NX-1    
         rhs(i,l) = u(i,l) - 0.5*k0(i)*(u(i+1,l)-u(i-1,l)) + dt*q(x(i),s(l),t(ZeitSchritt+1)) !calculation of the right hand side vector
      end do!i
   enddo   !l
   !$OMP End Parallel DO    
   
   !$OMP Parallel DO 
   do l = 1,NS-1   
      call Matrix_mal_Vektor(A_invers(:,:),rhs(:,l),u(:,l),NX+1) !Call again the matrix-vector method, the result now is u(:,l)

   end do !l 
  !$OMP End Parallel DO  
!-----------------------------------------------------------------------------------
! Do the second half-time step in s-terms, for an overview over the details, see the explanation above 
  !$OMP Parallel DO  
   do i=1,NX-1
      
      do l=1,NS-1
         rhs_sad(l,i) = u(i,l) + 0.5 * (u(i,l+1) - u(i,l-1))*dt/2*konvektionderv(x(i))/3/(2*ds)
      enddo !l
      rhs_sad(0,i) = s_bl(t(ZeitSchritt+1)) 
      rhs_sad(NS,i) = s_bh(t(ZeitSchritt+1)) 
   enddo 
  !$OMP End Parallel DO

  !$OMP Parallel DO  
   do i=1,NX-1 
      call Matrix_mal_Vektor(B_invers(:,:,i),rhs_sad(:,i),u(i,:),NS+1)
      
   enddo!i
  !$OMP End Parallel DO
!-------------------------------------------------------------------------   
   if(modulo(Zeitschritt,500)==0) then  ! make an output file every 500 time steps, the name of the file is the timestep
      call Dateiausgabe(NX,NS,ZeitSchritt+1,t(ZeitSchritt+1),x,s,u) ! Here we are making a data output (Dateiausgabe = Data Output)
   endif   

           
end do! Zeitschritt

end subroutine
!-----------------------------------------------------------------------

Subroutine Dateiausgabe(NX,NS,ZeitSchritt,t,x,s,u) ! This is a method for making a data file
integer :: NX, NS, ZeitSchritt
integer :: i,l

real(kind = 8) :: t
real(kind = 8) ,dimension(0:NX) :: x
real(kind = 8) ,dimension(0:NS) :: s
real(kind = 8) ,dimension(0:NX,0:NS) :: u 
character(99) :: Pfad


write(Pfad,'(A9,I6.6,A4)') 'Daten2.0/', ZeitSchritt, '.dat'! Name of the data file

Open(99, file = Pfad)
write(99,*) '# t = ', t

do  l = 0, NS
   do i = 0, NX
         write(99,*) x(i), s(l), u(i,l) ! Three values are printed: x,s,u
   end do! i
end do   
close(99)


end Subroutine Dateiausgabe


!-----------------------------------------------------------------------
!--------------------InitialConditions----------------------------------
!-----------------------------------------------------------------------
real (kind = 8) function f_initial(x,s)
real(kind = 8) :: x,s,pi

    pi = 3.14159265
    f_initial = 0.0 ! Here the initial values of u is defined, we chose the standard value of 0, but it may be chosen in any way
                    ! the user wants it to be, as long as it is compatible with the boundaries 

end function f_initial

!-----------------------------------------------------------------------
!--------------------DiffusionCoefficients------------------------------
!-----------------------------------------------------------------------
real (kind=8) function d_plus(x) !upstream diffusion coefficient
implicit none
real(kind=8) :: x,s
  d_plus = 1.0 ! this is the diffusion coefficient upstram of the shock

end function d_plus

real (kind = 8) function d_minus(x) !downstream diffusion coefficient
implicit none
real(kind=8) :: x,s
  d_minus = 1.0! this is the diffusion coefficient upstram of the shock

end function d_minus
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

real (kind = 8) function konvektion(x) !convection function 
implicit none
real(kind=8) :: x,t,V1,V2,q
V1 = 1.0
q  = 4.0
V2 = V1/q

   konvektion = (V1+V2)/2.+(V1-V2)/2.*dtanh(-20.*x) ! Here the velocity profile is defined
   

end function konvektion

real (kind = 8) function konvektionderv(x) !convection function derivative
implicit none
real(kind=8) :: x,t,V1,V2,q
V1 = 1.0
q  = 4.0
V2 = V1/q

   konvektionderv = -20.*(V1-V2)/2.*(1.-(dtanh(-20.*x))**2.)! derivative of the velocity profile


end function konvektionderv

real (kind = 8) function q(x,s,t)  ! Sources/Sinks
implicit none
real(kind=8) :: x,s,t,pi
real(kind=8) :: x2,x3,x4
real(kind=8) :: tmx,tmx2,tmx3,tmx4
real(kind = 8) :: term1, term2, term3
  pi = 3.14159265
  q = sqrt(0.1)/sqrt(pi)*exp(-1.*(x/0.1)**2)*sqrt(0.1)/sqrt(pi)*exp(-((exp(s)-1.)/0.1)**2)
  ! Here we define the source of the system, it may be freely chosen by the user
end function q

!-----------------------------------------------------------------------
!--------------------Boundaries-----------------------------------------
!-----------------------------------------------------------------------

real (kind = 8) function x_bl(t) ! lower x-boundary
implicit none
real(kind=8) :: t
   x_bl = 0.0 !lower boundary in x
end function x_bl


real (kind = 8) function x_bh(t) ! upper x-boundary
implicit none
real(kind=8) :: t
    x_bh = 0.0 !
end function x_bh

real (kind = 8) function s_bl(t) ! lower s-boundary
implicit none
real(kind=8) :: t
   !x_bl = sin(t*3.141592654)
   !x_bl = 1.0/sqrt(t)*exp(-1.0/(8.0*t)*(4-3*(t-0.1))**2)
   s_bl = 0.0
end function s_bl


real (kind = 8) function s_bh(t) ! upper s-boundary
implicit none
real(kind=8) :: t
   !x_bh = sin(t*3.141592654)
  ! x_bh = 1.0/sqrt(t)*exp(-1.0/(8.0*t)*(-4-3*(t-0.1))**2)
    s_bh = 0.0
end function s_bh

! From hereon there are the methods for the linear algebra
! These are standard methods, if the user wants he can swap those for 
! custom made ones. He should, however, be aware of the fact, that in this 
! Code the first entry of a matrix denots its colummn and the second its 
! row, in contrast to the usual notation in textbooks. This has
! simply historical reasons, since the code started out as a one dimensional
! code


Subroutine Matrix_mal_Vektor (Matrix, Vektor ,Ergebnis , N) !Matrix times Vector
implicit none
integer :: N
integer,parameter :: dp = 8

real(kind=dp), dimension(1:N, 1:N) :: Matrix
real(kind=dp), dimension(1:N) :: Vektor, Ergebnis
integer :: i,k

   do i = 1,N
      Ergebnis(i) = 0.
      do k = 1, N
         Ergebnis(i) = Ergebnis(i) + Matrix(k,i)*Vektor(k)
      end do !k
   end do !i
End Subroutine Matrix_mal_Vektor
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
Subroutine Matrix_mal_Matrix (Matrix1, Matrix2 ,Ergebnis , N) !Matrix times Matrix
implicit none
integer, parameter :: dp = 8
integer :: N
real(kind = dp), dimension(1:N, 1:N) :: Matrix1, Matrix2, Ergebnis
integer :: i,j,k

do j = 1, N
   do i = 1,N
      Ergebnis(i,j) = 0.
      do k = 1, N
         Ergebnis(i,j) = Ergebnis(i,j) + Matrix1(k,j)*Matrix2(i,k)
      end do !k
   end do !i
end do !j

End Subroutine Matrix_mal_Matrix

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

Subroutine Matrix_plus_Matrix(M1, M2, E, N)!Matrix plus Matrix
implicit none
integer, parameter :: dp = 8
integer :: N
real(kind = dp), dimension(1:N,1:N) :: M1, M2, E
integer :: i, j

do j = 1, N
   do i = 1, N
      E(i,j) = M1(i,j) + M2(i,j)
   end do
end do

End Subroutine Matrix_plus_Matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
Subroutine Gauss(s_Matrix,rhs, N) !Gaussian inversion method for Matrices
implicit none
integer :: N
integer, parameter :: dp = 8
real(kind=dp) , dimension(1:N, 1:N) :: s_Matrix
real(kind=dp) , dimension (1:N+1, 1:N) :: Lsg
real(kind=dp) , dimension (1:n) :: rhs
real(kind=dp) :: norm
real(kind=dp) :: m 
real(kind=dp) :: Maximum
integer :: i,j,k,h

character(9999) :: Zeile, Inhalt1
!-----Matrixgleichung wird erstellt----------------Matrix equation is made ---------
do i = 1, N
   do j = 1, N
      Lsg(i,j) = s_Matrix(i,j)
   end do !j
end do !i

do j = 1,N
   Lsg(N+1,j) = rhs(j)
end do

!Open(99, file = 'Lsg.dat')
!Inhalt1 = ''
!do j = 1,N
!Zeile = ''
!   do i = 1, N+1
!      write(Inhalt1,'(F10.6)') Lsg(i,j)
!      Zeile = trim(Zeile) // ' ' // trim(Inhalt1)
!   end do !i
!   write(99,*) trim(Zeile)
!end do!j
!Close(99)


!Schritt1: System wird äquilibriert ----- equilibration of the system -----
do j = 1,N
   norm = 0.
   do i = 1, N
      norm = norm + abs(Lsg(i,j))
   end do!i
   
   do i = 1,N+1
      Lsg(i,j) = Lsg(i,j)/norm
   end do 
end do!j


!Schritt 2: Jetzt kann mit der Umformung zur oberen Dreiecksmatrix begonnen werden -- transformation to upper triangle matrix begins
do j = 1, N-1
   call Spaltenpivotisierung(Lsg,N,j)
   do i = (j+1),N
      Maximum = 0.
           
      if (Lsg(j,i)/=0) then
         m = Lsg(j,j)/Lsg(j,i)
         do k = 1,N+1
            Lsg(k,i) = m*Lsg(k,i)-Lsg(k,j)
            Lsg(k,i) = Lsg(k,i)/m
         end do!k  

      if (abs(Lsg(j,j)) <= 1.0E-6) then
         write(*,*) 'singulaere Matrix bei',j, Lsg(j,j)
         read(*,*)
         stop
      end if

      end if
   end do!i 
end do!j

!Auf Diagonalelemente eindampfen---- reduce to diagonal elements

do j = N,2, -1
   do i = 1,j-1
      if (Lsg(j,i) /= 0) then
         m = Lsg(j,j)/Lsg(j,i)
         do k = 1,N+1         
            Lsg(k,i) = m*Lsg(k,i) - Lsg(k,j)
            Lsg(k,i) = Lsg(k,i)/m
         end do !k
      end if
   end do !i
end do!j

do j = 1, N
      rhs(j) = Lsg(N+1,j)/Lsg(j,j)
end do!j
End Subroutine Gauss


Subroutine Spaltenpivotisierung(Lsg,N,k)
implicit none
integer :: N
integer :: k
integer,parameter :: dp = 8
real(kind=dp), dimension(1:N+1,1:N) :: Lsg
real(kind=dp) :: Maximum
real(kind=dp) :: h
integer :: i,j
integer :: Max_Zeile

!Finde den betragsgroessten Eintrag in der Spalte------find the largest (aboslute) entry in the column
Maximum = 0.
do j = k,N
   if (abs(Lsg(k,j))>Maximum) then
      Maximum = abs(Lsg(k,j))
      Max_Zeile = j
   end if
end do

!Tausche aktuelle Spalte mit der des betragsgroessten Eintrages----swap the current column with the absolute largest value one

if (Max_Zeile /= k) then
   do i = 1,N+1
         h = Lsg(i,k)
         Lsg(i,k) = Lsg(i,Max_Zeile)
         Lsg(i,Max_Zeile) = h   
   end do!i
end if
End Subroutine Spaltenpivotisierung
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
Subroutine Mat_inv(s_Matrix,s_inv_Matrix, N)  ! inverse Matrix method
implicit none
integer :: N
integer, parameter :: dp = 8
real(kind=dp) , dimension(1:N, 1:N) :: s_Matrix
real(kind=dp) , dimension (1:2*N, 1:N) :: Lsg
real(kind=dp) , dimension (1:N,1:N) :: s_inv_Matrix
real(kind=dp) :: norm
real(kind=dp) :: m 
real(kind=dp) :: Maximum
integer :: i,j,k,h

character(9999) :: Zeile, Inhalt1
!-----Matrixgleichung wird erstellt----------------Matrix equation is made--------
do i = 1, N
   do j = 1, N
      Lsg(i,j) = s_Matrix(i,j)
   end do !j
end do !i

do i = N+1,2*N
   do j = 1,N
      if ((i-N) == j) then         
         Lsg(i,j) = 1.0
      else
         Lsg(i,j) = 0.0
      end if
   end do
end do

!Open(99, file = 'Lsg.dat')
!Inhalt1 = ''
!do j = 1,N
!Zeile = ''
!   do i = 1, 2*N
!      write(Inhalt1,'(F10.6)') Lsg(i,j)
!      Zeile = trim(Zeile) // ' ' // trim(Inhalt1)
!   end do !i
!   write(99,*) trim(Zeile)
!end do!j
!Close(99)


!Schritt1: System wird äquilibriert ------ equilibration of the System ----
do j = 1,N
   norm = 0.
   do i = 1, N
      norm = norm + abs(Lsg(i,j))
   end do!i
   
   do i = 1,2*N
      Lsg(i,j) = Lsg(i,j)/norm
   end do 
end do!j

!Open(99, file = 'Lsg_aeq.dat')
!Inhalt1 = ''
!do j = 1,N
!Zeile = ''
!   do i = 1, 2*N
!      write(Inhalt1,'(F10.6)') Lsg(i,j)
!      Zeile = trim(Zeile) // ' ' // trim(Inhalt1)
!   end do !i
!   write(99,*) trim(Zeile)
!end do!j
!Close(99)

!Schritt 2: Jetzt kann mit der Umformung zur oberen Dreiecksmatrix begonnen werden-----begin transformation for triangle matrix
do j = 1, N-1
   call Spaltenpivotisierung_Matrix(Lsg,N,j)
   do i = (j+1),N   
          
      if (Lsg(j,i)/=0) then
         m = Lsg(j,j)/Lsg(j,i) 

         do k = 1,2*N              
            Lsg(k,i) = m*Lsg(k,i)-Lsg(k,j)
            Lsg(k,i) = Lsg(k,i)/m
         end do!k  

      if (abs(Lsg(j,j)) <= 1.0E-6) then
         write(*,*) 'singulaere Matrix bei',j, Lsg(j,j)
         read(*,*)
         stop
      end if

      end if
   end do!i 
end do!j

!Open(99, file = 'ODM.dat')
!Inhalt1 = ''
!do j = 1,N
!Zeile = ''
!   do i = 1, 2*N
!      write(Inhalt1,'(F10.6)') Lsg(i,j)
!      Zeile = trim(Zeile) // ' ' // trim(Inhalt1)
!   end do !i
!   write(99,*) trim(Zeile)
!end do!j
!Close(99)

!Auf Diagonalelemente eindampfen ----- reduce to diagonal elements

do j = N,2, -1
   do i = 1,j-1
      if (Lsg(j,i) /= 0) then
         m = Lsg(j,j)/Lsg(j,i)
         do k = 1,2*N         
            Lsg(k,i) = m*Lsg(k,i) - Lsg(k,j)
            Lsg(k,i) = Lsg(k,i)/m
         end do !k
      end if
   end do !i
end do!j

do j = 1,N
   do i = N+1,2*N
      s_inv_Matrix(i-N,j) = Lsg(i,j)/Lsg(j,j)
   end do!i
end do!j
End Subroutine Mat_inv


Subroutine Spaltenpivotisierung_Matrix(Lsg,N,k)
implicit none
integer :: N
integer :: k
integer,parameter :: dp = 8
real(kind=dp), dimension(1:2*N,1:N) :: Lsg
real(kind=dp) :: Maximum
real(kind=dp) :: h
integer :: i,j
integer :: Max_Zeile

character(9999) :: Zeile, Inhalt1

!Finde den betragsgroessten Eintrag in der Spalte------find the largest (aboslute) entry in the column
Maximum = 0.
do j = k,N
   if (abs(Lsg(k,j))>Maximum) then
      Maximum = abs(Lsg(k,j))
      Max_Zeile = j
   end if
end do

!Tausche aktuelle Spalte mit der des betragsgroessten Eintrages----swap the current column with the absolute largest value one

if (Max_Zeile /= k) then
   do i = 1,2*N
         h = Lsg(i,k)
         Lsg(i,k) = Lsg(i,Max_Zeile)
         Lsg(i,Max_Zeile) = h   
   end do!i
end if

!Open(99, file = 'Spaltenpivotisierung.dat')
!   Inhalt1 = ''
!   do j = 1,N
!   Zeile = ''
!   do i = 1, 2*N
!      write(Inhalt1,'(F10.6)') Lsg(i,j)
!      Zeile = trim(Zeile) // ' ' // trim(Inhalt1)
!   end do !i
!   write(99,*) trim(Zeile)
!   end do!j
!Close(99)

End Subroutine Spaltenpivotisierung_Matrix


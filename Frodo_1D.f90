Program FRODO

use omp_lib ! import the library for parallel computing

implicit none
integer, parameter :: dp = 8

real(kind = dp), parameter :: alpha = 2.0 ! fractional parameter

real(kind = dp), parameter :: x_min = -40.0 !lower x-boundary of the domain
real(kind = dp), parameter :: x_max = 40.0  !upper x-boundary of the domain


real(kind = dp), parameter :: t_min = 0.0   !starting time of the computation
real(kind = dp), parameter :: t_max = 400.0 !ending time of the computation

integer :: NX,NS, NT  !points in x,s and t

write(*,*) 'x-axis = [',x_min,',',x_max,']'
write(*,*) 'Discretisation of the x-axis:'
!read(*,*) NX
NX = 2000  !see line 19


write(*,*) ' '

write(*,*) 't-axis = [',t_min,',',t_max,']'
write(*,*) 'Discretisation of the t-axis:'
!read(*,*) NT
NT = 100001 !see line 19

call Hauptprogramm(NX,NT,x_min,x_max,t_min,t_max,alpha)  ! Here, the main Programm is called and the relevant parameters are delivered
End Program

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine Hauptprogramm(NX,NT,x_min,x_max,t_min,t_max,alpha)  ! main routine of the programm
implicit none
integer, parameter :: dp = 8
! Here, the main parameters of the system are defined again, for an explanation see the first paragraph of the code
real(kind = dp):: alpha

real(kind = dp):: x_min 
real(kind = dp):: x_max

real(kind = dp):: t_min 
real(kind = dp):: t_max

integer :: NX, NT

!---------- passed on values ------------------------

!functions
real(kind=dp)   :: d_plus
real(kind=dp)   :: d_minus
real(kind=dp)   :: q
real(kind=dp)   :: x_bl, x_bh
real(kind=dp)   :: f_initial
real(kind=dp)   :: konvektion
real(kind=dp)   :: konvektionderv

!Here the necessary arrays for the main computation takes place. Explanations
! for the individual arrays below:

real(kind=dp), dimension (0:NX) :: u ! main array, denoting our solution f

real(kind=dp), dimension (0:NX) :: rhs !A*u = u_rhs
real(kind=dp), dimension (0:NX) :: rhs_pr  

real(kind=dp), dimension (0:NX,0:NX) :: A_i !Matrix for the x advancements

real(kind=dp), dimension (0:NX,0:NX) ::Proxy_NX

real(kind=dp), dimension (0:NX,0:NX) :: A_invers


!x and t arrays
real(kind=dp), dimension (0:NX) :: x
real(kind=dp), dimension (0:NT) :: t

!Grünwald-factors
real(kind=dp), dimension (0:NX) :: g_alpha

!Matrixeinträge
!real(kind=dp), dimension (0:NX) :: spe, sme !spe=sigma_plus_explizit , m = minus
real(kind=dp), dimension (0:NX) :: spi, smi ! i = implizit
real(kind=dp), dimension (0:NX) :: k,k0,k1
real(kind=dp), dimension (0:NX) :: adi


!step-widths
real(kind=dp) :: dx, dy, dt 

!loop variables
integer :: ZeitSchritt
integer :: i,j,l,m

!fractional-factor
real(kind=dp) :: fracFac ! this factor will later be defined to take care of the pre-factors of the 
                         ! Riesz derivative in comparison with the Liouville derivative

!--------------------------------------------------------------


dx = (x_max-x_min)/NX ! assign a value to the step-widths
dt = (t_max-t_min)/NT

fracFac = -0.5/dcos(alpha/2.*3.14159265) ! assign a value to the fractional factor

!assign the value to the x-array
do i = 0, NX
   x(i) = x_min + i*dx
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
u(0)  = x_bl(t(0))
u(NX) = x_bh(t(0))


!initial value assignment
do i = 1, NX-1 ! here all the u-values of the grid except for the boundaries will be assigned
   u(i) = f_initial(x(i)) 
end do

 

call Dateiausgabe(NX,0,t(0),x,u) !(Dateiausgabe = data-output) Calling this method will produce an output-data file 
                                      ! with all given values of u at time t   
! Here, the Matrices A and B will be filled and inverted. As long as none
! of the coefficients are time-Dependent, we we only need to do this once and 
! be done with it for the rest of code.
!.................................................................................
! First A will be filled
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
   
   
  
!--------------------------------------------------------------------------------------
do ZeitSchritt = 0, NT-1  !(Zeitschritt = time-step) Now we start with the core calculation
! Recipe: First do half a time step in the adiabatic 'direction' , then a full one in the spatial 'direction' and 
! then again half a step in the adiabatic 'direction'

   write(*,*) ZeitSchritt+1 ! write the current timestep to the screen (optional)

!-------------------------------------------------------------------------   
! We now do a full time step with the spatial terms
! Again, we do parallel computing
   !$OMP Parallel DO 
      rhs(0)  = x_bl(t(ZeitSchritt+1)) ! boundary conditions for the left and right boundaries in x for the right hand side term
      rhs(NX) = x_bh(t(ZeitSchritt+1))
      do i = 1 , NX-1    
         rhs(i) = u(i) - 0.5*k0(i)*(u(i+1)-u(i-1)) + dt*q(x(i),t(ZeitSchritt+1)) !calculation of the right hand side vector
      end do!i
   !$OMP End Parallel DO    
   
   !$OMP Parallel DO 
      call Matrix_mal_Vektor(A_invers(:,:),rhs(:),u(:),NX+1) !Call again the matrix-vector method, the result now is u(:)

  !$OMP End Parallel DO  
!-----------------------------------------------------------------------------------
!-------------------------------------------------------------------------   
   if(modulo(Zeitschritt,500)==0) then  ! make an output file every 500 time steps, the name of the file is the timestep
      call Dateiausgabe(NX,ZeitSchritt+1,t(ZeitSchritt+1),x,u)! Here we are making a data output (Dateiausgabe = Data Output)
   endif   

           
end do! Zeitschritt

end subroutine
!-----------------------------------------------------------------------

Subroutine Dateiausgabe(NX,ZeitSchritt,t,x,u)! This is a method for making a data file
integer :: NX, NS, ZeitSchritt
integer :: i,l

real(kind = 8) :: t
real(kind = 8) ,dimension(0:NX) :: x
real(kind = 8) ,dimension(0:NX) :: u 
character(99) :: Pfad


write(Pfad,'(A9,I6.6,A4)') 'Daten2.0/', ZeitSchritt, '.dat'! Name of the data file

Open(99, file = Pfad)
write(99,*) '# t = ', t

   do i = 0, NX
         write(99,*) x(i),  u(i) ! Three values are printed: x,s,u
   end do! i
close(99)


end Subroutine Dateiausgabe


!-----------------------------------------------------------------------
!--------------------InitialConditions----------------------------------
!-----------------------------------------------------------------------
real (kind = 8) function f_initial(x)
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
real(kind=8) :: x
  d_plus = 1.0 ! this is the diffusion coefficient upstram of the shock

end function d_plus

real (kind = 8) function d_minus(x) !downstream diffusion coefficient
implicit none
real(kind=8) :: x
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

real (kind = 8) function q(x,t)  ! Sources/Sinks
implicit none
real(kind=8) :: x,t,pi
real(kind=8) :: x2,x3,x4
real(kind=8) :: tmx,tmx2,tmx3,tmx4
real(kind = 8) :: term1, term2, term3
  pi = 3.14159265
  q = sqrt(0.1)/sqrt(pi)*exp(-1.*(x/0.1)**2)
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

Subroutine Matrix_plus_Matrix(M1, M2, E, N)
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

End Subroutine Matrix_plus_Matrix!Matrix plus Matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
Subroutine Gauss(s_Matrix,rhs, N)!Gaussian inversion method for Matrices
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
Subroutine Mat_inv(s_Matrix,s_inv_Matrix, N)
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


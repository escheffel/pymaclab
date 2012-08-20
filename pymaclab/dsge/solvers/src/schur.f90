! Fortran 90 code for solving first and second order approximations using
! Schur method.
! 
! Code written by Paul Gomme, based on Matlab by Paul Klein
!
! For details see Gomme and Klein: Second-order approximations of dynamic
!                                  models without the use of tensors
!
! Requires LAPACK and LAPACK95; both are included in Intel's MKL.
!=============================================================================
module TOOLS
  interface
     function KRON(x,y)
       use PRECISION
       implicit none
       real(long), dimension(:,:), intent(in) :: x, y
       real(long), dimension(size(x,1)*size(y,1), size(x,2)*size(y,2)) :: KRON
     end function KRON
  end interface

  interface
     function EYE(n)
       use PRECISION
       integer, intent(in) :: n
       real(long), dimension(n,n) :: EYE
     end function eye
  end interface

  interface
     function TRANSP(x)
       use PRECISION
       real(long), dimension(:,:), intent(in) :: x
       real(long), dimension(size(x,1), size(x,2)) :: TRANSP
     end function TRANSP
  end interface

  interface
     function TRACEM(x)
       use PRECISION
       real(long), dimension(:,:), intent(in) :: x
       real(long), dimension(size(x,1)/size(x,2)) :: TRACEM
     end function TRACEM
  end interface
  
  interface
     function INV(x)
       use PRECISION
       real(long), dimension(:,:), intent(in) :: x
       real(long), dimension(size(x,1),size(x,2)) :: INV
     end function INV
  end interface

  interface
     subroutine SCHUR_FIRST_ORDER(FUNC, x_vec, P, F, Q, NX, NY, info)
       use PRECISION
       real(long), dimension(:), intent(in)    :: x_vec
       real(long), dimension(:,:), intent(out) :: P, F, Q
       integer, intent(in)                     :: NX, NY
       integer, intent(out)                    :: info

       interface
          function FUNC(x)
            use precision
            real(long), dimension(:), intent(in) :: x
            real(long), dimension(size(x)/2)     :: FUNC
          end function FUNC
       end interface
     end subroutine SCHUR_FIRST_ORDER
  end interface

  interface
     subroutine SCHUR_SECOND_ORDER(FUNC, x_vec, Q_mat, F, P, NS, NC, E_mat, &
          & G_mat, sigma, kx, ky, dev)
       use PRECISION
       implicit none
       real(long), dimension(:), intent(in)    :: x_vec
       real(long), dimension(:,:), intent(in)  :: Q_mat, F, P, sigma
       real(long), dimension(:,:), intent(out) :: E_mat, G_mat
       real(long), dimension(:), intent(out)   :: kx, ky
       integer, intent(in)                     :: NS, NC
       real(long), optional, intent(in)        :: dev

       interface
          function FUNC(x)
            use precision
            real(long), dimension(:), intent(in) :: x
            real(long), dimension(size(x)/2)     :: FUNC
          end function FUNC
       end interface
     end subroutine SCHUR_SECOND_ORDER
  end interface
end module TOOLS
!=============================================================================
pure logical function Z_CHOOSE(s,t)
  integer, parameter :: DPC = KIND((1.0D0,1.0D0))
  complex(DPC), intent(in) :: s, t

  if (abs(s) > abs(t)) then
     z_choose = .true.
  else
     z_choose = .false.
  end if
end function Z_CHOOSE
!=============================================================================
subroutine SCHUR_SOLVE(FUNC, x_vec, F_mat, P_mat, NX, NY, sigma, E_mat, G_mat, k_x, k_y, dev, error_code)
  integer, parameter :: long = kind(0.0d0)
  use TOOLS, only : KRON, EYE, TRANSP, TRACEM, INV, SCHUR_FIRST_ORDER, SCHUR_SECOND_ORDER
  use mkl95_lapack, only : GGES, GESV, TGSYL
  implicit none
  real(long), dimension(:), intent(in)    :: x_vec
  real(long), dimension(:,:), intent(out) :: F_mat, P_mat
  real(long), dimension(:,:), optional, intent(in) :: sigma
  real(long), dimension(:,:), optional, intent(out) :: E_mat, G_mat
  real(long), dimension(:), optional, intent(out)   :: k_x, k_y
  integer, intent(in)                     :: NX, NY
  integer, intent(out), optional          :: error_code
  real(long), optional, intent(in)        :: dev

  interface
     function FUNC(x)
       use precision
       real(long), dimension(:), intent(in) :: x
       real(long), dimension(size(x)/2)     :: FUNC
     end function FUNC
  end interface

  real(long), dimension(:,:), allocatable :: Q_mat
  real(long) :: eps
  integer :: info

  allocate(Q_mat(NX+NY,2*(NX+NY)))

  call SCHUR_FIRST_ORDER(FUNC, x_vec, P_mat, F_mat, Q_mat, NX, NY, info)

  if (present(error_code)) error_code = info
  if (info /= 0) return

  if (present(dev)) then
     eps = dev
  else
     eps = 1d-4
  end if

  if (present(sigma) .and. present(E_mat) .and. present(G_mat) &
       & .and. present(k_x) .and. present(k_y)) then
     call SCHUR_SECOND_ORDER(FUNC, x_vec, Q_mat, F_mat, P_mat, NX, NY, &
          & E_mat, G_mat, sigma, k_x, k_y, eps)
  elseif (present(sigma) .or. present(E_mat) .or. present(G_mat) &
       & .or. present(k_x) .or. present(k_y)) then
     print *, 'Second-order solution requires specifying all of:'
     print *, 'E_mat, G_mat, sigma, k_x and k_y.'
  end if

end subroutine SCHUR_SOLVE
!=============================================================================
subroutine SCHUR_FIRST_ORDER(FUNC, x_vec, P, F, Q, NX, NY, error_code)
  integer, parameter :: long = kind(0.0d0)
  use mkl95_lapack, only : GGES, GESV
  implicit none
  real(long), dimension(:), intent(in) :: x_vec
  real(long), dimension(:,:), intent(out) :: P, F, Q
  integer, intent(in) :: NX, NY
  integer, intent(out) :: error_code

  interface
     function FUNC(x)
       use precision
       real(long), dimension(:), intent(in) :: x
       real(long), dimension(size(x)/2)     :: FUNC
     end function FUNC
  end interface

  interface
     pure logical function Z_CHOOSE(s,t)
       integer, parameter :: DPC = KIND((1.0D0,1.0D0))
       complex(DPC), intent(in) :: s, t
     end function Z_CHOOSE
  end interface

  integer, parameter :: DPC = KIND((1.0D0,1.0D0))

  complex(DPC), dimension(:,:), allocatable :: A_C, B_C, vsl, vsr, &
       & S11, T11, Z11, Z21, work1, work2, work3
  complex(DPC), dimension(:), allocatable :: s_vec, t_vec
  real(long), dimension(:), allocatable   :: x1, x2
  real(long)                              :: dx
  real(long), parameter                   :: eps=1d-5
  integer                                 :: i, N, sdim, info
  
  N = NX+NY

  allocate(x1(2*N), x2(2*N))
  allocate(A_C(N,N), B_C(N,N), vsl(N,N), vsr(N,N), s_vec(N), t_vec(N), &
       & S11(NX,NX), T11(NX,NX), Z11(NX,NX), Z21(NY,NX))

  do i = 1, 2*N
     x1 = x_vec
     x2 = x_vec
     x1(i) = x1(i)*(1d0+eps)+eps
     x2(i) = x2(i)*(1d0-eps)-eps
     dx = x1(i)-x2(i)
     Q(:,i) = (FUNC(x1)-FUNC(x2))/dx
  end do
  B_C = -Q(:,N+1:2*N)
  A_C = Q(:,1:N)

  error_code = 0

  call GGES(A_C, B_C, s_vec, t_vec, vsl, vsr, z_choose, sdim, info)  

  if (info == N+2) then
     print *, 'LAPACK: Problem with precision after reordering.'
  else if (info /= 0) then
     print *, 'LAPACK: Problem in Schur decomposition.'
     error_code = 1
     return
  end if

  if (sdim /= NX) then
     write(6,10) NX, sdim, info
     error_code = 2
     return
  end if

10 format(' Economic problem in Schur decomposition' / &
        & ' Number of states:', t40, i5 / &
        & ' Number of correctly ordered elements:', t40, i5 / &
        & ' Value of info from GGES:', t40, i5)

20 format(2(1x,e20.10e3))

  Z11 = vsr(1:NX,1:NX)
  S11 = A_C(1:NX,1:NX)
  T11 = B_C(1:NX,1:NX)
  Z21 = vsr(NX+1:N,1:NX)
  
  deallocate(A_C, B_C, vsr, vsl, s_vec, t_vec)
  allocate(work1(NX,NX), work2(NX,NX), work3(NX,NX))

  work1 = S11
  work2 = T11
! compute inv(S11)*T11 -- stored in work2
  call GESV(work1, work2, info=info)

! compute inv(Z11) -- stored in work3
  work1 = Z11
  work3 = 0d0
  do i=1,NX
     work3(i,i) = 1d0
  end do
  call GESV(work1, work3, info=info)
  
  F = DBLE(matmul(Z21,work3))

  P = DBLE(matmul(matmul(Z11,work2),work3))
  
  deallocate(S11, T11, Z11, Z21)


  if (info /= 0) then
     print *, 'Schur solver: failed.'
  end if
end subroutine SCHUR_FIRST_ORDER
!=============================================================================
! Subroutine: SCHUR_SECOND_ORDER
!
! Purpose: Solves for the second-order accurate approximation 
!          of the solution to a dynamic model
!
! The model is E_t[f(x_{t+1},y_{t+1},x_t,y_t)]=0
!          
! For the definition of the gradient and the Hessian, see
! Magnus and Neudecker: Matrix Differential Calculus with 
! Applications in Statistics and Econometrics.
!
! Outputs: 
!   E_mat, G_mat ... matrices of second-order terms
!   kx, ky ......... vectors adjusting constant term
!
! Inputs:
!   FUNC ...... which evaluates the Euler equations and constraints
!   x_vec ..... point around which approximation is computed; consists of:
!               future state variables, future control variables,
!               current state variables, current control variables
!               (in that order, in logarithms)
!   Q_mat ... matrix of first-order partial derivatives
!   F, P .... matrices characterizing the first-order solutions
!        (Q_mat, F and P are computed by first calling SCHUR_FIRST_ORDER)
!   NX ...... number of state variables
!   NY ...... number of control variables
!   sigma ... variance-covariance matrix for the innovations to the
!             state variables 
!   dev ..... (optional) the "deviations" used to compute the Hessian matrix
!
! x_{t+1} = kx + Px_t + (1/2)*kron(eye(nx),x'_t)Gx_t + eps_{t+1}
!
! y_{t}   = ky + Fx_t + (1/2)*kron(eye(ny),x'_t)Ex_t
!
! Calls: KRON, EYE, TRANSP, TRACEM, INV, GGES, TGSYL
!
! Original code written by Paul Gomme, based on Matlab by Paul Klein
! Later modifications by Paul Gomme
!
! For details see Gomme and Klein: Second-order approximations of dynamic
!                                  models without the use of tensors
! 
subroutine SCHUR_SECOND_ORDER(FUNC, x_vec, Q_mat, F, P, NX, NY, E_mat, &
     & G_mat, sigma, kx, ky, dev)
  integer, parameter :: long = kind(0.0d0)
  use TOOLS, only : KRON, EYE, TRANSP, TRACEM, INV
  use mkl95_lapack, only : GGES, TGSYL
  implicit none
  integer, parameter :: DPC = KIND((1.0D0,1.0D0))
  real(long), dimension(:), intent(in)    :: x_vec
  real(long), dimension(:,:), intent(in)  :: Q_mat, F, P, sigma
  real(long), dimension(:,:), intent(out) :: E_mat, G_mat
  real(long), dimension(:), intent(out)   :: kx, ky
  integer, intent(in)                     :: NX, NY
  real(long), optional, intent(in)        :: dev

  interface
     function FUNC(x)
       use precision
       real(long), dimension(:), intent(in) :: x
       real(long), dimension(size(x)/2)     :: FUNC
     end function FUNC
  end interface

  real(long), dimension(:), allocatable   :: x_prime, feval
  real(long), dimension(:,:), allocatable :: f1, f2, f3, f4
  real(long), dimension(:,:), allocatable :: A_1, B_1, B_2, B_4, C_1, C_2
  real(long), dimension(:,:), allocatable :: Inx, Iny, Im, M_mat
  real(long), dimension(:,:), allocatable :: A_tilde, B_tilde, C_tilde
  real(long), dimension(:,:), allocatable :: D_tilde, E_tilde, F_tilde
  real(long), dimension(:,:), allocatable :: P_1, Q_1, U_1, V_1
  real(long), dimension(:), allocatable   :: alpha_r, alpha_i, beta
  real(long), dimension(:,:), allocatable :: Hess, ma, eyeff
  real(long), dimension(:), allocatable   :: y, kxy, ve
  real(long)                              :: deviation
  integer                                 :: i, j, k, N, M, info

  M = NX+NY
  N = 2*M

  allocate(x_prime(N), feval(M))
  allocate(f1(M,NX), f2(M,NY), f3(M,NX), f4(M,NY))
  allocate(A_1(M*NX,NX), B_1(M*NX,NX*NX))
  allocate(B_2(M*NX,NX*NY), B_4(M*NX,NX*NY))
  ! 2011-04-13: Corrected dimension of C_2
  allocate(C_1(NX*NY,NX*NY), C_2(NX*NY,NX*NX), M_mat(N,NX))
  allocate(A_tilde((NX+NY)*NX,(NX+NY)*NX), D_tilde((NX+NY)*NX,(NX+NY)*NX))
  allocate(B_tilde(NX,NX), E_tilde(NX,NX))
  allocate(C_tilde((NX+NY)*NX,NX), F_tilde((NX+NY)*NX,NX))
  ! 2011-04-13: Added declaration of Im
  allocate(Inx(NX,NX), Iny(NY,NY), Im(M,M))
  allocate(Hess(M*N,N), y(M))
  allocate(ma(M,M), ve(M), kxy(M))

  if (present(dev)) then
     deviation = dev
  else
     deviation = 1d-4
  end if

  feval = FUNC(x_vec)

  do i = 1, N
     x_prime = x_vec
     x_prime(i) = x_prime(i) + deviation
     y = FUNC(x_prime)

     x_prime = x_vec
     x_prime(i) = x_prime(i) - deviation
     y = (y + FUNC(x_prime) - 2d0*feval) / (deviation**2)

     do k = 1, M
        Hess(N*(k-1)+i,i) = y(k)
     end do

     do j = i+1, N
        x_prime = x_vec
        x_prime(i) = x_prime(i) + deviation
        x_prime(j) = x_prime(j) + deviation
        y = FUNC(x_prime)

        x_prime = x_vec
        x_prime(i) = x_prime(i) - deviation
        x_prime(j) = x_prime(j) - deviation
        y = y + FUNC(x_prime)

        x_prime = x_vec
        x_prime(i) = x_prime(i) + deviation
        x_prime(j) = x_prime(j) - deviation
        y = y - FUNC(x_prime)

        x_prime = x_vec
        x_prime(i) = x_prime(i) - deviation
        x_prime(j) = x_prime(j) + deviation
        y = (y - FUNC(x_prime)) / (4d0*deviation**2d0)
        do k = 1, M
           Hess(N*(k-1)+j,i) = y(k)
           Hess(N*(k-1)+i,j) = y(k)
        end do
     end do
  end do

  f1 = Q_mat(:,1:NX)
  f2 = Q_mat(:,NX+1:M)
  f3 = Q_mat(:,M+1:M+NX)
  f4 = Q_mat(:,M+NX+1:2*M)

  Inx = EYE(NX)
  Iny = EYE(NY)
  Im = EYE(M)

  M_mat(1:NX,:) = P
  M_mat(NX+1:M,:) = matmul(F,P)
  M_mat(M+1:M+NX,:) = Inx
  M_mat(M+NX+1:2*M,:) = F

  A_1 = MATMUL(KRON(Im,transpose(M_mat)), MATMUL(Hess,M_mat))

  B_1 = KRON(f1,Inx)
  B_2 = KRON(f2,Inx)
  B_4 = KRON(f4,Inx)

  C_1 = KRON(Iny,transpose(P))
  C_2 = KRON(F,Inx)

  ! Code to solve Sylvester equation (with lots of setup). 2010-04-20

  A_tilde(:,1:NX*NY) = B_4
  A_tilde(:,NX*NY+1:NX*(NX+NY)) = B_1 + matmul(B_2,C_2)
  B_tilde = -P
  C_tilde = -A_1
  D_tilde = 0d0
  D_tilde(:,1:NX*NY) = matmul(B_2,C_1)
  E_tilde = Inx
  F_tilde = 0d0

  ! Need to put (A_tilde,D_tilde) into generalized Schur form
  allocate(alpha_r((NX+NY)*NX), alpha_i((NX+NY)*NX), beta((NX+NY)*NX))
  allocate(P_1(M*NX,M*NX), Q_1(M*NX,M*NX))

  call GGES(A_tilde, D_tilde, alpha_r, alpha_i, beta, P_1, Q_1, info=info)
     
  deallocate(alpha_r,alpha_i,beta)

  ! Need to put (B_tilde,E_tilde) into generalized Schur form
  allocate(alpha_r(NX), alpha_i(NX), beta(NX))
  allocate(U_1(NX,NX), V_1(NX,NX))

  call GGES(B_tilde, E_tilde, alpha_r, alpha_i, beta, U_1, V_1, info=info)

  deallocate(alpha_r,alpha_i,beta)

  C_tilde = matmul(transpose(P_1), matmul(C_tilde,V_1))
  F_tilde = matmul(transpose(P_1), matmul(F_tilde,V_1))

  ! Solve the generalized Sylvester equation
  call TGSYL(A_tilde, B_tilde, C_tilde, D_tilde, E_tilde, F_tilde, info=info)

  C_tilde = matmul(matmul(Q_1,C_tilde), transpose(V_1))

  E_mat = C_tilde(1:NX*NY,:)
  G_mat = C_tilde(NX*NY+1:NX*NX*NY,:)

  ma(:,1:NX) = f1 + matmul(f2,F)
  ma(:,NX+1:M) = f2 + f4
  ma = 2d0*ma

  allocate(eyeff(2*M,NX)) ! 2009-09-21

  eyeff = 0d0 ! 2009-09-21
  eyeff(1:NX,:) = Inx ! 2009-09-21
  do i = 1, NY
     eyeff(NX+i,:) = F(i,:) ! 2009-09-21
  end do

! 2009-09-21
  ve = matmul(f2, TRACEM(matmul(KRON(Iny,sigma), E_mat))) &
       & + TRACEM(matmul(matmul(KRON(Im,transpose(eyeff)),Hess), &
       & matmul(eyeff,sigma)))

  kxy = -matmul(INV(ma), ve)
  kx = kxy(1:NX)
  ky = kxy(NX+1:NX+NY)

  deallocate(x_prime, feval)
  deallocate(f1, f2, f3, f4, A_1, B_1, B_2, B_4, C_1, C_2, M_mat)
  deallocate(A_tilde, B_tilde, C_tilde, D_tilde, E_tilde, F_tilde)
  deallocate(Inx, Iny)
  deallocate(Hess, y, ma, ve, kxy)
end subroutine SCHUR_SECOND_ORDER
!=============================================================================
function TRACEM(x)
  integer, parameter :: long = kind(0.0d0)
  use TOOLS, only : TRACE
  implicit none
  real(long), dimension(:,:), intent(in) :: x
  real(long), dimension(size(x,1)/size(x,2)) :: TRACEM
  integer :: n, m, i

  n = size(x,2)
  m = size(x,1)/n

  do i = 1, m
     TRACEM(i) = TRACE(x((n*(i-1)+1):i*N,:))
  end do
end function TRACEM
!=============================================================================
function TRACE(x)
  integer, parameter :: long = kind(0.0d0)
  implicit none
  real(long), dimension(:,:), intent(in) :: x
  real(long) :: TRACE
  integer :: n, i

  n = size(x,1)

  TRACE = 0d0
  do i = 1, n
     TRACE = TRACE + x(i,i)
  end do
end function TRACE
!=============================================================================
function TRANSP(x)
  integer, parameter :: long = kind(0.0d0)
  implicit none
  real(long), dimension(:,:), intent(in) :: x
  real(long), dimension(size(x,1), size(x,2)) :: TRANSP
  integer :: i, n, m

  n = size(x,2)
  m = size(x,1)/n

  do i = 1, m
     TRANSP((n*(i-1)+1):i*n,1:n) = transpose(x((n*(i-1)+1):i*n,1:n))
  end do
end function TRANSP
!=============================================================================
function DET(q)
  integer, parameter :: long = kind(0.0d0)
  use mkl95_lapack, only : GETRF
  implicit none
  real(long), dimension(:,:), intent(in) :: q
!  real(long), dimension(1,1) :: DET
  real(long) :: DET
  real(long), dimension(:,:), allocatable :: a
  integer :: N, info, i
  integer, dimension(:), allocatable :: ipiv
  
  N = size(q,1)
  allocate(a(N,N), ipiv(N))
  a(:,:) = q(:,:)
  call GETRF(a, IPIV=ipiv, INFO=info)
  if (info .ne. 0) then
     print *, 'DET: problem in return from DGETRF.'
     DET = 0d0
  else
     DET = 1d0
     do i = 1, N
        DET = DET*a(i,i)
        if (i /= ipiv(i)) DET = -DET
     end do
  end if
end function DET
!=============================================================================
function EYE(n)
  integer, parameter :: long = kind(0.0d0)
  implicit none
  integer, intent(in) :: n
  real(long), dimension(n,n) :: EYE
  integer :: i

  EYE(:,:) = 0d0
  do i = 1, n
     EYE(i,i) = 1d0
  end do
end function EYE
!=============================================================================
function KRON(x,y)
  integer, parameter :: long = kind(0.0d0)
  implicit none
  real(long), dimension(:,:), intent(in) :: x, y
  real(long), dimension(size(x,1)*size(y,1), size(x,2)*size(y,2)) :: KRON
  integer :: NX1, NX2, NY1, NY2, ix1, ix2, iy1, iy2, iz1, iz2

  NX1 = size(x,1)
  NX2 = size(x,2)
  NY1 = size(y,1)
  NY2 = size(y,2)

  do ix1 = 1, NX1
     do ix2 = 1, NX2
        do iy1 = 1, NY1
           do iy2 = 1, NY2
              iz1 = (ix1-1)*NY1 + iy1
              iz2 = (ix2-1)*NY2 + iy2
              KRON(iz1,iz2) = x(ix1,ix2)*y(iy1,iy2)
           end do
        end do
     end do
  end do
end function KRON
!=============================================================================
function INV(x)
  integer, parameter :: long = kind(0.0d0)
  use mkl95_lapack, only : GESV
  implicit none
  real(long), dimension(:,:), intent(in) :: x
  real(long), dimension(size(x,1),size(x,2)) :: INV
  integer, dimension(:), allocatable :: indx
  real(long), dimension(:,:), allocatable :: y
  integer :: info, i, N

  if (size(x,1) .ne. size(x,2)) then
     print *, 'INV: Matrix must be square.'
     return
  end if

  N = size(x,1)
  allocate(indx(N), y(N,N))

  INV(:,:) = x(:,:)

  y(:,:) = 0d0
  do i = 1, N
     y(i,i) = 1d0
  end do

  call GESV(INV, y, INFO=info)

  if (info .ne. 0) then
     print *, 'INV: problem in call to DGESV'
  end if

  INV = y
end function INV

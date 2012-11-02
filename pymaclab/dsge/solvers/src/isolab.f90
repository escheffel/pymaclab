subroutine isolab(g,q,x,y,alpha,beta,m,l,retcon)

  use solab
  implicit none
  integer, intent(in) :: m,l
  integer, intent(out) :: retcon
  real(kind=8), intent(in)  :: x(m,m), y(m,m)
  real(kind=8), intent(out) :: g(m-l,l), q(l,l)
  complex(kind=8), intent(out) :: alpha(m)
  complex(kind=8), intent(out) :: beta(m)

  call zsolab(g,q,x,y,alpha,beta,m,l,retcon)

end subroutine isolab

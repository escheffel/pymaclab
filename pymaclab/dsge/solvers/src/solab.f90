module solab

contains

    subroutine zsolab(f,p,a,b,n,k,retco)

    implicit none
    integer, intent(in) :: n,k
    integer, intent(out) :: retco
!    integer i,j,l,m
!   the above aren't used
    integer lda,ldb,sdim,ldvsl,ldvsr,lwork,info,ipiv(k)
    logical bwork(n)
    real(kind=8), intent(in)  :: a(n,n), b(n,n)
    real(kind=8), intent(out) :: f(n-k,k), p(k,k)
    real(kind=8) rwork(8*n)
!    complex*16 s(n,n), t(n,n)
    complex(kind=8) s(n,n), t(n,n)
!    complex*16 alpha(n),beta(n),vsl(n,n),vsr(n,n)
    complex(kind=8) alpha(n),beta(n),vsl(n,n),vsr(n,n)
!    complex*16 tempwork(1)
    complex(kind=8) tempwork(1)
!    complex*16, allocatable :: work(:)
    complex(kind=8), allocatable :: work(:)
!    complex*16 z11work(8*k)
    complex(kind=8) z11work(8*k)
!    complex*16 z11(k,k), z11i(k,k)
    complex(kind=8) z11(k,k), z11i(k,k)
!    complex*16 s11i(k,k), t11(k,k), pp(k,k)
    complex(kind=8) s11i(k,k), t11(k,k), pp(k,k)
!    complex*16 z21(n-k,k), eyek(k,k)
!    complex(kind=8) z21(n-k,k), eyek(k,k)
    complex(kind=8) z21(n-k,k)

    retco = 0
    lwork = -1
    s = dcmplx(a)
    t = dcmplx(b)

    ldvsl = n
    ldvsr = n
    lda = n
    ldb = n

    call zgges('N','V','S',funcg,n,s,lda,t,ldb,sdim,alpha,beta,vsl,ldvsl,vsr,ldvsr,tempwork,lwork,rwork,bwork,info)
    lwork = int(tempwork(1))
    allocate(work(lwork))

    call zgges('N','V','S',funcg,n,s,lda,t,ldb,sdim,alpha,beta,vsl,ldvsl,vsr,ldvsr,work,lwork,rwork,bwork,info)


    if ((abs(alpha(k))<abs(beta(k))) .or. (abs(alpha(k+1))>abs(beta(k+1)))) then
     print*, 'Error: wrong number of stable eigenvalues'
     retco = 1
     return
    end if

    z11  = vsr(1:k,1:k)
    z11i = vsr(1:k,1:k)
    z21  = vsr(k+1:n,1:k)
    s11i = s(1:k,1:k)
    t11  = t(1:k,1:k)
    call zgetrf(k,k,z11i,k,ipiv,info)
    call zgetri(k,z11i,k,ipiv,z11work,8*k,info)
    call ztrtri('U','N',k,s11i,k,info)

    f  = real(matmul(z21,z11i),kind=8)
    pp = matmul(s11i,t11)
    pp = matmul(z11,pp)
    pp = matmul(pp,z11i)
    p  = real(pp,kind=8)

    end subroutine


    logical function funcg(sii,tii)
!     complex*16 sii,tii
     complex(kind=8) sii,tii
     funcg = abs(sii)>abs(tii)
    end function

end module solab

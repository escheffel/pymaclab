subroutine bkfilter(n, y, up, dn, k, ybp, yt)

! Aubhik 20.02.2005
!
! Baxter and King (1999) Band-Pass filter
!
! y is a data vector of length n
! ybp is the band pass filtered series, with 0.0 for the first and last k observations
! n is the length of the data vector y
! up is the lower frequency
!	for example, the business cycle band will have up = 6 for quarterly data
!	and up = 1.5 for annual data
! dn is the higher frequency
!	for example, the business cycle band will have dn = 32 for quarterly data
!	and dn = 8 for annual data
! k is the number of leads and lags used by the symmetric moving-average 
!	representation of the band pass filter.  Baxter and King (1999) find
!	that k = 12 is a practical choice.
! 
! Algorithm provided by Bob King.  This subroutine is the Fortran 90 version of 
! bpf.m and filtk.m written by Baxter and King (1999).  The band pass filter is 
! a moving average involving 2k+1 terms.  The weights are symmetric and calculated
! in akvec.  

implicit none

integer n, up, dn, k
real y(n), ybp(n), yt(n)

integer j
real pi, omlbar, omubar, kr, kj, theta
real akvec(k+1), avec(2*k+1), tvec(2*k+1)

intent(in):: n, up, dn, k, y
intent(out):: ybp, yt

pi = 2.0*acos(0.0)

omlbar = 2.0*pi/dble(dn)
omubar = 2.0*pi/dble(up)

akvec = 0.0
avec = 0.0

kr = dble(k)

akvec(1) = (omubar - omlbar)/pi

do j = 1, k
    kj = dble(j)
    akvec(j+1) = (dsin(dble(kj*omubar)) - dsin(dble(kj*omlbar)))/(kj*pi)
end do

theta = akvec(1) + 2.0*sum(akvec(2:k+1))
theta = -1.0*(theta/(2.0*kr + 1.0))

akvec = akvec + theta

avec(k+1) = akvec(1)

do j = 1, k
    avec(k+1 - j) = akvec(j+1)
    avec(k+1 + j) = akvec(j+1)
end do

do j = 1, k*2+1
    tvec(j) = 1.0 - avec(j)
end do

ybp = 0.0
yt = 0.0

do j = k + 1, n - k
    ybp(j) = dot_product(avec, y(j - k: j + k))
    yt(j) = dot_product(tvec, y(j - k: j + k))
end do

end subroutine bkfilter

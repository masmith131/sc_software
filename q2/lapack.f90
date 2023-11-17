! This file contains some explicit interfaces to LAPACK routines.
! Not every lapack routine has an interface here. Feel free to add some yourself.
module lapack

   interface

      subroutine sgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
         real alpha, beta
         integer k, lda, ldb, ldc, m, n
         character transa, transb
         real a( lda, * ), b( ldb, * ), c( ldc, * )
      end subroutine

      subroutine dgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
         double precision alpha, beta
         integer k, lda, ldb, ldc, m, n
         character transa, transb
         double precision a( lda, * ), b( ldb, * ), c( ldc, * )
      end subroutine

      subroutine zgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
         complex*16 alpha, beta
         integer k, lda, ldb, ldc, m, n
         character transa, transb
         complex*16 a( lda, * ), b( ldb, * ), c( ldc, * )
      end subroutine

      subroutine cgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
         complex alpha, beta
         integer k, lda, ldb, ldc, m, n
         character transa, transb
         complex a( lda, * ), b( ldb, * ), c( ldc, * )
      end subroutine

      subroutine dlacpy( uplo, m, n, a, lda, b, ldb )
         character          uplo
         integer            lda, ldb, m, n
         double precision   a( lda, * ), b( ldb, * )
      end subroutine

      subroutine slacpy( uplo, m, n, a, lda, b, ldb )
         character          uplo
         integer            lda, ldb, m, n
         real               a( lda, * ), b( ldb, * )
      end subroutine

      subroutine zlacpy( uplo, m, n, a, lda, b, ldb )
         character          uplo
         integer            lda, ldb, m, n
         complex*16         a( lda, * ), b( ldb, * )
      end subroutine

      subroutine clacpy( uplo, m, n, a, lda, b, ldb )
         character          uplo
         integer            lda, ldb, m, n
         complex            a( lda, * ), b( ldb, * )
      end subroutine

      subroutine slaset( uplo, m, n, alpha, beta, a, lda )
         character          uplo
         integer            lda, m, n
         real               alpha, beta
         real               a( lda, * )
      end subroutine

      subroutine dlaset( uplo, m, n, alpha, beta, a, lda )
         character          uplo
         integer            lda, m, n
         double precision   alpha, beta
         double precision   a( lda, * )
      end subroutine

      subroutine zlaset( uplo, m, n, alpha, beta, a, lda )
         character          uplo
         integer            lda, m, n
         complex*16         alpha, beta
         complex*16         a( lda, * )
      end subroutine

      subroutine claset( uplo, m, n, alpha, beta, a, lda )
         character          uplo
         integer            lda, m, n
         complex            alpha, beta
         complex            a( lda, * )
      end subroutine

      subroutine slaswp( n, a, lda, k1, k2, ipiv, incx )
         integer            incx, k1, k2, lda, n
         integer            ipiv( * )
         real               a( lda, * )
      end subroutine

      subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )
         integer            incx, k1, k2, lda, n
         integer            ipiv( * )
         double precision   a( lda, * )
      end subroutine

      subroutine zlaswp( n, a, lda, k1, k2, ipiv, incx )
         integer            incx, k1, k2, lda, n
         integer            ipiv( * )
         complex*16         a( lda, * )
      end subroutine

      subroutine claswp( n, a, lda, k1, k2, ipiv, incx )
         integer            incx, k1, k2, lda, n
         integer            ipiv( * )
         complex            a( lda, * )
      end subroutine

      subroutine sgetrf( m, n, a, lda, ipiv, info )
         integer            info, lda, m, n
         integer            ipiv( * )
         real               a( lda, * )
      end subroutine

      subroutine dgetrf( m, n, a, lda, ipiv, info )
         integer            info, lda, m, n
         integer            ipiv( * )
         double precision   a( lda, * )
      end subroutine

      subroutine zgetrf( m, n, a, lda, ipiv, info )
         integer            info, lda, m, n
         integer            ipiv( * )
         complex*16         a( lda, * )
      end subroutine

      subroutine cgetrf( m, n, a, lda, ipiv, info )
         integer            info, lda, m, n
         integer            ipiv( * )
         complex            a( lda, * )
      end subroutine

      subroutine dgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info )
         implicit none
         character          jobvl, jobvr
         integer            info, lda, ldvl, ldvr, lwork, n
         double precision   a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), wi( * ), work( * ), wr( * )
      end subroutine

      subroutine sgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info )
         implicit none
         character          jobvl, jobvr
         integer            info, lda, ldvl, ldvr, lwork, n
         real   a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), wi( * ), work( * ), wr( * )
      end subroutine

      subroutine zgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
         implicit none
         character          jobvl, jobvr
         integer            info, lda, ldvl, ldvr, lwork, n
         double precision   rwork( * )
         complex*16         a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), w( * ), work( * )
      end subroutine

      subroutine cgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
         implicit none
         character          jobvl, jobvr
         integer            info, lda, ldvl, ldvr, lwork, n
         real   rwork( * )
         complex         a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), w( * ), work( * )
      end subroutine

      subroutine sgetri( n, a, lda, ipiv, work, lwork, info )
         integer            info, lda, lwork, n
         integer            ipiv( * )
         real               a( lda, * ), work( * )
      end subroutine

      subroutine dgetri( n, a, lda, ipiv, work, lwork, info )
         integer            info, lda, lwork, n
         integer            ipiv( * )
         double precision   a( lda, * ), work( * )
      end subroutine

      subroutine zgetri( n, a, lda, ipiv, work, lwork, info )
         integer            info, lda, lwork, n
         integer            ipiv( * )
         complex*16         a( lda, * ), work( * )
      end subroutine

      subroutine cgetri( n, a, lda, ipiv, work, lwork, info )
         integer            info, lda, lwork, n
         integer            ipiv( * )
         complex            a( lda, * ), work( * )
      end subroutine

   end interface

end module

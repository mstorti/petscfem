      subroutine   dsyev2(n, a, w, work, lwork, info );
      implicit real*8 (a-h,o-z)
      call dsyev('V','U',n, a,lda, w, work, lwork, info);
      return
      end

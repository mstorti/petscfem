      subroutine   dsyev2(n, a, w, work, lwork, info );
      implicit real*8 (a-h,o-z)
      write (*,*) n,lwork
      call dsyev('V','U',n, a,n, w, work, lwork, info);
      call qq('V')
      return
      end

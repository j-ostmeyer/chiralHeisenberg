      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
      common/chi/chi(kvol,3),seed
      real*8 chidp(kvol,3)
      real*8 seed
c
      call sread_sp
      do ia=1,3
      do i=1,kvol
      chidp(i,ia)=chi(i,ia)
      write(6,*) ia, i, chidp(i,ia)
      enddo
      enddo
      call swrite_dp(chidp)
      stop
      end
c*****************************************************************
      subroutine sread_sp
      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
      common/chi/chi(kvol,3),seed
      real*8 seed
      open(unit=10,file='con',
     1     status='unknown',form='unformatted')
      read (10) chi,seed
      close(10)
      return
      end
c
      subroutine swrite_dp(chidp)
      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
      common/chi/chi(kvol,3),seed
      real*8 chidp(kvol,3)
      real*8 seed
      open(unit=31,file='condp',
     1     status='unknown',form='unformatted')
      write (31) chidp,seed
      close(31)
      return
      end


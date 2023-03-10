       program dwf3d
c*******************************************************************
c    Rational Hybrid Monte Carlo algorithm for 3d U(1) chiral Heisenberg model 
c    with N=1 Domain Wall Fermion
c
c    the "third" direction is actually indexed 4 in the code - sorry!
c      (should have used C I know)
c
c           { 1 - hermition mass term psibar psi
c    imass= { 3 - antih.    mass term i psibar gamma_3 psi
c     linear combinations of above require code modification
c
c    code exploits fact that gamma_3 is diagonal in Dirac basis used - speeds up 
c    evolution slightly, and measurement by factor of two. 
c
c    Fermion expectation values are measured using a noisy estimator.
c    on the Wilson matrix, which has dimension 2*4*kvol*kthird
c    inversions done using congrad, and matrix multiplies with dslash, dslashd
c
c    trajectory length is random with mean dt*iterl
c    The code runs for a fixed number iter2 of trajectories.
c
c    Phi: pseudofermion field 
c    chi: auxiliary scalar field with 3 real components living on domain walls
c    am: bare fermion mass  - not currently implemented
c    actiona: running average of total action
c
c                                               SJH January  2023
c*******************************************************************
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksizet*ksize*ksize)
      parameter(respbp=0.000001,rescgg=0.000001,rescga=0.000000001)
      parameter(itermax=1000)
      parameter(ndiag=25,ndiagg=12)
      parameter(One=1.0)
      common/chi/chi(kvol,3),seed
      common/trial/chit(kvol,3),pp(kvol,3)
      common/sforce/dSdchi(kvol,3)
      common /para/beta,am3,ibound
      common/remez2/anum2(0:ndiag),aden2(ndiag),
     &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4/anum4(0:ndiag),aden4(ndiag),
     &              bnum4(0:ndiag),bden4(ndiag)
      common/remez2g/anum2g(0:ndiagg),aden2g(ndiagg),
     &              bnum2g(0:ndiagg),bden2g(ndiagg)
      common /neighb/id(kvol,3),iu(kvol,3)
      common/param/ancg,ancgh,ancgf,ancgpf
      common/trans/tpi 
      common/dum1/ R(kthird,kvol,4,2),ps(kvol,2)
      common/vector/X1(kthird,kvol,4,2)
      common/ranseed/yran,idum
      common/v/v(97)
      complex Phi(kthird,kvol,4,2),X0(kthird,kvol,4,2)
      complex R,zi,qq,qbqb
      complex X1
      complex a,b
      real rano
      real*8 H0,H1,S0,S1,dH,dS,hg,hp
      real*8 seed
      real*8 anum2,aden2,bnum2,bden2
      real*8 anum4,aden4,bnum4,bden4
      real*8 anum2g,aden2g,bnum2g,bden2g
      real cond(3),susc(3),conda(3)
c*******************************************************************
c     input
c*******************************************************************
      ibound=-1
      istart=-1
      iread=1
      iwrite=1
      iprint=5
      iseed=0
      icheck=100
      zi=(0.0,1.0)
      tpi=2.0*acos(-1.0)
c*******************************************************************
c     end of input
c*******************************************************************
      open(unit=7,file='output',status='unknown')
      open(unit=25,file='midout',status='unknown')
      open(unit=98,file='control',status='unknown')
      open(unit=36,file='remez2',status='unknown')
      open(unit=37,file='remez4',status='unknown')
      open(unit=38,file='remez2g',status='unknown')
      if(iread.eq.1) then
      call sread
      endif
      read(25,*) dt,beta,am3,am,imass,iterl,iter2
      close(25)
c set a new seed by hand...
      if(iseed.ne.0)then
      seed=4587658465895.0
      endif
      write(7,*) 'seed: ', seed
      call ranset(seed)
      idum=-1
      y=rano(yran,idum)
      write(6,*) 'rano: ', y,yran,idum
c*******************************************************************
c     initialization
c     istart.lt.0 : start from tape
c     istart=0    : ordered start
c     istart=1    : random start
c*******************************************************************
      call init(istart)
c  read in Remez coefficients
      read(36,*) anum2(0)
      read(37,*) anum4(0)
      read(38,*) anum2g(0)
      do i=1,ndiag
      read(36,*) anum2(i),aden2(i)
      read(37,*) anum4(i),aden4(i)
      enddo
      do i=1,ndiagg
      read(38,*) anum2g(i),aden2g(i)
      enddo
      read(36,*) bnum2(0)
      read(37,*) bnum4(0)
      read(38,*) bnum2g(0)
      do i=1,ndiag
      read(36,*) bnum2(i),bden2(i)
      read(37,*) bnum4(i),bden4(i)
      enddo
      do i=1,ndiagg
      read(38,*) bnum2g(i),bden2g(i)
      enddo
c*******************************************************************
c     print heading
c*******************************************************************
      traj=iterl*dt
      proby=1.0/float(iterl)
c     write(6, 9001)ksize,ksizet,kthird,dt,traj,iter2,beta,am3,am,
c    & imass
      write(7, 9001)ksize,ksizet,kthird,dt,traj,iter2,beta,am3,am,
     & imass
9001  format(' ksize=',i3,' ksizet=',i3,/
     1 ,' kthird=',i3,/
     1 ,' time step: dt=',f6.4,' trajectory length=',f9.6,/
     1 ,' # trajectories=',i6,' beta=',f9.6,/
     1 ,' am3=',f6.4,' am=',f6.4/
     1 ,' imass=',i2)
c     write(6,9004) rescgg,rescga,respbp
      write(7,9004) rescgg,rescga,respbp
9004  format(' Stopping residuals: guidance: ',e9.4,' acceptance: ',
     &     e9.4,' estimator: ',e9.4)
      call ranget(seed)
c     write(6,*) 'seed: ', seed
      write(7,*) 'seed: ', seed
c*******************************************************************
c       initialize for averages
c*******************************************************************
      actiona=0.0
      vel2a=0.0
      condmata=0.0
      susca=0.0
      do ia=1,3
      conda(ia)=0.0
      enddo
      pbpa=0.0
      ancg=0.0
      ancgh=0.0
      ancgf=0.0
      ancgpf=0.0
      ancgma=0.0
      yav=0.0
      yyav=0.0 
      naccp=0
      ipbp=0
      itot=0
c*******************************************************************
c     start of classical evolution
c*******************************************************************
      do 601 isweep=1,iter2
c uncomment line below to go straight to measurement
c     goto 666
c*******************************************************************
c     initialise trial fields
c*******************************************************************
      do 2007 ia = 1,3
      do 2007 i=1,kvol
      chit(i,ia)=chi(i,ia)
2007  continue
c
c*******************************************************************
c  Pseudofermion fields: Phi =  (MdaggerM)^0.25  * R, where R is gaussian
c*******************************************************************
      do iflavor=1,2
      do idirac=1,4
      do ithird=1,kthird
      call gauss0(ps)
      do i=1,kvol
      R(ithird,i,idirac,iflavor)=cmplx(ps(i,1),ps(i,2))
      enddo
      enddo
      enddo
      enddo
c  For Phi = {MdaggerM}^0.25 * R
c
      call qmrherm(R,rescga,itercg,am,imass,anum4,aden4,ndiag,0)
      ancgpf=ancgpf+float(itercg)
c
      do iflavor=1,2
      do idirac=1,4
      do i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac,iflavor)=X1(ithird,i,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
c*******************************************************************
c     heatbath for p  (Bit clunky for now) 
c*******************************************************************
      call gaussp(ps)
      do i=1,kvol
      pp(i,1)=ps(i,1)
      pp(i,2)=ps(i,2)
      enddo
      call gaussp(ps)
      do i=1,kvol
      pp(i,3)=ps(i,1)
      enddo
c*******************************************************************
c  call to Hamiltonian
c     
      call hamilton(Phi,
     &      H0,hg,hp,S0,rescga,isweep,am,imass)
      if(isweep.eq.1) then
      action=S0/kvol
      gaction=hg/kvol
      paction=hp/kvol
      endif 
c     goto 501
c*******************************************************************
c      half-step forward for p
c*******************************************************************
      call force(Phi,rescgg,am,imass)
      d=dt*0.5
      do 2004 ia=1,3
      do 2004 i=1,kvol
      pp(i,ia)=pp(i,ia)-dSdchi(i,ia)*d
2004  continue
c*******************************************************************
c     start of main loop for classical time evolution
c*******************************************************************
      do 500 iter=1,itermax
c
c  step (i) st(t+dt)=st(t)+p(t+dt/2)*dt;
c
      d=dt
      do 2001 ia=1,3
      do 2001 i=1,kvol
      chit(i,ia)=chit(i,ia)+d*pp(i,ia)
2001  continue
c
c  step (ii)  p(t+3dt/2)=p(t+dt/2)-dSds(t+dt)*dt (1/2 step on last iteration)
c
      call force(Phi,rescgg,am,imass)
c
c test for end of random trajectory
c 
c     ytest=ranf()
      ytest=rano(yrun,idum)
      if(ytest.lt.proby)then
      d=dt*0.5
      do 2005 ia=1,3
      do 2005 i=1,kvol
      pp(i,ia)=pp(i,ia)-d*dSdchi(i,ia)
2005  continue
      itot=itot+iter 
      goto 501
      else
      d=dt
      do 3005 ia=1,3
      do 3005 i=1,kvol
      pp(i,ia)=pp(i,ia)-d*dSdchi(i,ia)
3005  continue
      endif
c 
500   continue
c**********************************************************************
c  Monte Carlo step: accept new fields with probability=
c              min(1,exp(H0-H1))
c**********************************************************************
501   continue 
      call hamilton(Phi,
     &       H1,hg,hp,S1,rescga,isweep,am,imass)
      dH=H0-H1
      dS=S0-S1
      write(98,*) dH,dS
      y=exp(dH)      
      yav=yav+y 
      yyav=yyav+y*y 
      if(dH.lt.0.0)then
c     x=ranf()
      x=rano(yrun,idum)
      if(x.gt.y)goto 600
      endif
c
c     step accepted: set s=st
c
      do 2006 ia=1,3
      do 2006 i=1,kvol
      chi(i,ia)=chit(i,ia)
2006  continue      
      naccp=naccp+1
      action=S1/kvol
      gaction=hg/kvol
      paction=hp/kvol
600   continue
      write(11,*) isweep,gaction,paction
      actiona=actiona+action 
      vel2=0.0
      do 457 ia=1,3
      do 457 i=1,kvol
      vel2=vel2+pp(i,ia)*pp(i,ia)
457   continue
      vel2=vel2/kvol
      vel2a=vel2a+vel2
c
c     uncomment to disable measurements
c     goto 601
666   if((isweep/iprint)*iprint.eq.isweep)then
      do 2066 ia=1,3
      do 2066 i=1,kvol
      chit(i,ia)=chi(i,ia)
2066  continue      
c
c  Here are the observables related to the auxiliary scalar field 
c    PLEASE SANITY CHECK
      susctot=0.0
      condmod=0.0
      do ia=1,3
        cond(ia)=0.0
        susc(ia)=0.0
        do i=1,kvol
          cond(ia)=cond(ia)+chi(i,ia)
        enddo
        susc(ia)=cond(ia)*cond(ia)
        cond(ia)=cond(ia)/kvol
        condmod=condmod+cond(ia)*cond(ia)
        susc(ia)=susc(ia)/kvol
        susctot=susctot+susc(ia)
      enddo
      condmod=sqrt(condmod)
      condmoda=condmoda+condmod
      susca=susca+susctot
      do ia=1,3
      conda(ia)=conda(ia)+cond(ia)
      enddo
      write(12,*) (cond(ia), ia=1,3)
      write(13,*) (susc(ia), ia=1,3)
      write(14,*) condmod, susctot
c
c  I don't think we're ever going to use this.....
c     call measure(pbp,respbp,ancgm,am,imass)
      pbpa=pbpa+pbp
      ancgma=ancgma+ancgm
      ipbp=ipbp+1
c     write(11,*) pbp
      write(6,*) isweep,':  ',condmod,susctot
      endif
c
      if((isweep/icheck)*icheck.eq.isweep)then
      call ranget(seed)
      call swrite
      flush(12)
      flush(13)
      flush(14)
c     flush(100)
c     flush(200)
c     write(7,9023) seed
      endif
c
601   continue
c*******************************************************************
c     end of main loop
c*******************************************************************
      actiona=actiona/iter2 
      vel2a=vel2a/iter2 
      condmoda=condmoda/ipbp
      susca=susca/ipbp
      do ia=1,3
      conda(ia)=conda(ia)/ipbp
      enddo
      pbpa=pbpa/ipbp
      ancg=ancg/itot
      ancgf=ancgf/itot
      ancgh=ancgh/(2*iter2)
      ancgpf=ancgpf/iter2
      ancgma=ancgma/ipbp
      yav=yav/iter2
      yyav=yyav/iter2-yav*yav 
      yyav=sqrt(yyav/(iter2-1)) 
      atraj=dt*itot/iter2 
c*******************************************************************
c     print global averages
c*******************************************************************
c     write(6, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgh,ancgma,
c    & pbpa,vel2a,actiona
      write(7, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgh,
     & ancgf,ancgpf,ancgma,condmoda,susca,conda(1),conda(2),conda(3),
     & pbpa,vel2a,actiona
9022  format(' averages for last ',i6,' trajectories',/ 
     & 1x,' # of acceptances: ',i6,' average trajectory length= ',f8.3/
     & 1x,' <exp-dH>=',e11.4,' +/-',e10.3/
     2 1x,' av. # QMR itr. guidance: ',f9.3,'; acceptance: ',f9.3/
     & 1x,' force: ',f9.3,'; pseudofermion: ',f9.3//
     & 1x,' measurement: ',f9.3//
     1 1x,' condmod=',e11.3,' susc=',e11.3/
     & 1x 'cond components=',3e11.3/
     1 1x,' psibarpsi=',e11.3/
     & 1x,' mean square velocity=',e10.3,'; action per site=',e10.3//)
      write(7, 9024)
      write(7, 9024)
9024  format(1x)
c
      close(11)
c
      if(iwrite.eq.1) then
      call ranget(seed)
      call swrite
      write(7,*) 'seed: ', seed
      endif
c
      stop
      end
c******************************************************************
c   calculate dSds for sigma fields at each intermediate time
c******************************************************************
      subroutine force(Phi,res1,am,imass)
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksize*ksize*ksizet)
      parameter(kferm=kvol*kthird*4*2)
      parameter(ndiagg=12,ndiag=ndiagg)
      common/trial/chi(kvol,3),pp(kvol,3)
      common/neighb/id(kvol,3),iu(kvol,3)
      common/sforce/dSdchi(kvol,3)
      common/para/beta,am3,ibound
      common/param/ancg,ancgh,ancgf,ancgpf
      common/vector/X1(kferm)
c     common/dum1/ R(kthird,kvol,4,2),ps(kvol,2)
      common/remez2g/anum(0:ndiagg),aden(ndiagg),
     &              bnum(0:ndiagg),bden(ndiagg)
      complex Phi(kferm), X2(kferm)
      complex X1
      real*8 anum,aden,bnum,bden
c
c     write(6,111)
111   format(' Hi from force')
c
c     factor of 2 because fermions are 8 component objects
c     another factor of 2 from differentiating chi squared
c     and a factor of 0.5 due to N=1
      do ia=1,3
      do i=1,kvol
         dSdchi(i,ia)=2*beta*chi(i,ia)
      enddo
      enddo
c  pseudofermion action is
c   Phi^dagger {MdaggerM(m)})^-1/2 Phi
c
      do i=1,kferm
      X2(i)=Phi(i)
      enddo
      call qmrherm(X2,res1,itercg,am,imass,bnum,bden,ndiag,0)
      ancg=ancg+float(itercg)
      do i=1,kferm
      X2(i)=X1(i)
      enddo
c
c  evaluates X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
      call qmrherm(X2,res1,itercg,am,imass,anum,aden,ndiag,1)
      ancgf=ancgf+float(itercg)
c
      return
      end
c******************************************************************
c   Evaluation of Hamiltonian function
c******************************************************************
      subroutine hamilton(Phi,
     &       h,hg,hp,s,res2,isweep,am,imass)
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*2*kvol*kthird)
      parameter(ndiag=25)
      parameter(One=1.0)
      common/trial/chi(kvol,3),pp(kvol,3)
      common/remez2/anum(0:ndiag),aden(ndiag),bnum(0:ndiag),bden(ndiag)
      common /para/beta,am3,ibound
      common/param/ancg,ancgh,ancgf,ancgpf
      common/vector/ X1(kferm)      
      common/dum1/R(kferm),ps(kvol,2)
      complex Phi(kferm)
      complex X1,R
      real*8 hp,hg,hf,h,s
      real*8 anum,aden,bnum,bden
c     write(6,111)
111   format(' Hi from hamilton')
c
      hp=0.0
      hg=0.0
      hf=0.0
c
      do 22 ia=1,3
      do 22 i=1,kvol
      hp=hp+pp(i,ia)*pp(i,ia)
22    continue
c
      hp=hp*0.5
c
      do ia=1,3
      do i=1,kvol
      hg=hg+chi(i,ia)*chi(i,ia)
      enddo
      enddo
c     hg=2*beta*hg
      hg=beta*hg
      h=hg+hp
c 
c uncomment these lines to quench the fermions!
c     write(6,*) isweep,':  hg', hg,'   hp', hp,'   h',h
c     return
c
c  pseudofermion action is
c   Phi^dagger {MdaggerM(m)})^-1/2 Phi
c
      do i=1,kferm
      R(i)=Phi(i)
      enddo
      call qmrherm(R,res2,itercg,am,imass,bnum,bden,ndiag,0)
c
      ancgh=ancgh+float(itercg)

c         
      do 4 iferm=1,kferm
      hf=hf+conjg(Phi(iferm))*X1(iferm)
4     continue
c
      h=hg+hp+hf
      write(6,*) isweep,':  hg', hg,'   hp', hp,'   hf', hf,
     &   '   h',h
      s=hg+hf
c
      return
      end              
c******************************************************************
c    multisolver matrix inversion via Lanczos technique
c  eg. Golub & van Loan "Matrix Computations" 9.3.1
c       solves (MdaggerM+diag)*x=Phi for ndiag different values of diag
c   eg. Frommer et al hep-lat/9504020
c   iflag=0: simply evaluates X = {MdaggerM}^p * Phi
c   can be interchanged with congrad for p=-1
c   iflag=1: evaluates DWF force term
c*****************************************************************m
      subroutine qmrherm(Phi,res,itercg,am,imass,anum,aden,ndiag,iflag)
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*2*kthird*kvol)
      parameter(niterc=2000)
      common/trial/chi(kvol,3),pp(kvol,3)
      common/sforce/dSdchi(kvol,3)
      common/para/bbb,am3,ibound
      common/vector/x(kferm)
      complex Phi(kferm)
      complex x
      real*8 anum(0:ndiag),aden(ndiag)
      complex vtild(kferm),q(kferm),pm1(kferm,ndiag)
      complex qm1(kferm),p(kferm,ndiag),x3(kferm),r(kferm)
      complex x1(kferm,ndiag),x2(kferm)
      real alpha(ndiag),beta
      real amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
c
c     write(6,111)
111   format(' Hi from qmrherm')
c
      resid=sqrt(kferm*res*res)
c     write(6,*) resid
      itercg=0
c
c   initialise r=Phi
c
      do i=1,kferm
         r(i)=Phi(i)
         qm1(i)=(0.0,0.0)
         x(i)=anum(0)*Phi(i)
      enddo
      beta=0.0
      do i=1,kferm
         beta=beta + real(r(i))**2 + aimag(r(i))**2
      enddo
      beta=sqrt(beta)
      phimod=beta
c     write(6,*) '|| Phi || = ', phimod
c
      do iter=1,niterc
      itercg=itercg+1
c
c  Lanczos steps
c
      do i=1,kferm
         q(i)=r(i)/beta
      enddo
c
      call dslash(vtild,q,chi,am,imass)
      call dslashd(x3,vtild,chi,am,imass)
c
      alphatild=0.0
      do i=1,kferm
         alphatild=alphatild+conjg(q(i))*x3(i)
      enddo
c
      do i=1,kferm
         r(i)=x3(i)-alphatild*q(i)-beta*qm1(i)
         qm1(i)=q(i)
      enddo
c
      beta0=beta
      beta=0.0
      do i=1,kferm
         beta=beta + real(r(i))**2 + aimag(r(i))**2
      enddo
      beta=sqrt(beta)
c
      do idiag=1,ndiag
      alpha(idiag)=alphatild+aden(idiag)
      enddo
c
      if(iter.eq.1)then
           do idiag=1,ndiag
             d(idiag)=alpha(idiag)
             do i=1,kferm
              p(i,idiag)=q(i)
              pm1(i,idiag)=p(i,idiag)
             enddo
             rho(idiag)=beta0/alpha(idiag)
             rhom1(idiag)=rho(idiag)
             do i=1,kferm
              x1(i,idiag)=rho(idiag)*q(i)
             enddo
           enddo
      else
           do idiag=1,ndiag
             amu(idiag)=beta0/d(idiag)
             dm1(idiag)=d(idiag)
             d(idiag)=alpha(idiag)-beta0*amu(idiag)
             do i=1,kferm
              p(i,idiag)=q(i)-amu(idiag)*pm1(i,idiag)
              pm1(i,idiag)=p(i,idiag)
             enddo
             rho(idiag)=-amu(idiag)*dm1(idiag)*rhom1(idiag)/d(idiag)
c  Convergence criterion (a bit ad hoc for now...)
             if(idiag.eq.1)then
               rhomax=abs(phimod*rho(idiag))
             else
               if(abs(phimod*rho(idiag)).gt.rhomax)
     &             rhomax=abs(phimod*rho(idiag))
             endif
             rhom1(idiag)=rho(idiag)
             do i=1,kferm
               x1(i,idiag)=x1(i,idiag)+rho(idiag)*p(i,idiag)
             enddo
           enddo
c  check to see whether the residual is acceptable for all ndiag....
c  criterion is a bit ad hoc
c  relaxing by a factor arelax improves stability and leads to quicker
c  convergence
       arelax=2.0
           if(rhomax.lt.arelax*resid) then
c     call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiag)
             goto 8
           endif
      endif
c
c   end of loop over iter
      enddo
8     continue
c
      if(iflag.eq.0)then
c
c  evaluate solution x=(MdaggerM)^{-0.5}*Phi
      do idiag=1,ndiag
      do i=1,kferm
      x(i)=x(i)+anum(idiag)*x1(i,idiag)
      enddo
      enddo
c
      else
c
c compute contribution to force...
c
      do idiag=1,ndiag
c
c  X2 = M*X1
      do i=1,kferm
      R(i)=x1(i,idiag)
      enddo
      call dslash(X2,R,chi,am,imass)
      call derivs(R,X2,anum(idiag),imass)
c
      enddo
c
      endif
c
      return
      end
c**********************************************************************
      subroutine derivs(R,X2,anum,imass)
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksizet*ksize*ksize)
      common/sforce/dSdchi(kvol,3)
      complex R(kthird,kvol,4,2),X2(kthird,kvol,4,2)
      real*8 anum
c
c     write(6,111)
111   format(' Hi from derivs')

c     dSdchi=dSdchi-alpha_i*Re(Rdagger *(d(Mdagger)dp)* X2) 
c     Cf. Montvay & Muenster (7.215)
c     I don't yet understand the -sign :-(   SJH 31/1/23
      tzi=-2*anum
c     factor of 2 picks up second term in M&M (7.215)
c
      if(imass.eq.1)then
c
      do idirac=1,2
      do i=1,kvol
      dSdchi(i,1)=dSdchi(i,1)-real(tzi*
     &(conjg(R(kthird,i,idirac,1))*
     &        X2(1,i,idirac,2)
     &+conjg(R(kthird,i,idirac,2))*
     &        X2(1,i,idirac,1)))
      dSdchi(i,2)=dSdchi(i,2)-aimag(tzi*
     &(conjg(R(kthird,i,idirac,1))*
     &        X2(1,i,idirac,2)
     &-conjg(R(kthird,i,idirac,2))*
     &        X2(1,i,idirac,1)))
      dSdchi(i,3)=dSdchi(i,3)-real(tzi*
     &(conjg(R(kthird,i,idirac,1))*
     &        X2(1,i,idirac,1)
     &-conjg(R(kthird,i,idirac,2))*
     &        X2(1,i,idirac,2)))
      enddo
      enddo
c
      do idirac=3,4
      do i=1,kvol
      dSdchi(i,1)=dSdchi(i,1)-real(tzi*
     &(conjg(R(1,i,idirac,1))*
     &        X2(kthird,i,idirac,2)
     &+conjg(R(1,i,idirac,2))*
     &        X2(kthird,i,idirac,1)))
      dSdchi(i,2)=dSdchi(i,2)-aimag(tzi*
     &(conjg(R(1,i,idirac,1))*
     &        X2(kthird,i,idirac,2)
     &-conjg(R(1,i,idirac,2))*
     &        X2(kthird,i,idirac,1)))
      dSdchi(i,3)=dSdchi(i,3)-real(tzi*
     &(conjg(R(1,i,idirac,1))*
     &        X2(kthird,i,idirac,1)
     &-conjg(R(1,i,idirac,2))*
     &        X2(kthird,i,idirac,2)))
      enddo
      enddo
c
      elseif(imass.eq.3)then
c
      do idirac=1,2
      do i=1,kvol
      dSdchi(i,1)=dSdchi(i,1)-aimag(tzi*
     &(conjg(R(kthird,i,idirac,1))*
     &        X2(1,i,idirac,2)
     &+conjg(R(kthird,i,idirac,2))*
     &        X2(1,i,idirac,1)))
      dSdchi(i,2)=dSdchi(i,2)+real(tzi*
     &(conjg(R(kthird,i,idirac,1))*
     &        X2(1,i,idirac,2)
     &-conjg(R(kthird,i,idirac,2))*
     &        X2(1,i,idirac,1)))
      dSdchi(i,3)=dSdchi(i,3)-aimag(tzi*
     &(conjg(R(kthird,i,idirac,1))*
     &        X2(1,i,idirac,1)
     &-conjg(R(kthird,i,idirac,2))*
     &        X2(1,i,idirac,2)))
      enddo
      enddo
c
      do idirac=3,4
      do i=1,kvol
      dSdchi(i,1)=dSdchi(i,1)+aimag(tzi*
     &(conjg(R(1,i,idirac,1))*
     &        X2(kthird,i,idirac,2)
     &+conjg(R(1,i,idirac,2))*
     &        X2(kthird,i,idirac,1)))
      dSdchi(i,2)=dSdchi(i,2)-real(tzi*
     &(conjg(R(1,i,idirac,1))*
     &        X2(kthird,i,idirac,2)
     &-conjg(R(1,i,idirac,2))*
     &        X2(kthird,i,idirac,1)))
      dSdchi(i,3)=dSdchi(i,3)+aimag(tzi*
     &(conjg(R(1,i,idirac,1))*
     &        X2(kthird,i,idirac,1)
     &-conjg(R(1,i,idirac,2))*
     &        X2(kthird,i,idirac,2)))
      enddo
      enddo

      endif
c
      return
      end
c******************************************************************
c    tests inversion for each term in partial fraction expansion
c    accuracy severely curtailed by truncation/rounding - use for testing purposes only.
c******************************************************************
      subroutine testinv(Phi,resmax,itercg,am,imass,x,aden,ndiag)
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*2*kthird*kvol)
      common/trial/chi(kvol,3),pp(kvol,3)
      common/para/bbb,am3,ibound
      complex Phi(kferm)
      complex vtild(kferm)
      complex x3(kferm)
      complex x(kferm,ndiag),x1(kferm),x2(kferm)
      real*8 aden(ndiag)
c
c     write(6,111)
111   format(' Hi from testinv')
c
      resmax=0.0
c
      do idiag=1,ndiag
      residual=0.0
      do i=1,kferm
      x3(i)=x(i,idiag)
      enddo
      call dslash(x2,x3,chi,am,imass)
      call dslashd(x1,x2,chi,am,imass)
      do i=1,kferm
      vtild(i)=x1(i)+aden(idiag)*x3(i)-Phi(i)
      residual=residual + real(vtild(i))**2 + aimag(vtild(i))**2
      enddo
      if(residual.gt.resmax) resmax=residual
c
c     write(6,*) idiag, 'itercg = ',itercg, ' residual = ',residual
      enddo
c
      resmax=sqrt(resmax)
c
      return
      end
c******************************************************************
c    matrix inversion via conjugate gradient algorithm
c       solves (Mdagger)Mx=Phi, 
c           NB. no even/odd partitioning
c******************************************************************
      subroutine congrad(Phi,res,itercg,am,imass)
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*2*kthird*kvol)
      parameter(niterc=kthird*kvol)
      common/trial/chi(kvol,3),pp(kvol,3)
      common/para/bbb,am3,ibound
      common/vector/x(kferm)
      complex Phi(kferm)
      complex x
      complex x1(kferm),x2(kferm),p(kferm),r(kferm)
c     write(6,111)
111   format(' Hi from congrad')
c
      resid=kferm*res*res
      itercg=0
c
      do 1 nx=1,niterc
      itercg=itercg+1
      if(nx.gt.1) goto 51
c
c   initialise p=x, r=Phi(na)
c
      do 2 i=1,kferm
      p(i)=x(i)
      r(i)=Phi(i)
2     continue
      betad=1.0
      alpha=1.0
51    alphad=0.0
c
c  x1=Mp
c
      call dslash(x1,p,chi,am,imass)
c
      if(nx.eq.1) goto 201
c
c   alpha=(r,r)/(p,(Mdagger)Mp)
c 
      alphad=0.0
      do 31 i=1,kferm
      alphad=alphad + real(x1(i))**2 + aimag(x1(i))**2
31    continue
      alpha=alphan/alphad
c      
c   x=x+alpha*p
c
      do 4 i=1,kferm
      x(i)=x(i)+alpha*p(i)
4     continue
201   continue
c     
c   x2=(Mdagger)x1=(Mdagger)Mp
c
      call dslashd(x2,x1,chi,am,imass)
c
c   r=r-alpha*(Mdagger)Mp
c
      do 6 i=1,kferm
      r(i)=r(i)-alpha*x2(i)
6     continue
c
c   beta=(r_k+1,r_k+1)/(r_k,r_k)
c
      betan=0.0 
      do 61 i=1,kferm
      betan=betan + real(r(i))**2 + aimag(r(i))**2 
61    continue 
      beta=betan/betad
      betad=betan
      alphan=betan
c
      if(nx.eq.1) beta=0.0
c
c   p=r+beta*p
c
      do 7 i=1,kferm
      p(i)=r(i)+beta*p(i)
7     continue
      if(betan.lt.resid) goto 8
1     continue
      write(6,1000)
      write(7,1000)
1000  format(' # iterations of congrad exceeds niterc')
8     continue
c     write(6,*) itercg
      return
      end      
c******************************************************************
c   Calculate fermion expectation values via a noisy estimator
c   -matrix inversion via conjugate gradient algorithm
c       solves Mx=x1
c     (Numerical Recipes section 2.10 pp.70-73)   
c*******************************************************************
      subroutine measure(psibarpsi,res,aviter,am,imass)
      parameter(ksize=12,ksizet=12,kthird=8,kvol=ksizet*ksize*ksize)
      parameter(akappa=0.5)
      parameter(knoise=10)
      common/trial/chi(kvol,3),pp(kvol,3)
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      common /neighb/id(kvol,3),iu(kvol,3)
      common/vector/xi(kthird,kvol,4,2)
      complex x(kvol,4,2), Phi(kthird,kvol,4,2)
      complex xi,gamval
      complex psibarpsi1,psibarpsi2
      real ps(kvol,2),pt(kvol,2)
      real pbp(knoise)
      integer gamin
c     write(6,*) 'hi from measure'
c
      iter=0
      pbp=0.0
      do inoise=1,knoise
c
c     set up noise
      call gauss0(ps)
      psibarpsi1=(0.0,0.0)
      call gauss0(pt)
      psibarpsi2=(0.0,0.0)
c
      do idsource=1,2
c
c  source on domain wall at ithird=1
c   
      if(imass.ne.5)then
         do iflavor=1,2
         do idirac=1,4
           if(idirac.eq.idsource)then
             do i=1,kvol
               x(i,idirac,iflavor)=cmplx(ps(i,1),ps(i,2))
             enddo
           else
             do i=1,kvol
               x(i,idirac,iflavor)=(0.0,0.0)
             enddo
           endif
         enddo
         enddo
      else
         do iflavor=1,2
         do idirac=1,4
           if(idirac.eq.idsource+2)then
             do i=1,kvol
               x(i,idirac,iflavor)=cmplx(ps(i,1),ps(i,2))
             enddo
           else
             do i=1,kvol
               x(i,idirac,iflavor)=(0.0,0.0)
             enddo
           endif
         enddo
         enddo
      endif
c
c
      if(imass.ne.5)then
         do ithird=1,kthird
         if(ithird.eq.1)then
            do iflavor=1,2
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac,iflavor)=x(i,idirac,iflavor)
            enddo
            enddo
            enddo
         else
            do iflavor=1,2
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac,iflavor)=(0.0,0.0)
            enddo
            enddo
            enddo
         endif
         enddo
      else
c   xi = 0.5(1+gamma_4)*gamma_5*eta on DW at ithird=1
         do ithird=1,kthird
         if(ithird.eq.1)then
            do iflavor=1,2
            do idirac=1,2
            do i=1,kvol
               xi(ithird,i,idirac,iflavor)=-x(i,idirac+2,iflavor)
               xi(ithird,i,idirac+2,iflavor)=(0.0,0.0)
            enddo
            enddo
            enddo
         else
c   xi = 0.5(1+gamma_4)*eta on DW at ithird=1
            do iflavor=1,2
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac,iflavor)=(0.0,0.0)
            enddo
            enddo
            enddo
         endif
         enddo
      endif
c
c Phi= Mdagger*xi
c
      call dslashd(Phi,xi,chi,am,imass)
c
      do 302 iflavor=1,2
      do 302 idirac=1,4
      do 302 i=1,kvol
      do 302 ithird=1,kthird
      xi(ithird,i,idirac,iflavor)=Phi(ithird,i,idirac,iflavor)
302   continue
c  
c xi= (MdaggerM)**-1 * Phi 
c
      call congrad(Phi,res,itercg,am,imass) 
      iter=iter+itercg
c
      if(imass.ne.5)then
c  pbp1 = x^dagger (0.5(1+gamma_4)) xi(kthird)
         do iflavor=1,2
         do i=1,kvol
         psibarpsi1=psibarpsi1
     &    +conjg(x(i,idsource,iflavor))*xi(kthird,i,idsource,iflavor)
         enddo
         enddo
      else
c  pbp1 = x^dagger (0.5(1-gamma_4)) xi(1)
         do iflavor=1,2
         do i=1,kvol
         psibarpsi1=psibarpsi1
     &    +conjg(x(i,idsource+2,iflavor))*xi(1,i,idsource+2,iflavor)
         enddo
         enddo
      endif
c
c
c  source on domain wall at ithird=kthird
      idsource2=idsource+2
c
      if(imass.ne.5)then
         do iflavor=1,2
         do idirac=1,4
         do i=1,kvol
         if(idirac.eq.idsource2)then
             x(i,idirac,iflavor)=cmplx(pt(i,1),pt(i,2))
         else
             x(i,idirac,iflavor)=(0.0,0.0)
         endif
         enddo
         enddo
         enddo
      else
         do iflavor=1,2
         do idirac=1,4
         do i=1,kvol
         if(idirac.eq.idsource2-2)then
             x(i,idirac,iflavor)=cmplx(pt(i,1),pt(i,2))
         else
             x(i,idirac,iflavor)=(0.0,0.0)
         endif
         enddo
         enddo
         enddo
      endif
c
      if(imass.ne.5)then
c   xi = 0.5(1-gamma_4)*eta on DW at ithird=kthird
         do ithird=1,kthird
         if(ithird.eq.kthird)then
            do iflavor=1,2
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac,iflavor)=x(i,idirac,iflavor)
            enddo
            enddo
            enddo
         else
            do iflavor=1,2
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac,iflavor)=(0.0,0.0)
            enddo
            enddo
            enddo
         endif
         enddo
      else
c   xi = 0.5(1-gamma_4)*gamma_5*eta on DW at ithird=kthird
         do ithird=1,kthird
         if(ithird.eq.kthird)then
            do iflavor=1,2
            do idirac=1,2
            do i=1,kvol
               xi(ithird,i,idirac+2,iflavor)=-x(i,idirac,iflavor)
               xi(ithird,i,idirac,iflavor)=(0.0,0.0)
            enddo
            enddo
            enddo
         else
            do iflavor=1,2
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac,iflavor)=(0.0,0.0)
            enddo
            enddo
            enddo
         endif
         enddo
      endif
c
c Phi= Mdagger*xi
c
      call dslashd(Phi,xi,chi,am,imass)
c
      do 305 iflavor=1,2
      do 305 idirac=1,4
      do 305 i=1,kvol
      do 305 ithird=1,kthird
      xi(ithird,i,idirac,iflavor)=Phi(ithird,i,idirac,iflavor)
305   continue
c  
c xi= (MdaggerM)**-1 * Phi 
c
      call congrad(Phi,res,itercg,am,imass) 
      iter=iter+itercg
c
      if(imass.ne.5)then
c pbp2= - x^dagger (0.5(1-gamma_4)) xi(1)
         do iflavor=1,2
         do i=1,kvol
         psibarpsi2=psibarpsi2
     &     +conjg(x(i,idsource2,iflavor))*xi(1,i,idsource2,iflavor)
         enddo
         enddo
      else
c pbp2= - x^dagger (0.5(1-gamma_4)) xi(kthird)
         do iflavor=1,2
         do i=1,kvol
         psibarpsi2=psibarpsi2
     &     +conjg(x(i,idsource,iflavor))*xi(kthird,i,idsource,iflavor)
         enddo
         enddo
      endif
c
c  end trace on Dirac indices....
      enddo
c
      if(imass.eq.1)then
         psibarpsi1=psibarpsi1/kvol
         psibarpsi2=psibarpsi2/kvol
         pbp(inoise)=psibarpsi1+psibarpsi2
      elseif(imass.eq.3)then
         psibarpsi1=(0.0,-1.0)*psibarpsi1/kvol
         psibarpsi2=(0.0,+1.0)*psibarpsi2/kvol
         pbp(inoise)=psibarpsi1+psibarpsi2
      elseif(imass.eq.5)then
         psibarpsi1=(0.0,-1.0)*psibarpsi1/kvol
         psibarpsi2=(0.0,-1.0)*psibarpsi2/kvol
         pbp(inoise)=psibarpsi1+psibarpsi2
      endif
c        write(6,*) real(psibarpsi1),aimag(psibarpsi1),
c    &       real(psibarpsi2),aimag(psibarpsi2)
         write(100,*) real(psibarpsi1),aimag(psibarpsi1),
     &       real(psibarpsi2),aimag(psibarpsi2)
c
c end loop on noise
      enddo
c
      psibarpsi=0.0
      susclsing=0.0
c
      do inoise=1,knoise
      psibarpsi=psibarpsi+pbp(inoise)
         do jnoise=knoise,inoise+1,-1
         susclsing=susclsing+pbp(inoise)*pbp(jnoise)
         enddo
      enddo
      psibarpsi=psibarpsi/knoise
      susclsing=2*kvol*susclsing/(knoise*(knoise-1))
         write(200,*) psibarpsi,susclsing
      aviter=float(iter)/(4*knoise)
      return
      end
c*****************************************************************
      subroutine sread
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
      subroutine swrite
      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
      common/chi/chi(kvol,3),seed
      real*8 seed
      open(unit=31,file='con',
     1     status='unknown',form='unformatted')
      write (31) chi,seed
      close(31)
      return
      end
c
      subroutine init(nc)
c*******************************************************************
c     sets initial values
c     nc=0 cold start
c     nc=1 hot start
c     nc<0 no initialization
c*******************************************************************
      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
      parameter(akappa=0.5)
      common/chi/chi(kvol,3),seed
      common /neighb/id(kvol,3),iu(kvol,3)
      common/dirac/gamval(6,4),gamin(6,4)
      common/ranseed/yran,idum
      common/v/v(97)
      complex gamval,one,zi
      real rano
      real*8 seed
      integer gamin
c
c
      one=(1.0,0.0)
      zi=(0.0,1.0)
c*******************************************************************
c  calculate constants
c*******************************************************************
      call addrc
c*******************************************************************
c    setup Dirac algebra
c*******************************************************************
c
c     gamma_1
c
      gamval(1,1)=-zi
      gamval(1,2)=-zi
      gamval(1,3)= zi
      gamval(1,4)= zi
c
      gamin(1,1)=4
      gamin(1,2)=3
      gamin(1,3)=2
      gamin(1,4)=1
c
c     gamma_2
c
      gamval(2,1)=-one
      gamval(2,2)= one
      gamval(2,3)= one
      gamval(2,4)=-one
c
      gamin(2,1)=4
      gamin(2,2)=3
      gamin(2,3)=2
      gamin(2,4)=1
c
c     gamma_3
c
      gamval(3,1)=-zi
      gamval(3,2)= zi
      gamval(3,3)= zi
      gamval(3,4)=-zi
c
      gamin(3,1)=3
      gamin(3,2)=4
      gamin(3,3)=1
      gamin(3,4)=2
c
c     gamma_4
c
      gamval(4,1)= one
      gamval(4,2)= one
      gamval(4,3)= -one
      gamval(4,4)= -one
c
      gamin(4,1)=1
      gamin(4,2)=2
      gamin(4,3)=3
      gamin(4,4)=4
c
c     gamma_5 = gamma_1 * gamma_2 * gamma_3 * gamma_4
c
      gamval(5,1)=-one
      gamval(5,2)=-one
      gamval(5,3)=-one
      gamval(5,4)=-one
c
      gamin(5,1)=3
      gamin(5,2)=4
      gamin(5,3)=1
      gamin(5,4)=2
c
c     gamma_4 * gamma_5 (called gamma_3 gamma_5 in notes)
      gamval(6,1)=-one
      gamval(6,2)=-one
      gamval(6,3)= one
      gamval(6,4)= one
c
      gamin(6,1)=3
      gamin(6,2)=4
      gamin(6,3)=1
      gamin(6,4)=2
c
c
      do idirac=1,4
      do mu=1,6
      gamval(mu,idirac)=gamval(mu,idirac)*akappa
      enddo
      enddo
c
      if(nc.lt.0) return
c
c     initialize chi fields
c
      if(nc .eq. 1)goto 40
c     (else cold start)
      do 10 ia=1,3
      do 10 ind=1,kvol
      chi(ind,ia)=0.0
10    continue
      do 11 ind=1,kvol
      chi(ind,1)=0.3
11    continue
      return
c
40    continue
      g=0.2
      do 61 ia=1,3
      do 61 ind=1,kvol
      chi(ind,ia)=0.7+g*(2.0*rano(yran,idum)-1.0)
61    continue
      return
      end
c
      subroutine addrc
      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
c*******************************************************************
c
c     loads the addresses required during the update
c
c*******************************************************************
      common /neighb/id(kvol,3),iu(kvol,3)
      do 30 j3=1,ksizet
      do 30 j2=1,ksize
      do 30 j1=1,ksize
      ic=((j3-1)*ksize+(j2-1))*ksize+j1
      call ia(j1-1,j2,j3,id(ic,1))
      call ia(j1+1,j2,j3,iu(ic,1))
      call ia(j1,j2-1,j3,id(ic,2))
      call ia(j1,j2+1,j3,iu(ic,2))
      call ia(j1,j2,j3-1,id(ic,3))
      call ia(j1,j2,j3+1,iu(ic,3))
  30  continue
      return
      end
c
      subroutine ia(i1,i2,i3,nnn)
      parameter(ksize=12,ksizet=12)
c*******************************************************************
c
c     address calculator
c
c*******************************************************************
      n1=i1
      n2=i2
      n3=i3 
      if(n1) 2,2,3
   2  n1=n1+ksize
      go to 4
   3  if(n1-ksize) 4,4,5
   5  n1=n1-ksize
   4  if(n2) 6,6,7
   6  n2=n2+ksize
      go to 8
   7  if(n2-ksize) 8,8,9
   9  n2=n2-ksize
   8  if(n3) 10,10,11
  10  n3=n3+ksizet 
      go to 12
  11  if(n3-ksizet) 12,12,13
  13  n3=n3-ksizet   
  12  nnn=((n3-1)*ksize+(n2-1))*ksize+n1
      return
      end
c**********************************************************************
c calculate vector of gaussian random numbers with unit variance
c to refresh momenta
c   Numerical Recipes pp.203
c**********************************************************************
      subroutine gaussp(ps)
      parameter(ksize=12,ksizet=12,kvol=ksize*ksize*ksizet)
      common/trans/tpi 
      common/ranseed/yran,idum
      common/v/v(97)
      real rano
      real ps(kvol,2)
c     write(6,1)
1     format(' Hi from gaussp')
      do 1000 il=1,kvol
1000  ps(il,2)=sqrt(-2.0*log(rano(yran,idum)))
      do 1001 il=1,kvol
      theta=tpi*rano(yran,idum)
      ps(il,1)=ps(il,2)*sin(theta)
      ps(il,2)=ps(il,2)*cos(theta)
1001  continue 
      return
      end      
c**********************************************************************
c calculate vector of gaussian random numbers with unit variance
c to generate pseudofermion fields R
c   Numerical Recipes pp.203
c**********************************************************************
      subroutine gauss0(ps)
      parameter(ksize=12,ksizet=12,kvol=ksize*ksize*ksizet)
      common/trans/tpi 
      common/ranseed/yran,idum
      common/v/v(97)
      real rano
      real ps(kvol,2)
c     write(6,1)
1     format(' Hi from gauss0')
      do 1000 il=1,kvol
1000  ps(il,2)=sqrt(-log(rano(yran,idum)))
      do 1001 il=1,kvol
      theta=tpi*rano(yran,idum)
      ps(il,1)=ps(il,2)*sin(theta)
      ps(il,2)=ps(il,2)*cos(theta)
1001  continue 
      return
      end      
c*****************************************
c  Random number generator Numerical recipes 7.1
c
          real function rano(y,idum)
          common/v/v(97)
c
          if(idum.lt.0)then
               idum=1
               do j=1,97
                  dum=ranf()
               enddo
               do j=1,97
                  v(j)=ranf()
               enddo
               y=ranf()
          endif
c
          j=1+int(97.0*y)
          if(j.gt.97) j=97
          if(j.lt.1) j=1
c         write(6,*) j,y
c         write(6,*) 'problems with rano'
c         stop
c         endif
          y=v(j)
          rano=y
          v(j)=ranf()
          return
          end
C========================================================================
C
          SUBROUTINE RANGET(SEED)
          DOUBLE PRECISION    SEED,     G900GT,   G900ST,   DUMMY
          SEED  =  G900GT()
          RETURN
          ENTRY RANSET(SEED)
          DUMMY  =  G900ST(SEED)
          RETURN
          END
          REAL FUNCTION RANF()
          DOUBLE PRECISION    DRANF,    G900GT,   G900ST
          DOUBLE PRECISION    DS(2),    DM(2),    DSEED
          DOUBLE PRECISION    DX24,     DX48
          DOUBLE PRECISION    DL,       DC,       DU,       DR
          LOGICAL             SINGLE
          DATA      DS     /  1665 1885.D0, 286 8876.D0  /
          DATA      DM     /  1518 4245.D0, 265 1554.D0  /
          DATA      DX24   /  1677 7216.D0  /
          DATA      DX48   /  281 4749 7671 0656.D0  /
          SINGLE  =  .TRUE.
          GOTO 10
          ENTRY DRANF()
          SINGLE  =  .FALSE.
  10      DL  =  DS(1) * DM(1)
          DC  =  DINT(DL/DX24)
          DL  =  DL - DC*DX24
          DU  =  DS(1)*DM(2) + DS(2)*DM(1) + DC
          DS(2)  =  DU - DINT(DU/DX24)*DX24
          DS(1)  =  DL
          DR     =  (DS(2)*DX24 + DS(1)) / DX48
          IF(SINGLE)  THEN
             RANF  =  SNGL(DR)
          ELSE
             DRANF  =  DR
          ENDIF
          RETURN
          ENTRY G900GT()
          G900GT  =  DS(2)*DX24 + DS(1)
          RETURN
          ENTRY G900ST(DSEED)
          DS(2)  =  DINT(DSEED/DX24)
          DS(1)  =  DSEED - DS(2)*DX24
          G900ST =  DS(1)
          RETURN
          END
c***********************************************************************
      subroutine dslash(Phi,R,chi,am,imass)
c
c     calculates Phi = M*R
c
      parameter(ksize=12,ksizet=12,kthird=8)
      parameter(akappa=0.5)
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      real chi(ksize,ksize,ksizet,3)
      complex spm(ksize,ksize,ksizet)
      complex Phi(kthird,ksize,ksize,ksizet,4,2),
     &   R(kthird,ksize,ksize,ksizet,4,2)
      complex gamval
      integer gamin
c     write(6,*) 'hi from dslash'
c
c     diagonal term
      diag=(3.0-am3)+1.0
      do iflavor=1,2
      do idirac=1,4
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &  diag*R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c     Wilson term
      do iflavor=1,2
      do idirac=1,4
c     mu=1
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &    Phi(ithird,ix,iy,it,idirac,iflavor)
     &    -akappa*(R(ithird,mod(ix,ksize)+1,iy,it,idirac,iflavor)
     &    +R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,idirac,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=2
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &    -akappa*(R(ithird,ix,mod(iy,ksize)+1,it,idirac,iflavor)
     &    +R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,idirac,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=0
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &    -akappa*(R(ithird,ix,iy,mod(it,ksizet)+1,idirac,iflavor)
     &    +R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,idirac,iflavor))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
      do iflavor=1,2
      do idirac=1,4
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,ksizet,idirac,iflavor)=
     &   Phi(ithird,ix,iy,ksizet,idirac,iflavor)
     &    +R(ithird,ix,iy,1,idirac,iflavor)
      enddo
      enddo
      enddo
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,1,idirac,iflavor)=
     &   Phi(ithird,ix,iy,1,idirac,iflavor)
     &    +R(ithird,ix,iy,ksizet,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c
c     Dirac term
      do iflavor=1,2
      do idirac=1,4
c     mu=1
      mu=1
      igork=gamin(mu,idirac)
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &  +gamval(mu,idirac)*
     &    ( R(ithird,mod(ix,ksize)+1,iy,it,igork,iflavor)
     &    -R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,igork,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=2
      mu=2
      igork=gamin(mu,idirac)
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &  +gamval(mu,idirac)*
     &    ( R(ithird,ix,mod(iy,ksize)+1,it,igork,iflavor)
     &    -R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,igork,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=0
      mu=3
      igork=gamin(mu,idirac)
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &  +gamval(mu,idirac)*
     &    ( R(ithird,ix,iy,mod(it,ksizet)+1,igork,iflavor)
     &    -R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,igork,iflavor))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
      do iflavor=1,2
      do idirac=1,4
      igork=gamin(3,idirac)
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,ksizet,idirac,iflavor)=
     &   Phi(ithird,ix,iy,ksizet,idirac,iflavor)
     &    -2*gamval(3,idirac)*R(ithird,ix,iy,1,igork,iflavor)
      enddo
      enddo
      enddo
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,1,idirac,iflavor)=
     &   Phi(ithird,ix,iy,1,idirac,iflavor)
     &    +2*gamval(3,idirac)*R(ithird,ix,iy,ksizet,igork,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c
c  s-like term exploiting projection
      do iflavor=1,2
      do idirac=3,4
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &   -R(ithird+1,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      do idirac=1,2
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird+1,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird+1,ix,iy,it,idirac,iflavor)
     &    -R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c
c  coupling to chi (couples the two walls)
      if(imass.eq.1)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,3),0.0)
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
cc
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,1),-chi(ix,iy,it,2))
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      elseif(imass.eq.3)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(0.0,+chi(ix,iy,it,3))
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
cc
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,2),chi(ix,iy,it,1))
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      endif
c
      return
      end
c***********************************************************************
      subroutine dslashd(Phi,R,chi,am,imass)
c
c     calculates Phi = Mdagger*R
c
      parameter(ksize=12,ksizet=12,kthird=8)
      parameter(akappa=0.5)
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      real chi(ksize,ksize,ksizet,3)
      complex spm(ksize,ksize,ksizet)
      complex Phi(kthird,ksize,ksize,ksizet,4,2),
     &   R(kthird,ksize,ksize,ksizet,4,2)
      complex gamval
      integer gamin
c     write(6,*) 'hi from dslashd'
c
c     diagonal term
      diag=(3.0-am3)+1.0
      do iflavor=1,2
      do idirac=1,4
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   diag*R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c     Wilson term
      do iflavor=1,2
      do idirac=1,4
c     mu=1
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &    -akappa*(R(ithird,mod(ix,ksize)+1,iy,it,idirac,iflavor)
     &    +R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,idirac,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=2
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &    -akappa*(R(ithird,ix,mod(iy,ksize)+1,it,idirac,iflavor)
     &    +R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,idirac,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=0
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &    -akappa*(R(ithird,ix,iy,mod(it,ksizet)+1,idirac,iflavor)
     &   +R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,idirac,iflavor))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
      do iflavor=1,2
      do idirac=1,4
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,ksizet,idirac,iflavor)=
     &   Phi(ithird,ix,iy,ksizet,idirac,iflavor)
     &    +R(ithird,ix,iy,1,idirac,iflavor)
      enddo
      enddo
      enddo
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,1,idirac,iflavor)=
     &   Phi(ithird,ix,iy,1,idirac,iflavor)
     &    +R(ithird,ix,iy,ksizet,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c
c     Dirac term
      do iflavor=1,2
      do idirac=1,4
c     mu=1
      mu=1
      igork=gamin(mu,idirac)
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &-gamval(mu,idirac)*
     &    ( R(ithird,mod(ix,ksize)+1,iy,it,igork,iflavor)
     &     -R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,igork,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=2
      mu=2
      igork=gamin(mu,idirac)
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &-gamval(mu,idirac)*
     &    ( R(ithird,ix,mod(iy,ksize)+1,it,igork,iflavor)
     &     -R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,igork,iflavor))
      enddo
      enddo
      enddo
      enddo
c     mu=0
      mu=3
      igork=gamin(mu,idirac)
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &-gamval(mu,idirac)*
     &    ( R(ithird,ix,iy,mod(it,ksizet)+1,igork,iflavor)
     &     -R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,igork,iflavor))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
      do iflavor=1,2
      do idirac=1,4
      igork=gamin(3,idirac)
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,ksizet,idirac,iflavor)=
     &   Phi(ithird,ix,iy,ksizet,idirac,iflavor)
     &    +2*gamval(3,idirac)*R(ithird,ix,iy,1,igork,iflavor)
      enddo
      enddo
      enddo
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,1,idirac,iflavor)=
     &   Phi(ithird,ix,iy,1,idirac,iflavor)
     &    -2*gamval(3,idirac)*R(ithird,ix,iy,ksizet,igork,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c
c  s-like term exploiting projection
      do iflavor=1,2
      do idirac=1,2
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &   -R(ithird+1,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      do idirac=3,4
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird+1,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird+1,ix,iy,it,idirac,iflavor)
     &    -R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  chi term (couples the two walls)
      if(imass.eq.1)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,3),0.0)
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
cc
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,1),-chi(ix,iy,it,2))
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      elseif(imass.eq.3)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(0.0,-chi(ix,iy,it,3))
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(-chi(ix,iy,it,2),-chi(ix,iy,it,1))
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      endif
c
      return
      end
c***********************************************************************
      subroutine dslash_v2(Phi,R,chi,am,imass)
c
c     calculates Phi = M*R
c
      parameter(ksize=12,ksizet=12,kthird=8)
      parameter(akappa=0.5)
      common/para/beta,am3,ibound
      real chi(ksize,ksize,ksizet,3)
      complex spm(ksize,ksize,ksizet)
      complex Phi(kthird,ksize,ksize,ksizet,4,2),
     &   R(kthird,ksize,ksize,ksizet,4,2)
      complex zi,a,b
c     write(6,*) 'hi from dslash'
c
      zi=cmplx(0.0,1.0)
      zik=akappa*zi
c
c     diagonal term - Phi scaled by factor of two for now to save flops
      diag=2*((3.0-am3)+1.0)
      do iflavor=1,2
      do idirac=1,4
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &  diag*R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c     Wilson + Dirac term using specific gamma matrix representation in init
      do iflavor=1,2
c     mu=1
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,mod(ix,ksize)+1,iy,it,1,iflavor)+
     &   zi*R(ithird,mod(ix,ksize)+1,iy,it,4,iflavor)
      b=R(ithird,mod(ix,ksize)+1,iy,it,2,iflavor)+
     &   zi*R(ithird,mod(ix,ksize)+1,iy,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)+zi*b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)+zi*a
      enddo
      enddo
      enddo
      enddo
c
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,1,iflavor)-
     &   zi*R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,4,iflavor)
      b=R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,2,iflavor)-
     &   zi*R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)-zi*b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)-zi*a
c
      enddo
      enddo
      enddo
      enddo
c     mu=2
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,mod(iy,ksize)+1,it,1,iflavor)+
     &   R(ithird,ix,mod(iy,ksize)+1,it,4,iflavor)
      b=R(ithird,ix,mod(iy,ksize)+1,it,2,iflavor)-
     &   R(ithird,ix,mod(iy,ksize)+1,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)+b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)-a
      enddo
      enddo
      enddo
      enddo
c
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,1,iflavor)-
     &   R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,4,iflavor)
      b=R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,2,iflavor)+
     &   R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)-b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)+a
      enddo
      enddo
      enddo
      enddo
c
c     mu=0
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,iy,mod(it,ksizet)+1,1,iflavor)+
     &   zi*R(ithird,ix,iy,mod(it,ksizet)+1,3,iflavor)
      b=R(ithird,ix,iy,mod(it,ksizet)+1,2,iflavor)-
     &   zi*R(ithird,ix,iy,mod(it,ksizet)+1,4,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)+zi*a
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)-zi*b
      enddo
      enddo
      enddo
      enddo
c
c
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      a=R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,1,iflavor)-
     &   zi*R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,3,iflavor)
      b=R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,2,iflavor)+
     &   zi*R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,4,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)-zi*a
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)+zi*b
c
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  Now define Phi with correct normalisation
      do iflavor=1,2
      do idirac=1,4
      do it=1,ksizet
      do iy=1,ksize
      do ix = 1,ksize
      do ithird=1,kthird
         Phi(ithird,ix,iy,it,idirac,iflavor)=
     &    akappa*Phi(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
      do iflavor=1,2
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,iy,1,1,iflavor)+
     &   zi*R(ithird,ix,iy,1,3,iflavor)
      b=R(ithird,ix,iy,1,2,iflavor)-
     &   zi*R(ithird,ix,iy,1,4,iflavor)
      Phi(ithird,ix,iy,ksizet,1,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,1,iflavor)+a
      Phi(ithird,ix,iy,ksizet,2,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,2,iflavor)+b
      Phi(ithird,ix,iy,ksizet,3,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,3,iflavor)-zi*a
      Phi(ithird,ix,iy,ksizet,4,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,4,iflavor)+zi*b
c
      enddo
      enddo
      enddo
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,iy,ksizet,1,iflavor)-
     &   zi*R(ithird,ix,iy,ksizet,3,iflavor)
      b=R(ithird,ix,iy,ksizet,2,iflavor)+
     &   zi*R(ithird,ix,iy,ksizet,4,iflavor)
      Phi(ithird,ix,iy,1,1,iflavor)=
     &    Phi(ithird,ix,iy,1,1,iflavor)+a
      Phi(ithird,ix,iy,1,2,iflavor)=
     &    Phi(ithird,ix,iy,1,2,iflavor)+b
      Phi(ithird,ix,iy,1,3,iflavor)=
     &    Phi(ithird,ix,iy,1,3,iflavor)+zi*a
      Phi(ithird,ix,iy,1,4,iflavor)=
     &    Phi(ithird,ix,iy,1,4,iflavor)-zi*b
c
      enddo
      enddo
      enddo
      enddo
      endif
c
c  s-like term exploiting projection
      do iflavor=1,2
      do idirac=3,4
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &   -R(ithird+1,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      do idirac=1,2
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird+1,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird+1,ix,iy,it,idirac,iflavor)
     &    -R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c
c  coupling to chi (couples the two walls)
      if(imass.eq.1)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,3),0.0)
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
cc
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,1),-chi(ix,iy,it,2))
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      elseif(imass.eq.3)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(0.0,+chi(ix,iy,it,3))
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
cc
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,2),chi(ix,iy,it,1))
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      endif
c
      return
      end
c***********************************************************************
      subroutine dslashd_v2(Phi,R,chi,am,imass)
c
c     calculates Phi = Mdagger*R
c
      parameter(ksize=12,ksizet=12,kthird=8)
      parameter(akappa=0.5)
      common/para/beta,am3,ibound
      real chi(ksize,ksize,ksizet,3)
      complex spm(ksize,ksize,ksizet)
      complex Phi(kthird,ksize,ksize,ksizet,4,2),
     &   R(kthird,ksize,ksize,ksizet,4,2)
      complex zi,a,b
c     write(6,*) 'hi from dslashd'
c
      zi = cmplx(0.0,1.0)
      zik=akappa*zi
c
c     diagonal term - Phi is rescaled by factor two for now
      diag=2*((3.0-am3)+1.0)
      do iflavor=1,2
      do idirac=1,4
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   diag*R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c     Wilson + Dirac term in one exploiting r=1 term
      do iflavor=1,2
c     mu=1
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,mod(ix,ksize)+1,iy,it,1,iflavor)-
     &   zi*R(ithird,mod(ix,ksize)+1,iy,it,4,iflavor)
      b=R(ithird,mod(ix,ksize)+1,iy,it,2,iflavor)-
     &   zi*R(ithird,mod(ix,ksize)+1,iy,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)-zi*b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)-zi*a
      enddo
      enddo
      enddo
      enddo
c
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      a=R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,1,iflavor)+
     &   zi*R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,4,iflavor)
      b=R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,2,iflavor)+
     &   zi*R(ithird,mod(ix+ksize-2,ksize)+1,iy,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)+zi*b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)+zi*a
c
      enddo
      enddo
      enddo
      enddo
c     mu=2
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,mod(iy,ksize)+1,it,1,iflavor)-
     &   R(ithird,ix,mod(iy,ksize)+1,it,4,iflavor)
      b=R(ithird,ix,mod(iy,ksize)+1,it,2,iflavor)+
     &   R(ithird,ix,mod(iy,ksize)+1,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)-b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)+a
      enddo
      enddo
      enddo
      enddo
c
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      a=R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,1,iflavor)+
     &   R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,4,iflavor)
      b=R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,2,iflavor)-
     &   R(ithird,ix,mod(iy+ksize-2,ksize)+1,it,3,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)+b
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)-a
c
      enddo
      enddo
      enddo
      enddo
c     mu=0
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,iy,mod(it,ksizet)+1,1,iflavor)-
     &   zi*R(ithird,ix,iy,mod(it,ksizet)+1,3,iflavor)
      b=R(ithird,ix,iy,mod(it,ksizet)+1,2,iflavor)+
     &   zi*R(ithird,ix,iy,mod(it,ksizet)+1,4,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)-zi*a
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)+zi*b
      enddo
      enddo
      enddo
      enddo
c
      do  it=1,ksizet
      do  iy=1,ksize
      do  ix=1,ksize
      do ithird=1,kthird
      a=R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,1,iflavor)+
     &   zi*R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,3,iflavor)
      b=R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,2,iflavor)-
     &   zi*R(ithird,ix,iy,mod(it+ksizet-2,ksizet)+1,4,iflavor)
      Phi(ithird,ix,iy,it,1,iflavor)=
     &    Phi(ithird,ix,iy,it,1,iflavor)-a
      Phi(ithird,ix,iy,it,2,iflavor)=
     &    Phi(ithird,ix,iy,it,2,iflavor)-b
      Phi(ithird,ix,iy,it,3,iflavor)=
     &    Phi(ithird,ix,iy,it,3,iflavor)+zi*a
      Phi(ithird,ix,iy,it,4,iflavor)=
     &    Phi(ithird,ix,iy,it,4,iflavor)-zi*b
c
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  Now define Phi with correct normalisation
      do iflavor=1,2
      do idirac=1,4
      do it=1,ksizet
      do iy=1,ksize
      do ix = 1,ksize
      do ithird=1,kthird
         Phi(ithird,ix,iy,it,idirac,iflavor)=
     &    akappa*Phi(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
      do iflavor=1,2
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,iy,1,1,iflavor)-
     &   zi*R(ithird,ix,iy,1,3,iflavor)
      b=R(ithird,ix,iy,1,2,iflavor)+
     &   zi*R(ithird,ix,iy,1,4,iflavor)
      Phi(ithird,ix,iy,ksizet,1,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,1,iflavor)+a
      Phi(ithird,ix,iy,ksizet,2,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,2,iflavor)+b
      Phi(ithird,ix,iy,ksizet,3,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,3,iflavor)+zi*a
      Phi(ithird,ix,iy,ksizet,4,iflavor)=
     &    Phi(ithird,ix,iy,ksizet,4,iflavor)-zi*b
c
      enddo
      enddo
      enddo
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird
c
      a=R(ithird,ix,iy,ksizet,1,iflavor)+
     &   zi*R(ithird,ix,iy,ksizet,3,iflavor)
      b=R(ithird,ix,iy,ksizet,2,iflavor)-
     &   zi*R(ithird,ix,iy,ksizet,4,iflavor)
      Phi(ithird,ix,iy,1,1,iflavor)=
     &    Phi(ithird,ix,iy,1,1,iflavor)+a
      Phi(ithird,ix,iy,1,2,iflavor)=
     &    Phi(ithird,ix,iy,1,2,iflavor)+b
      Phi(ithird,ix,iy,1,3,iflavor)=
     &    Phi(ithird,ix,iy,1,3,iflavor)-zi*a
      Phi(ithird,ix,iy,1,4,iflavor)=
     &    Phi(ithird,ix,iy,1,4,iflavor)+zi*b
c
      enddo
      enddo
      enddo
      enddo
      endif
c
c  s-like term exploiting projection
      do iflavor=1,2
      do idirac=1,2
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird,ix,iy,it,idirac,iflavor)
     &   -R(ithird+1,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      do idirac=3,4
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      do ithird=1,kthird-1
      Phi(ithird+1,ix,iy,it,idirac,iflavor)=
     &   Phi(ithird+1,ix,iy,it,idirac,iflavor)
     &    -R(ithird,ix,iy,it,idirac,iflavor)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  chi term (couples the two walls)
      if(imass.eq.1)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,3),0.0)
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
cc
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(chi(ix,iy,it,1),-chi(ix,iy,it,2))
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      elseif(imass.eq.3)then
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(0.0,-chi(ix,iy,it,3))
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,1)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,1)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         enddo
         enddo
         enddo
         enddo
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         spm(ix,iy,it)=cmplx(-chi(ix,iy,it,2),-chi(ix,iy,it,1))
         enddo
         enddo
         enddo
         do idirac=1,2
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(kthird,ix,iy,it,idirac,1)=Phi(kthird,ix,iy,it,idirac,1)
     &    +spm(ix,iy,it)*R(1,ix,iy,it,idirac,2)
         Phi(kthird,ix,iy,it,idirac,2)=Phi(kthird,ix,iy,it,idirac,2)
     &    -conjg(spm(ix,iy,it))*R(1,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
         do idirac=3,4
         do it=1,ksizet
         do iy=1,ksize
         do ix=1,ksize
         Phi(1,ix,iy,it,idirac,1)=Phi(1,ix,iy,it,idirac,1)
     &    -spm(ix,iy,it)*R(kthird,ix,iy,it,idirac,2)
         Phi(1,ix,iy,it,idirac,2)=Phi(1,ix,iy,it,idirac,2)
     &    +conjg(spm(ix,iy,it))*R(kthird,ix,iy,it,idirac,1)
         enddo
         enddo
         enddo
         enddo
      endif
c
      return
      end

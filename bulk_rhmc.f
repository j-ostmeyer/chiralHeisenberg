       program dwf3d
c*******************************************************************
c    Rational Hybrid Monte Carlo algorithm for bulk Thirring Model with Domain Wall
c         fermions
c
c    requires operation of QMR on complex vectors to determine
c    (Mdagger M)**-1  Phi 
c
c    requires input partial fraction coefficients from Remez algorithm
c    of Clarke & Kennedy
c
c    the "third" direction is actually indexed 4 in the code - sorry!
c      (should have used C I know)
c
c           { 1 - hermition mass term psibar psi
c    imass= { 3 - antih.    mass term i psibar gamma_3 psi
c           { 5 - antih.    mass term i psibar gamma_5 psi
c     linear combinations of above require code modification
c
c    code exploits fact that gamma_3 is diagonal in Dirac basis used - speeds up 
c    evolution slightly, and measurement by factor of two. 
c
c    Fermion expectation values are measured using a noisy estimator.
c    on the Wilson matrix, which has dimension 4*kvol*kthird
c    inversions done using congrad, and matrix multiplies with dslash, dslashd
c
c    Pauli-Villars determinant defined using a hermitian mass m_h=One
c
c    trajectory length is random with mean dt*iterl
c    The code runs for a fixed number iter2 of trajectories.
c
c    Phi: pseudofermion field 
c    am: bare fermion mass 
c    actiona: running average of total action
c
c                                               SJH February 2017
c*******************************************************************
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(respbp=0.000001,rescgg=0.000001,rescga=0.000000001)
      parameter(rescgm=0.000000001)
      parameter(itermax=1000)
      parameter(One=1.0)
      parameter(ndiag=25,ndiagg=12,Nf=1)
      common/gauge/theta(kvol,3),seed
      common/trial/ut(kvol,3),thetat(kvol,3),pp(kvol,3)
      common /para/beta,am3,ibound
      common/remez2/anum2(0:ndiag),aden2(ndiag),
     &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4/anum4(0:ndiag),aden4(ndiag),
     &              bnum4(0:ndiag),bden4(ndiag)
      common/remez2g/anum2g(0:ndiagg),aden2g(ndiagg),
     &              bnum2g(0:ndiagg),bden2g(ndiagg)
      common/remez4g/anum4g(0:ndiagg),aden4g(ndiagg),
     &              bnum4g(0:ndiagg),bden4g(ndiagg)
      common/gforce/dSdpi(kvol,3)
      common /neighb/id(kvol,3),iu(kvol,3)
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common/trans/tpi 
      common/dum1/ R(kthird,kvol,4),ps(kvol,2)
      common/vector/X1(kthird,kvol,4)
      common/ranseed/yran,idum
      common/v/v(97)
c     complex Phi(kthird,kvol,4,Nf),X0(kthird,kvol,4)
c     complex R,zi,qq,qbqb
c     complex u,ut,X1
c     complex a,b
      complex*16 Phi(kthird,kvol,4,Nf),X0(kthird,kvol,4)
      complex*16 R,zi,qq,qbqb
      complex*16 u,ut,X1
      complex*16 a,b
      real rano
      real*8 H0,H1,S0,S1,dH,dS,hg,hp
      real*8 seed
      real*8 anum2,aden2,bnum2,bden2
      real*8 anum4,aden4,bnum4,bden4
      real*8 anum2g,aden2g,bnum2g,bden2g
      real*8 anum4g,aden4g,bnum4g,bden4g
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
      open(unit=39,file='remez4g',status='unknown')
      if(iread.eq.1) then
      call sread
      endif
      read(25,*) dt,beta,am3,am,imass,iterl,iter2
      close(25)
c set a new seed by hand...
      if(iseed.ne.0)then
      seed=4139764973254.0
      endif
      write(7,*) 'seed: ', seed
      call rranset(seed)
      idum=-1
      y=rano(yran,idum)
c     write(6,*) 'rano: ', y,yran,idum
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
      read(39,*) anum4g(0)
      do i=1,ndiag
      read(36,*) anum2(i),aden2(i)
      read(37,*) anum4(i),aden4(i)
      enddo
      do i=1,ndiagg
      read(38,*) anum2g(i),aden2g(i)
      read(39,*) anum4g(i),aden4g(i)
      enddo
      read(36,*) bnum2(0)
      read(37,*) bnum4(0)
      read(38,*) bnum2g(0)
      read(39,*) bnum4g(0)
      do i=1,ndiag
      read(36,*) bnum2(i),bden2(i)
      read(37,*) bnum4(i),bden4(i)
      enddo
      do i=1,ndiagg
      read(38,*) bnum2g(i),bden2g(i)
      read(39,*) bnum4g(i),bden4g(i)
      enddo
c*******************************************************************
c     print heading
c*******************************************************************
      traj=iterl*dt
      proby=1.0/float(iterl)
c     write(6, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg,
c    & iter2,beta,am3,am,imass
      write(7, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg,
     & iter2,beta,am3,am,imass
9001  format(' ksize=',i3,' ksizet=',i3,/
     1 ,' kthird=',i3,/
     1 ,' Nf =',i3,/
     1 ,' time step: dt=',f6.4,' trajectory length=',f9.6,/
     1 ,' Remez ndiag: action =',i3' guidance=',i3,/
     1 ,' # trajectories=',i6,' beta=',f9.6,/
     1 ,' am3=',f6.4,' am=',f6.4/
     1 ,' imass=',i2)
c     write(6,9004) rescgg,rescga,respbp
      write(7,9004) rescgg,rescga,respbp
9004  format(' Stopping residuals: guidance: ',e9.4,' acceptance: ',
     &     e9.4,' estimator: ',e9.4)
c     write(6,9044) rescgm
      write(7,9044) rescgm
9044  format(' Stopping residuals: meson: ',e9.4)
      call rranget(seed)
c     write(6,*) 'seed: ', seed
      write(7,*) 'seed: ', seed
c*******************************************************************
c       initialize for averages
c*******************************************************************
      actiona=0.0
      vel2a=0.0
      pbpa=0.0
      ancg=0.0
      ancgh=0.0
      ancgf=0.0
      ancgpf=0.0
      ancgpv=0.0
      ancghpv=0.0
      ancgfpv=0.0
      ancgpfpv=0.0
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
      do 2007 mu=1,3
      do 2007 i=1,kvol
      thetat(i,mu)=theta(i,mu)
2007  continue
c
      call coef(ut,thetat)
c*******************************************************************
c  Pseudofermion fields: Phi = {MdaggerM(1)}^-1/4 * {MdaggerM(m)}^1/4 * R, where
c   R is gaussian
c*******************************************************************
      do ia=1,Nf
c
      do idirac=1,4
      do ithird=1,kthird
      call gauss0(ps)
      do i=1,kvol
      R(ithird,i,idirac)=cmplx(ps(i,1),ps(i,2))
      enddo
      enddo
      enddo
c
c  For now Phi = {MdaggerM}^0.25 * R
c
      call qmrherm(R,rescga,itercg,am,imass,anum4,aden4,ndiag,
     &             0,isweep,0)
      ancgpf=ancgpf+float(itercg)
      do idirac=1,4
      do i=1,kvol
      do ithird=1,kthird
      R(ithird,i,idirac)=X1(ithird,i,idirac)
      enddo
      enddo
      enddo
      call qmrherm(R,rescga,itercg,One,1,bnum4,bden4,ndiag,
     &             0,isweep,0)
      ancgpfpv=ancgpfpv+float(itercg)
c
      do idirac=1,4
      do i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac,ia)=X1(ithird,i,idirac)
      enddo
      enddo
      enddo
c
      enddo
c*******************************************************************
c     heatbath for p 
c*******************************************************************
c  for some occult reason this write statement is needed to ensure compatibility with earlier versions
c     write(6,*) yran,idum
c     write(98,*) yran,idum
      do mu=1,3
      call gaussp(ps)
      do i=1,kvol
      pp(i,mu)=ps(i,1)
      enddo
      enddo
c     write(6,*) yran,idum
c*******************************************************************
c  call to Hamiltonian
c     
      call hamilton(Phi,
     &      H0,hg,hp,S0,rescga,isweep,0,am,imass)
      if(isweep.eq.1) then
      action=S0/kvol
      gaction=hg/kvol
      paction=hp/kvol
      endif 
c     goto 501
c*******************************************************************
c      half-step forward for p
c*******************************************************************
      call force(Phi,rescgg,am,imass,isweep,0)
      d=dt*0.5
      do 2004 mu=1,3
      do 2004 i=1,kvol
      pp(i,mu)=pp(i,mu)-dSdpi(i,mu)*d
2004  continue
c*******************************************************************
c     start of main loop for classical time evolution
c*******************************************************************
      do 500 iter=1,itermax
c
c  step (i) st(t+dt)=st(t)+p(t+dt/2)*dt;
c
      d=dt
      do 2001 mu=1,3
      do 2001 i=1,kvol
      thetat(i,mu)=thetat(i,mu)+d*pp(i,mu)
2001  continue
c
c  step (ii)  p(t+3dt/2)=p(t+dt/2)-dSds(t+dt)*dt (1/2 step on last iteration)
c
      call coef(ut,thetat)
      call force(Phi,rescgg,am,imass,isweep,iter)
c
c test for end of random trajectory
c 
      ytest=rano(yran,idum)
      if(ytest.lt.proby)then
      d=dt*0.5
      do 2005 mu=1,3
      do 2005 i=1,kvol
      pp(i,mu)=pp(i,mu)-d*dSdpi(i,mu)
2005  continue
      itot=itot+iter 
      goto 501
      else
      d=dt
      do 3005 mu=1,3
      do 3005 i=1,kvol
      pp(i,mu)=pp(i,mu)-d*dSdpi(i,mu)
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
     &       H1,hg,hp,S1,rescga,isweep,-1,am,imass)
      dH=H0-H1
      dS=S0-S1
      write(98,*) dH,dS
      y=exp(dH)      
      yav=yav+y 
      yyav=yyav+y*y 
      if(dH.lt.0.0)then
      x=rano(yran,idum)
      if(x.gt.y)goto 600
      endif
c
c     step accepted: set s=st
c
      do 2006 mu=1,3
      do 2006 i=1,kvol
      theta(i,mu)=thetat(i,mu)
2006  continue      
      naccp=naccp+1
      action=S1/kvol
      gaction=hg/kvol
      paction=hp/kvol
600   continue
      write(11,*) isweep,gaction,paction
      actiona=actiona+action 
      vel2=0.0
      do 457 mu=1,3
      do 457 i=1,kvol
      vel2=vel2+pp(i,mu)*pp(i,mu)
457   continue
      vel2=vel2/(3*kvol)
      vel2a=vel2a+vel2
c
c     uncomment to disable measurements
c     goto 601
666   if((isweep/iprint)*iprint.eq.isweep)then
      do 2066 mu=1,3
      do 2066 i=1,kvol
      thetat(i,mu)=theta(i,mu)
2066  continue      
      call coef(ut,thetat)
      call measure(pbp,respbp,ancgm,am,imass)
c     call meson(rescgm,itercg,ancgm,am,imass)
      pbpa=pbpa+pbp
      ancgma=ancgma+ancgm
      ipbp=ipbp+1
c     write(11,*) pbp
c     write(6,*) isweep,':  ',pbp,ancgm
      endif
c
      if((isweep/icheck)*icheck.eq.isweep)then
      call rranget(seed)
      call swrite
      flush(100)
      flush(200)
c     flush(302)
c     flush(400)
c     flush(500)
c     flush(501)
      if(imass.ne.1)then
c     flush(401)
c     flush(402)
c     flush(403)
      endif
c     write(7,9023) seed
      endif
c
601   continue
c*******************************************************************
c     end of main loop
c*******************************************************************
      actiona=actiona/iter2 
      vel2a=vel2a/iter2 
      pbpa=pbpa/ipbp
      ancg=ancg/(Nf*itot)
      ancgh=ancgh/(2*Nf*iter2)
      ancgpf=ancgpf/(Nf*iter2)
      ancgpv=ancgpv/(Nf*itot)
      ancgf=ancgf/(Nf*itot)
      ancgfpv=ancgfpv/(Nf*itot)
      ancghpv=ancghpv/(2*Nf*iter2)
      ancgpfpv=ancgpfpv/(iter2*Nf)
      ancgma=ancgma/ipbp
      yav=yav/iter2
      yyav=yyav/iter2-yav*yav 
      yyav=sqrt(yyav/(iter2-1)) 
      atraj=dt*itot/iter2 
c*******************************************************************
c     print global averages
c*******************************************************************
c     write(6, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgpv,ancgh,ancghpv,ancgf,
c    & ancgfpv,ancgpf,ancgpfpv,pbpa,vel2a,actiona
      write(7, 9022) iter2,naccp,atraj,yav,yyav,
     & ancg,ancgpv,ancgh,ancghpv,ancgf,ancgfpv,ancgpf,ancgpfpv,
     & pbpa,ancgma,vel2a,actiona
9022  format(' averages for last ',i6,' trajectories',/ 
     & 1x,' # of acceptances: ',i6,' average trajectory length= ',f8.3/
     & 1x,' <exp-dH>=',e11.4,' +/-',e10.3/
     2 1x,' av. # QMR itr.'/
     & 1x,'     guidance: DWF  ',f9.3,'; PV  ',f9.3/
     & 1x,'   acceptance: DWF  ',f9.3,'; PV  ',f9.3/
     & 1x,'        force: DWF  ',f9.3,'; PV  ',f9.3/
     & 1x,'pseudofermion: DWF  ',f9.3,'; PV  ',f9.3/
     1 1x,' psibarpsi=',e11.3/
     2 1x,' av. # QMR itr.',f9.3//
     & 1x,' mean square velocity=',e10.3,'; action per site=',e10.3//)
      write(7, 9024)
      write(7, 9024)
9024  format(1x)
c
      close(11)
c
      if(iwrite.eq.1) then
      call rranget(seed)
      call swrite
      write(7,*) 'seed: ', seed
      endif
c
      stop
      end
c******************************************************************
c   calculate dSds for gauge fields at each intermediate time
c******************************************************************
      subroutine force(Phi,res1,am,imass,isweep,iter)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksize*ksize*ksizet)
      parameter(kvol2=ksize*ksize)
      parameter(kferm=4*kvol*kthird)
      parameter(ndiagg=12,ndiag=ndiagg)
      parameter(One=1.0)
      parameter(Nf=1)
      common/remez2g/anum2(0:ndiag),aden2(ndiag),
     &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4g/anum4(0:ndiag),aden4(ndiag),
     &              bnum4(0:ndiag),bden4(ndiag)
      common/phi0/Phi0(kferm,25)
      common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
      common/para/beta,am3,ibound
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common/gforce/dSdpi(kvol,3)
      common/vector/X1(kferm)
c     complex Phi(kferm,Nf),X2(kferm)
c     complex X1,u,Phi0
      complex*16 Phi(kferm,Nf),X2(kferm)
      complex*16 X1,u,Phi0
      real*8 anum2,aden2,bnum2,bden2
      real*8 anum4,aden4,bnum4,bden4
c
c     write(6,111)
111   format(' Hi from force')
c
      do mu=1,3
      do i=1,kvol
         dSdpi(i,mu)=0.0
      enddo
      enddo
c
c uncomment this line to quench the fermions!
c     return
c pseudofermion action is
c   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
c
      do ia=1,Nf
c
      do i=1,kferm
      X2(i)=Phi(i,ia)
      enddo
      call qmrherm(X2,res1,itercg,One,1,anum4,aden4,ndiag,
     &             1,isweep,iter)
      ancgpv=ancgpv+float(itercg)
      do i=1,kferm
      X2(i)=X1(i)
      enddo
c
      call qmrherm(X2,res1,itercg,am,imass,bnum2,bden2,ndiag,
     &             0,isweep,iter)
      ancg=ancg+float(itercg)
c     write(111,*) itercg
      do i=1,kferm
      X2(i)=X1(i)
      enddo
c
c  evaluates -X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
      call qmrherm(X2,res1,itercg,am,imass,anum2,aden2,ndiag,
     &             2,isweep,iter)
      ancgf=ancgf+float(itercg)
c     write(113,*) itercg
c  evaluates +2Re{Phidagger * d/dpi[{MdaggerM(1)}^1/4] * X2}
      call qmrherm(X2,res1,itercg,One,1,anum4,aden4,ndiag,
     &             3,isweep,iter)
      ancgfpv=ancgfpv+float(itercg)
c
      enddo
c
      if(ibound.eq.-1)then
      ioffset=(ksizet-1)*kvol2
      do i=ioffset+1,kvol
      dSdpi(i,3)=-dSdpi(i,3)
      enddo
      endif
c
      b=beta*Nf
      do mu=1,3
      do i=1,kvol
         dSdpi(i,mu)=dSdpi(i,mu)+b*theta(i,mu)
      enddo
      enddo
c
      return
      end
c******************************************************************
c   Evaluation of Hamiltonian function
c******************************************************************
      subroutine hamilton(Phi,
     &       h,hg,hp,s,res2,isweep,iflag,am,imass)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(kmom=3*kvol,kferm=4*kvol*kthird)
      parameter(ndiag=25)
      parameter(One=1.0)
      parameter(Nf=1)
      common/trial/u(kvol,3),theta(kmom),pp(kmom)
      common/remez2/anum2(0:ndiag),aden2(ndiag),
     &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4/anum4(0:ndiag),aden4(ndiag),
     &              bnum4(0:ndiag),bden4(ndiag)
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common /para/beta,am3,ibound
      common/vector/ X1(kferm)      
      common/dum1/R(kferm),ps(kvol,2)
c     complex Phi(kferm,Nf)
c     complex X1,R
c     complex u
      complex*16 Phi(kferm,Nf)
      complex*16 X1,R
      complex*16 u
      real*8 hp,hg,hf,h,s
      real*8 anum2,aden2,bnum2,bden2
      real*8 anum4,aden4,bnum4,bden4
c     write(6,111)
111   format(' Hi from hamilton')
c
      hp=0.0
      hg=0.0
      hf=0.0
c
      do 22 imom=1,kmom
      hp=hp+pp(imom)*pp(imom)
22    continue
c
      hp=hp*0.5
c
      do imom=1,kmom
      hg=hg+theta(imom)*theta(imom)
      enddo
      hg=0.5*Nf*beta*hg
      h=hg+hp
c 
c uncomment these lines to quench the fermions!
c     write(6,*) isweep,':  hg', hg,'   hp', hp,'   h',h
c     return
c         
c  pseudofermion action is
c   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
c
      do ia=1,Nf
c
      do i=1,kferm
      R(i)=Phi(i,ia)
      enddo
      call qmrherm(R,res2,itercg,One,1,anum4,aden4,ndiag,
     &             0,isweep,iflag)
      ancghpv=ancghpv+float(itercg)
      do i=1,kferm
      R(i)=X1(i)
      enddo
c
      call qmrherm(R,res2,itercg,am,imass,bnum2,bden2,ndiag,
     &             0,isweep,iflag)
      ancgh=ancgh+float(itercg)
c
      do 4 iferm=1,kferm
      hf=hf+conjg(R(iferm))*X1(iferm)
4     continue
c
      enddo
c
      h=hg+hp+hf
c     write(6,*) isweep,':  hg', hg,'   hp', hp,'   hf', hf,
c    &   '   h',h
      s=hg+hf
c
      return
      end              
c******************************************************************
c    multisolver matrix inversion via Lanczos technique
c  eg. Golub & van Loan "Matrix Computations" 9.3.1
c       solves (MdaggerM+diag)*x=Phi for ndiag different values of diag
c   iflag=0: simply evaluates X = {MdaggerM}^p * Phi
c   can be interchanged with congrad for p=-1
c   iflag=1: in addition updates Phi0 register needed for PV force term
c   iflag=2: evaluates DWF force term
c   iflag=3: evaluates PV force term
c*****************************************************************m
      subroutine qmrherm(Phi,res,itercg,am,imass,anum,aden,ndiag,
     &                   iflag,isweep,iter)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*kthird*kvol)
c     parameter(niterc=10*kferm)
      parameter(niterc=7500)
      common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
      common/para/bbb,am3,ibound
      common/vector/x(kferm)
      common/gforce/dSdpi(kvol,3)
      common/phi0/Phi0(kferm,25)
c     complex Phi(kferm)
c     complex x,u,Phi0
      complex*16 Phi(kferm)
      complex*16 x,u,Phi0
      real*8 anum(0:ndiag),aden(ndiag),coeff
c     complex vtild(kferm),q(kferm),pm1(kferm,ndiag)
c     complex qm1(kferm),p(kferm,ndiag),x3(kferm),R(kferm)
c     complex x1(kferm,ndiag),x2(kferm)
c     real alpha(ndiag),beta
c     real amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
      complex*16 vtild(kferm),q(kferm),pm1(kferm,ndiag)
      complex*16 qm1(kferm),p(kferm,ndiag),x3(kferm),R(kferm)
      complex*16 x1(kferm,ndiag),x2(kferm)
      real*8 alpha(ndiag),beta
      real*8 amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
c
c     write(6,111)
111   format(' Hi from qmrherm')
c
      resid=sqrt(kferm*res*res)
c     write(6,*) iflag, resid
      itercg=0
c
c   initialise r=Phi
c
      do i=1,kferm
         R(i)=Phi(i)
         qm1(i)=(0.0,0.0)
         x(i)=anum(0)*Phi(i)
      enddo
      beta=0.0
      do i=1,kferm
         beta=beta+conjg(r(i))*r(i)
      enddo
      beta=sqrt(beta)
      phimod=beta
c     write(6,*) '|| Phi || = ', phimod
c
      do niter=1,niterc
      itercg=itercg+1
c
c  Lanczos steps
c
      do i=1,kferm
         q(i)=R(i)/beta
      enddo
c
      call dslash(vtild,q,u,am,imass)
      call dslashd(x3,vtild,u,am,imass)
c
      alphatild=0.0
      do i=1,kferm
         alphatild=alphatild+conjg(q(i))*x3(i)
      enddo
c
      do i=1,kferm
         R(i)=x3(i)-alphatild*q(i)-beta*qm1(i)
         qm1(i)=q(i)
      enddo
c
      beta0=beta
      beta=0.0
      do i=1,kferm
         beta=beta+conjg(R(i))*R(i)
      enddo
      beta=sqrt(beta)
c
      do idiag=1,ndiag
      alpha(idiag)=alphatild+aden(idiag)
      enddo
c
      if(niter.eq.1)then
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
c  criterion is a bit ad hoc -- relaxing by a factor arelax improves code
c  stability and leads to quicker convergence
           arelax=2.0
           if(rhomax.lt.arelax*resid) then
c          if(rhomax.lt.resid) then
c          call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiag)
c  convergence based on || residual || not working well in single precision...
c          if(resmax.lt.resid) goto 8
             goto 8
           endif
      endif
c
c   end of loop over iter
      enddo
      write(7,*) 'QMRniterc!, isweep,iter,iflag,imass,anum,ndiag = '
     &        ,isweep,iter,iflag,imass,anum(0),ndiag
8     continue
c
      if(iflag.lt.2)then
c  Now evaluate solution x=(MdaggerM)^p * Phi
      do idiag=1,ndiag
      do i=1,kferm
      x(i)=x(i)+anum(idiag)*x1(i,idiag)
      enddo
      enddo
c
c  update phi0 block if required...
      if(iflag.eq.1) then
      do idiag=1,ndiag
      do i=1,kferm
      Phi0(i,idiag)=x1(i,idiag)
      enddo
      enddo
      endif
c
      else
c
      do idiag=1,ndiag
c
c  X2 = M*X1
      do i=1,kferm
      R(i)=x1(i,idiag)
      enddo
      call dslash(X2,R,u,am,imass)
c
      if(iflag.eq.2)then
          coeff=anum(idiag)
          call derivs(R,X2,coeff,0)
      else
          coeff=-anum(idiag)
          do i=1,kferm
          R(i)=Phi0(i,idiag)
          enddo
          call derivs(R,X2,coeff,0)
c
          call dslash(X2,R,u,am,imass)
          do i=1,kferm
          R(i)=x1(i,idiag)
          enddo
          call derivs(X2,R,coeff,1)
      endif
c
      enddo
      endif
c
c
      return
      end
c**********************************************************************
c  iflag = 0 : evaluates Rdagger*(Mdagger)'*X2
c  iflag = 1 : evaluates Rdagger*(M)'*X2
c**********************************************************************
      subroutine derivs(R,X2,anum,iflag)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(akappa=0.5)
      common/neighb/id(kvol,3),iu(kvol,3)
      common/dirac/gamval(6,4),gamin(6,4)
      common/gforce/dSdpi(kvol,3)
c     complex R(kthird,kvol,4),X2(kthird,kvol,4)
      complex*16 R(kthird,kvol,4),X2(kthird,kvol,4)
c     complex gamval
      complex*16 gamval
      real*8 anum
      complex*16 tzi
      integer gamin
c
c     write(6,111)
111   format(' Hi from derivs')

c     dSdpi=dSdpi-Re(Rdagger *(d(Mdagger)dp)* X2)
c     Cf. Montvay & Muenster (7.215)
      tzi=cmplx(0.0,2*anum)
c     factor of 2 picks up second term in M&M (7.215)
c
      do mu=1,3
      do idirac=1,4
      do ithird=1,kthird
c
      do i=1,kvol
      dSdpi(i,mu)=dSdpi(i,mu)-akappa*real(tzi*
     &(conjg(R(ithird,i,idirac))*
     & X2(ithird,iu(i,mu),idirac)
     &-conjg(R(ithird,iu(i,mu),idirac))*
     &  X2(ithird,i,idirac)))
      enddo
c
      igork1=gamin(mu,idirac)
      if(iflag.eq.0)then
      do i=1,kvol
      dSdpi(i,mu)=dSdpi(i,mu)-real(tzi*gamval(mu,idirac)*
     &(conjg(R(ithird,i,idirac))*
     &        X2(ithird, iu(i,mu),igork1)
     &+conjg(R(ithird,iu(i,mu),idirac))*
     &             X2(ithird,i,igork1)))
      enddo
      else
      do i=1,kvol
      dSdpi(i,mu)=dSdpi(i,mu)+real(tzi*gamval(mu,idirac)*
     &(conjg(R(ithird,i,idirac))*
     &        X2(ithird,iu(i,mu),igork1)
     &+conjg(R(ithird,iu(i,mu),idirac))*
     &             X2(ithird,i,igork1)))
      enddo
      endif
c
      enddo
      enddo
      enddo
c
      return
      end
c******************************************************************
c   Calculates residual for testing purposes....
c   needs to run with double precision vectors to be useful.....
c******************************************************************
      subroutine testinv(Phi,resmax,itercg,am,imass,x,aden,ndiag)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*kthird*kvol)
      common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
      common/para/bbb,am3,ibound
      complex Phi(kferm)
      complex*16 vtild(kferm)
      complex*16 x3(kferm)
      complex*16 x(kferm,ndiag),x1(kferm),x2(kferm)
      complex*16 u
c     complex vtild(kferm)
c     complex x3(kferm)
c     complex x(kferm,ndiag),x1(kferm),x2(kferm)
c      complex u
      real*8 residual
      real*8 aden(ndiag)
c
      write(6,111)
111   format(' Hi from testinv')
c
      resmax=0.0
c
      do idiag=1,ndiag
      residual=0.0
      do i=1,kferm
      x3(i)=x(i,idiag)
      enddo
      call dslash(x2,x3,u,am,imass)
      call dslashd(x1,x2,u,am,imass)
      do i=1,kferm
      vtild(i)=x1(i)+aden(idiag)*x3(i)-Phi(i)
      residual=residual+conjg(vtild(i))*vtild(i)
      enddo
c     residual=sqrt(residual)
      if(residual.gt.resmax) resmax=residual
c
      write(6,*) idiag, 'itercg = ',itercg, ' residual = ',residual
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
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*kthird*kvol)
      parameter(niterc=kthird*kvol)
      common/trial/u(kvol,3),thetas(kvol,3),pp(kvol,3)
      common/para/bbb,am3,ibound
      common/vector/x(kferm)
c     complex Phi(kferm)
c     complex x,u
c     complex x1(kferm),x2(kferm),p(kferm),r(kferm)
      complex*16 Phi(kferm)
      complex*16 x,u
      complex*16 x1(kferm),x2(kferm),p(kferm),r(kferm)
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
      call dslash(x1,p,u,am,imass)
c
      if(nx.eq.1) goto 201
c
c   alpha=(r,r)/(p,(Mdagger)Mp)
c 
      alphad=0.0
      do 31 i=1,kferm
      alphad=alphad+conjg(x1(i))*x1(i)
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
      call dslashd(x2,x1,u,am,imass)
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
      betan=betan+conjg(r(i))*r(i) 
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
c     write(6,1000)
      write(7,1000)
1000  format(' # iterations of congrad exceeds niterc')
8     continue
c     write(6,*) itercg
      return
      end      
c*****************************************************************
c   Calculate fermion expectation values via a noisy estimator
c   -matrix inversion via conjugate gradient algorithm
c       solves Mx=x1
c     (Numerical Recipes section 2.10 pp.70-73)   
c*******************************************************************
      subroutine measure(psibarpsi,res,aviter,am,imass)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(akappa=0.5)
      parameter(knoise=10)
      common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      common /neighb/id(kvol,3),iu(kvol,3)
      common/vector/xi(kthird,kvol,4)
      common/ranseed/yran,idum
      common/v/v(97)
c     complex x(kvol,4), Phi(kthird,kvol,4)
c     complex xi,gamval
c     complex psibarpsi1,psibarpsi2
c     complex u
      complex*16 x(kvol,4), Phi(kthird,kvol,4)
      complex*16 xi,gamval
      complex*16 psibarpsi1,psibarpsi2
      complex*16 u
      real*8 cnum(0:1),cden(1)
      real ps(kvol,2),pt(kvol,2)
      real*8 pbp(knoise)
      integer gamin
c     write(6,*) 'hi from measure'
c
      iter=0
c     pbp=0.0
      cnum(0)=0.0
      cnum(1)=1.0
      cden(1)=0.0
c
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
         do idirac=1,4
           if(idirac.eq.idsource)then
             do i=1,kvol
               x(i,idirac)=cmplx(ps(i,1),ps(i,2))
             enddo
           else
             do i=1,kvol
               x(i,idirac)=(0.0,0.0)
             enddo
           endif
         enddo
      else
         do idirac=1,4
           if(idirac.eq.idsource+2)then
             do i=1,kvol
               x(i,idirac)=cmplx(ps(i,1),ps(i,2))
             enddo
           else
             do i=1,kvol
               x(i,idirac)=(0.0,0.0)
             enddo
           endif
         enddo
      endif
c
c
      if(imass.ne.5)then
         do ithird=1,kthird
         if(ithird.eq.1)then
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac)=x(i,idirac)
            enddo
            enddo
         else
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac)=(0.0,0.0)
            enddo
            enddo
         endif
         enddo
      else
c   xi = 0.5(1+gamma_4)*gamma_5*eta on DW at ithird=1
         do ithird=1,kthird
         if(ithird.eq.1)then
            do idirac=1,2
            do i=1,kvol
               xi(ithird,i,idirac)=-x(i,idirac+2)
               xi(ithird,i,idirac+2)=(0.0,0.0)
            enddo
            enddo
         else
c   xi = 0.5(1+gamma_4)*eta on DW at ithird=1
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac)=(0.0,0.0)
            enddo
            enddo
         endif
         enddo
      endif
c
c Phi= Mdagger*xi
c
      call dslashd(Phi,xi,u,am,imass)
c     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
      call congrad(Phi,res,itercg,am,imass)
      iter=iter+itercg
c
      if(imass.ne.5)then
c  pbp1 = x^dagger (0.5(1+gamma_4)) xi(kthird)
         do i=1,kvol
         psibarpsi1=psibarpsi1
     &    +conjg(x(i,idsource))*xi(kthird,i,idsource)
         enddo
      else
c  pbp1 = x^dagger (0.5(1-gamma_4)) xi(1)
         do i=1,kvol
         psibarpsi1=psibarpsi1
     &    +conjg(x(i,idsource+2))*xi(1,i,idsource+2)
         enddo
      endif
c
c
c  source on domain wall at ithird=kthird
      idsource2=idsource+2
c
      if(imass.ne.5)then
         do idirac=1,4
         do i=1,kvol
         if(idirac.eq.idsource2)then
             x(i,idirac)=cmplx(pt(i,1),pt(i,2))
         else
             x(i,idirac)=(0.0,0.0)
         endif
         enddo
         enddo
      else
         do idirac=1,4
         do i=1,kvol
         if(idirac.eq.idsource2-2)then
             x(i,idirac)=cmplx(pt(i,1),pt(i,2))
         else
             x(i,idirac)=(0.0,0.0)
         endif
         enddo
         enddo
      endif
c
      if(imass.ne.5)then
c   xi = 0.5(1-gamma_4)*eta on DW at ithird=kthird
         do ithird=1,kthird
         if(ithird.eq.kthird)then
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac)=x(i,idirac)
            enddo
            enddo
         else
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac)=(0.0,0.0)
            enddo
            enddo
         endif
         enddo
      else
c   xi = 0.5(1-gamma_4)*gamma_5*eta on DW at ithird=kthird
         do ithird=1,kthird
         if(ithird.eq.kthird)then
            do idirac=1,2
            do i=1,kvol
               xi(ithird,i,idirac+2)=-x(i,idirac)
               xi(ithird,i,idirac)=(0.0,0.0)
            enddo
            enddo
         else
            do idirac=1,4
            do i=1,kvol
               xi(ithird,i,idirac)=(0.0,0.0)
            enddo
            enddo
         endif
         enddo
      endif
c
c Phi= Mdagger*xi
c
      call dslashd(Phi,xi,u,am,imass)
c
c xi= (M)**-1 * Phi
c
c     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
      call congrad(Phi,res,itercg,am,imass)
      iter=iter+itercg
c
      if(imass.ne.5)then
c pbp2= - x^dagger (0.5(1-gamma_4)) xi(1)
         do i=1,kvol
         psibarpsi2=psibarpsi2
     &     +conjg(x(i,idsource2))*xi(1,i,idsource2)
         enddo
      else
c pbp2= - x^dagger (0.5(1-gamma_4)) xi(kthird)
         do i=1,kvol
         psibarpsi2=psibarpsi2
     &     +conjg(x(i,idsource))*xi(kthird,i,idsource)
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
c******************************************************************
c   Calculate meson correlators using point sources on domain walls
c   -matrix inversion via conjugate gradient algorithm
c       solves Mx=x1
c     (Numerical Recipes section 2.10 pp.70-73)   
c*******************************************************************
      subroutine meson(res,itercg,aviter,am,imass)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(ksize2=ksize*ksize)
      parameter(akappa=0.5)
      common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      common /neighb/id(kvol,3),iu(kvol,3)
      common/vector/xi(kthird,kvol,4)
      common/ranseed/yran,idum
      common/v/v(97)
c     complex x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
c     complex xi,gamval
c     complex prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
      complex*16 x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
      complex*16 xi,gamval
      complex*16 prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
c     complex prop00n(kvol,3:4,1:2),prop0Ln(kvol,3:4,3:4)
      real cpm(0:ksizet-1),cmm(0:ksizet-1)
c     complex cpmn(0:ksizet-1),cmmn(0:ksizet-1)
c     complex cferm1(0:ksizet-1), cferm2(0:ksizet-1)
c     complex u
      complex*16 cferm1(0:ksizet-1), cferm2(0:ksizet-1)
      complex*16 u
      real ps(kvol,2)
      real rano
      integer gamin
c     write(6,*) 'hi from meson'
c      
      nsource=5
      nsmear=10
      c=0.25
      iter=0
c
      do it=0,ksizet-1
      cpm(it)=(0.0,0.0)
      cmm(it)=(0.0,0.0)
c     cpmn(it)=(0.0,0.0)
c     cmmn(it)=(0.0,0.0)
      cferm1(it)=(0.0,0.0)
      cferm2(it)=(0.0,0.0)
      enddo
c
c      susceptibility
      chim=0.0
      chip=0.0
c
      do ksource=1,nsource
c
c   random location for +m source
      ixxx=int(ksize*rano(yran,idum))+1
      iyyy=int(ksize*rano(yran,idum))+1
      ittt=int(ksizet*rano(yran,idum))+1
      isource=ixxx+ksize*((iyyy-1)+ksize*(ittt-1))
c     write(6,*) ixxx,iyyy,ittt, isource
c
c  loop over Dirac index of source...   
      do idsource=3,4
c
c  source on domain wall at ithird=1
      do 300 ithird=1,kthird
      if(ithird.eq.1)then
      do idirac=1,4
      do i=1,kvol
      x(i,idirac)=(0.0,0.0)
      enddo
c  wall source
c     ioff=(ittt-1)*ksize2
c     do i=1,ksize2
c     x(i+ioff,idirac)=(1.0,0.0)/ksize2
c     enddo
c
      enddo
c  point source at fixed site, spin...
      x(isource,idsource)=(1.0,0.0)
c
c now smear it.....
c
      do ismear=1,nsmear
      call dslash2d(x0,x,u)
      do idirac=1,4
      do i=1,kvol
      x(i,idirac)=(1.0-c)*x(i,idirac)+c*x0(i,idirac)
      enddo
      enddo
      enddo
c
      else
      do idirac=1,4
      do i=1,kvol
      xi(ithird,i,idirac)=(0.0,0.0)
      enddo
      enddo
      endif
300   continue
c
c   xi = x  on DW at ithird=1
c
      do idirac=1,4
      do i=1,kvol
      xi(1,i,idirac)=x(i,idirac)
      enddo
      enddo
c
c Phi= Mdagger*xi
c
      call dslashd(Phi,xi,u,am,imass)
c
c  preconditioning (no,really)
      call dslashd(xi,Phi,u,am,imass)
c  
c xi= (MdaggerM)**-1 * Phi 
c
c     call congrad(Phi,res,itercg,am,imass) 
      iter=iter+itercg
c
      do idsink=1,2
      idsink2=idsink+2
      do i=1,kvol
      prop00(i,idsource,idsink)=xi(1,i,idsink)
      prop0L(i,idsource,idsink2)=xi(kthird,i,idsink2)
      enddo
      enddo
c
c     if(imass.ne.1)then
c  now evaluate with sign of mass reversed (not needed for hermitian mass term)
c     am=-am
c  source on domain wall at ithird=1
c     do 400 ithird=1,kthird
c     if(ithird.eq.1)then
c     do idirac=1,4
c     do i=1,kvol
c     xi(1,i,idirac)=x(i,idirac)
c     enddo
c     enddo
c     else
c     do idirac=1,4
c     do i=1,kvol
c     xi(ithird,i,idirac)=(0.0,0.0)
c     enddo
c     enddo
c     endif
400   continue
c
c Phi= Mdagger*xi
c
c     call dslashd(Phi,xi,u,am,imass)
c
c     call dslashd(xi,Phi,u,am,imass)
c  
c xi= (MdaggerM)**-1 * Phi 
c
c     call congrad(Phi,res,itercg,am,imass) 
c     iter=iter+itercg
c
c     do idsink=1,2
c     idsink2=idsink+2
c     do i=1,kvol
c     prop00n(i,idsource,idsink)=xi(1,i,idsink)
c     prop0Ln(i,idsource,idsink2)=xi(kthird,i,idsink2)
c     enddo
c     enddo
c
c     am=-am
c     endif
c
c  end loop on source Dirac index....
      enddo
c
c  Now tie up the ends....
c
c  First C+-
c
c  now evaluate the trace (exploiting projection)
      do id1=3,4
      do id2=1,2
      do it=0,ksizet-1
      itt=mod((ittt+it-1),ksizet)+1
      ioff=(itt-1)*ksize2
      do i=1,ksize2
      cpm(it)=cpm(it)
     &        +prop00(i+ioff,id1,id2)*conjg(prop00(i+ioff,id1,id2))
      enddo
      enddo
      enddo
      enddo
c
c     if(imass.ne.1)then
c     do id1=3,4
c     do id2=1,2
c     do it=0,ksizet-1
c     itt=mod((ittt+it-1),ksizet)+1
c     ioff=(itt-1)*ksize2
c     do i=1,ksize2
c     cpmn(it)=cpmn(it)
c    &        +prop00(i+ioff,id1,id2)*conjg(prop00n(i+ioff,id1,id2))
c     enddo
c     enddo
c     enddo
c     enddo
c     endif
c
c  next C--
c  now evaluate the trace exploiting projection
      do id1=3,4
      do id2=3,4
      do it=0,ksizet-1
      itt=mod((ittt+it-1),ksizet)+1
      ioff=(itt-1)*ksize2
      do i=1,ksize2
      cmm(it)=cmm(it)
     &   +prop0L(i+ioff,id1,id2)*conjg(prop0L(i+ioff,id1,id2))
      enddo
      enddo
      enddo
      enddo
c
c     if(imass.ne.1)then
c     do id1=3,4
c     do id2=3,4
c     do it=0,ksizet-1
c     itt=mod((ittt+it-1),ksizet)+1
c     ioff=(itt-1)*ksize2
c     do i=1,ksize2
c     cmmn(it)=cmmn(it)
c    &   +prop0L(i+ioff,id1,id2)*conjg(prop0Ln(i+ioff,id1,id2))
c     enddo
c     enddo
c     enddo
c     enddo
c     endif
c
c    now the fermion propagator
c  = tr{ P_-*Psi(0,1)Psibar(x,Ls) + gamma_0*P_-*Psi(0,1)Psibar(x,1) }
      do idd=3,4
      do it=0,ksizet-1
      itt=mod((ittt+it-1),ksizet)+1
c correct for apbc
      if(itt.ge.ittt)then
        isign=1
      else
        isign=ibound
      endif
c
      ioff=(itt-1)*ksize2
      do i=1,ksize2
      cferm1(it)=cferm1(it)
     & +isign*akappa*prop0L(i+ioff,idd,idd)
      cferm2(it)=cferm2(it)
     & +isign*gamval(3,idd)*prop00(i+ioff,idd,gamin(3,idd))
      enddo
      enddo
      enddo
c
c  finish loop over sources
      enddo
c
      do it=0,ksizet-1
      cpm(it)=cpm(it)/nsource
      cmm(it)=cmm(it)/nsource
c  Cf. (54) of 1507.07717
      chim=chim+2*(cpm(it)+cmm(it))
      enddo
c     if(imass.ne.1)then
c       if(imass.eq.3)then
c         do it=0,ksizet-1
c           cpmn(it)=cpmn(it)/nsource
c           cmmn(it)=cmmn(it)/nsource
c  Cf. (54),(61) of 1507.07717
c           chip=chip-2*(cpmn(it)-cmmn(it))
c         enddo
c       else
c         do it=0,ksizet-1
c           cpmn(it)=cpmn(it)/nsource
c           cmmn(it)=cmmn(it)/nsource
c  Cf. (64),(65) of 1507.07717
c           chip=chip-2*(cpm(it)-cmm(it))
c         enddo
c       endif
c     endif
c
      do it=0,ksizet-1
      write(302,*) it, cpm(it), cmm(it)
      write(500,*) it, real(cferm1(it)), aimag(cferm1(it))
      write(501,*) it, real(cferm2(it)), aimag(cferm2(it))
      enddo
c     write(6,*) chim
      write(400,*) chim
c     if(imass.ne.1)then
c     do it=0,ksizet-1
c     write(402,*) it, real(cpmn(it)), real(cmmn(it))
c     write(403,*) it, aimag(cpmn(it)), aimag(cmmn(it))
c     enddo
c     write(401,*) chip
c     endif
c
c     if(imass.eq.1)then
      aviter=float(iter)/(2*nsource)
c     else
c     aviter=float(iter)/(4*nsource)
c     endif
c
      return
      end
c*******************************************************************
c
      subroutine sread
      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
      common/gauge/ theta(kvol,3),seed
      real*8 seed
      open(unit=10,file='con',
     1     status='unknown',form='unformatted')
      read (10) theta,seed
      close(10)
      return
      end
c
      subroutine swrite
      parameter(ksize=12,ksizet=12,kvol=ksizet*ksize*ksize)
      common/gauge/ theta(kvol,3),seed
      real*8 seed
      open(unit=31,file='con',
     1     status='unknown',form='unformatted')
      write (31) theta,seed
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
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(akappa=0.5)
      common/gauge/theta(kvol,3),seed
      common /neighb/id(kvol,3),iu(kvol,3)
      common/dirac/gamval(6,4),gamin(6,4)
      common/ranseed/yran,idum
      common/v/v(97)
c     complex gamval,one,zi
      complex*16 gamval,one,zi
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
c     initialize gauge fields
c
      if(nc .eq. 1)goto 40
c     (else cold start)
      do 10 mu=1,3
      do 10 ind=1,kvol
      theta(ind,mu)=0.0
10    continue
      return
c
40    continue
      g=0.05
      do 61 mu=1,3
      do 61 ind=1,kvol
c     theta(ind,mu)=2.0*g*rranf()-1.0
      theta(ind,mu)=2.0*g*rano(yran,idum)-1.0
61    continue
      return
      end
c
      subroutine addrc
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
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
c******************************************************************
c   calculate compact links from non-compact links
c******************************************************************
      subroutine coef(u,theta)
      parameter(ksize=12,ksizet=12,kvol2=ksize*ksize,kvol=kvol2*ksizet)
      common/para/beta,am3,ibound
c     complex u(kvol,3)
      complex*16 u(kvol,3)
      real theta(kvol,3)
c
      do mu=1,3
      do i=1,kvol
c        u(i,mu)=exp(cmplx(0.0,theta(i,mu)))
         u(i,mu)=(1.0+cmplx(0.0,theta(i,mu)))
      enddo
      enddo
c
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
      ioffset=(ksizet-1)*kvol2
      do i=1,kvol2
         ind=ioffset+i
         u(ind,3)=-u(ind,3)
      enddo
      endif
c
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
                  dum=rranf()
               enddo
               do j=1,97
                  v(j)=rranf()
               enddo
               y=rranf()
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
          v(j)=rranf()
          return
          end
C========================================================================
C
          SUBROUTINE RRANGET(SEED)
          DOUBLE PRECISION    SEED,     G900GT,   G900ST,   DUMMY
          SEED  =  G900GT()
          RETURN
          ENTRY RRANSET(SEED)
          DUMMY  =  G900ST(SEED)
          RETURN
          END
          REAL FUNCTION RRANF()
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
             RRANF  =  SNGL(DR)
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
      subroutine dslash(Phi,R,u,am,imass)
c
c     calculates Phi = M*R
c
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(akappa=0.5)
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      common /neighb/id(kvol,3),iu(kvol,3)
c     complex u(kvol,3)
c     complex Phi(kthird,kvol,4),R(kthird,kvol,4)
c     complex gamval
c     complex zkappa
      complex*16 u(kvol,3)
      complex*16 Phi(kthird,kvol,4),R(kthird,kvol,4)
      complex*16 gamval
      complex*16 zkappa
      integer gamin
c     write(6,*) 'hi from dslash'
c
c     diagonal term
      diag=(3.0-am3)+1.0
      do idirac=1,4
      do  i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac)=diag*R(ithird,i,idirac)
      enddo
      enddo
      enddo
c
c     Wilson term
      do mu=1,3
      do idirac=1,4
      do i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac)=Phi(ithird,i,idirac)
     &    -akappa*(      u(i,mu)*R(ithird,iu(i,mu),idirac)
     &         +conjg(u(id(i,mu),mu))*R(ithird,id(i,mu),idirac))
      enddo
      enddo
      enddo
      enddo
c
c     Dirac term
      do mu=1,3
      do idirac=1,4
      igork=gamin(mu,idirac)
      do i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac)=Phi(ithird,i,idirac)
     &+gamval(mu,idirac)*
     &    (          u(i,mu)*R(ithird,iu(i,mu),igork)
     &         -conjg(u(id(i,mu),mu))*R(ithird,id(i,mu),igork))
      enddo
      enddo
      enddo
      enddo
c
c  s-like term exploiting projection
      do idirac=3,4
      do i=1,kvol
      do ithird=1,kthird-1
      Phi(ithird,i,idirac)=Phi(ithird,i,idirac)
     &   -R(ithird+1,i,idirac)
      enddo
      enddo
      enddo
      do idirac=1,2
      do ithird=1,kthird-1
      do i=1,kvol
      Phi(ithird+1,i,idirac)=Phi(ithird+1,i,idirac)
     &    -R(ithird,i,idirac)
      enddo
      enddo
      enddo
c
c  Mass term (couples the two walls unless imass=5)
      if(imass.eq.1)then
         zkappa=cmplx(am,0.0)
         do idirac=3,4
         do i=1,kvol
         Phi(kthird,i,idirac)=Phi(kthird,i,idirac)+zkappa*R(1,i,idirac)
         enddo
         enddo
         do idirac=1,2
         do i=1,kvol
         Phi(1,i,idirac)=Phi(1,i,idirac)+zkappa*R(kthird,i,idirac)
         enddo
         enddo
      elseif(imass.eq.3)then
         zkappa=cmplx(0.0,-am)
         do idirac=3,4
         do i=1,kvol
         Phi(kthird,i,idirac)=Phi(kthird,i,idirac)-zkappa*R(1,i,idirac)
         enddo
         enddo
         do idirac=1,2
         do i=1,kvol
         Phi(1,i,idirac)=Phi(1,i,idirac)+zkappa*R(kthird,i,idirac)
         enddo
         enddo
      elseif(imass.eq.5)then
         zkappa=cmplx(0.0,-am)
         do idirac=3,4
         igork=gamin(5,idirac)
         do i=1,kvol
         Phi(kthird,i,idirac)=Phi(kthird,i,idirac)
     &                        -zkappa*R(kthird,i,idirac-2)
c        Phi(kthird,i,idirac)=Phi(kthird,i,idirac)
c    &           +2*zkappa*gamval(5,idirac)*R(kthird,i,igork)
         enddo
         enddo
         do idirac=1,2
         igork=gamin(5,idirac)
         do i=1,kvol
         Phi(1,i,idirac)=Phi(1,i,idirac)-zkappa*R(1,i,idirac+2)
c        Phi(1,i,idirac)=Phi(1,i,idirac)
c    &        +2*zkappa*gamval(5,idirac)*R(1,i,igork)
         enddo
         enddo
      endif
c
      return
      end
c***********************************************************************
      subroutine dslashd(Phi,R,u,am,imass)
c
c     calculates Phi = Mdagger*R
c
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(akappa=0.5)
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      common /neighb/id(kvol,3),iu(kvol,3)
c     complex u(kvol,3)
c     complex Phi(kthird,kvol,4),R(kthird,kvol,4)
c     complex gamval
c     complex zkappa
      complex*16 u(kvol,3)
      complex*16 Phi(kthird,kvol,4),R(kthird,kvol,4)
      complex*16 gamval
      complex*16 zkappa
      integer gamin
c     write(6,*) 'hi from dslashd'
c
c     diagonal term (hermitian)
      diag=(3.0-am3)+1.0
      do idirac=1,4
      do  i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac)=diag*R(ithird,i,idirac)
      enddo
      enddo
      enddo
c
c     Wilson term (hermitian)
      do mu=1,3
      do idirac=1,4
      do i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac)=Phi(ithird,i,idirac)
     &    -akappa*(      u(i,mu)*R(ithird,iu(i,mu),idirac)
     &         +conjg(u(id(i,mu),mu))*R(ithird,id(i,mu),idirac))
      enddo
      enddo
      enddo
      enddo
c
c     Dirac term (antihermitian)
      do mu=1,3
      do idirac=1,4
      igork=gamin(mu,idirac)
      do i=1,kvol
      do ithird=1,kthird
      Phi(ithird,i,idirac)=Phi(ithird,i,idirac)
     &-gamval(mu,idirac)*
     &    (          u(i,mu)*R(ithird,iu(i,mu),igork)
     &         -conjg(u(id(i,mu),mu))*R(ithird,id(i,mu),igork))
      enddo
      enddo
      enddo
      enddo
c
c  s-like term exploiting projection
      do idirac=1,2
      do i=1,kvol
      do ithird=1,kthird-1
      Phi(ithird,i,idirac)=Phi(ithird,i,idirac)
     &   -R(ithird+1,i,idirac)
      enddo
      enddo
      enddo
      do idirac=3,4
      do ithird=1,kthird-1
      do i=1,kvol
      Phi(ithird+1,i,idirac)=Phi(ithird+1,i,idirac)
     &   -R(ithird,i,idirac)
      enddo
      enddo
      enddo
c
c  Mass term (couples the two walls unless imass=5) 
      if(imass.eq.1)then
         zkappa=cmplx(am,0.0)
         do idirac=1,2
         do i=1,kvol
         Phi(kthird,i,idirac)=Phi(kthird,i,idirac)+zkappa*R(1,i,idirac)
         enddo
         enddo
         do idirac=3,4
         do i=1,kvol
         Phi(1,i,idirac)=Phi(1,i,idirac)+zkappa*R(kthird,i,idirac)
         enddo
         enddo
      elseif(imass.eq.3)then
         zkappa=cmplx(0.0,am)
         do idirac=1,2
         do i=1,kvol
         Phi(kthird,i,idirac)=Phi(kthird,i,idirac)+zkappa*R(1,i,idirac)
         enddo
         enddo
         do idirac=3,4
         do i=1,kvol
         Phi(1,i,idirac)=Phi(1,i,idirac)-zkappa*R(kthird,i,idirac)
         enddo
         enddo
      elseif(imass.eq.5)then
         zkappa=cmplx(0.0,am)
         do idirac=1,2
         do i=1,kvol
         Phi(kthird,i,idirac)=Phi(kthird,i,idirac)
     &                        -zkappa*R(kthird,i,idirac+2)
         enddo
         enddo
         do idirac=3,4
         do i=1,kvol
         Phi(1,i,idirac)=Phi(1,i,idirac)-zkappa*R(1,i,idirac-2)
         enddo
         enddo
      endif
c
      return
      end
c***********************************************************************
      subroutine dslash2d(Phi,R,u)
c
c     calculates Phi = M*R
c
      parameter(ksize=12,ksizet=12,kvol=ksize*ksize*ksizet)
      parameter(akappa=0.5)
      common/dirac/gamval(6,4),gamin(6,4)
      common /neighb/id(kvol,3),iu(kvol,3)
c     complex u(kvol,3)
c     complex Phi(kvol,4),R(kvol,4)
c     complex gamval
      complex*16 u(kvol,3)
      complex*16 Phi(kvol,4),R(kvol,4)
      complex*16 gamval
      integer gamin
c     write(6,*) 'hi from dslash2d'
c
c     diagonal term
      diag=2.0
      do idirac=1,4
      do  i=1,kvol
      Phi(i,idirac)=diag*R(i,idirac)
      enddo
      enddo
c
c     Wilson term
      do mu=1,2
      do idirac=1,4
      do i=1,kvol
      Phi(i,idirac)=Phi(i,idirac)
     &    -akappa*(      u(i,mu)*R(iu(i,mu),idirac)
     &         +conjg(u(id(i,mu),mu))*R(id(i,mu),idirac))
      enddo
      enddo
      enddo
c
c     Dirac term
      do mu=1,2
      do idirac=1,4
      igork=gamin(mu,idirac)
      do i=1,kvol
      Phi(i,idirac)=Phi(i,idirac)
     &+gamval(mu,idirac)*
     &    (          u(i,mu)*R(iu(i,mu),igork)
     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),igork))
      enddo
      enddo
      enddo
c
      return
      end

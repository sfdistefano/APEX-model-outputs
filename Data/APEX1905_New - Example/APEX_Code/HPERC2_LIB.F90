    MODULE HPERC2_lib
    ! This module solves, flux, and tri solve the Richards eqn (RE) and the
    ! advection-dispersion eqn (ADE) for water in a soil
    ! profile for a specified period. Subroutine solve solves the RE in any 
    ! surface pond to solve the ADE. The basic reference for the methods is:
    ! Ross, P.J. 2003. Modeling soil water and solute transport - fast, simplified
    ! numerical solutions. Agron. J. 95:1352-1361.
!
! Copyright P.J. Ross 2001-2010
! Read the file "Terms of Use.doc" distributed with this software
!
    !USE numbers
    use PARM
    !USE hypvgvg ! set to required props module
    IMPLICIT NONE
    PRIVATE
    SAVE
    PUBLIC::botbc ! bottom boundary conditions
    PUBLIC::h0max,qprecmax,hbot ! boundary params
    PUBLIC::dSmax,dSmaxr,dtmax,dtmin,dsmmax,nwsteps ! solution parameters
    PUBLIC::solve ! solution routine
    PUBLIC::params,fvars ! derived types
    PUBLIC::gf,hmin,hdry,hjoin,drypwr,S1,par ! soil water parameters
    PUBLIC::nKv ! information
    PUBLIC::allo,hyofS,hyofh,Sofh,hypar,weight,deallo ! routines !temp cdj added deallo
    PUBLIC::fbd,dis,isotype,isopar ! soil solute parameters
    PUBLIC::solpar,setiso ! routines
    TYPE params
      REAL(FRK)::the,thre,hg,m,n,Kxm,p
      REAL(FRK)::he,Ke,KSe,phie,phiSe
      REAL(FRK)::hr1,hrS1,hp,hpwr
      REAL(FRK)::a0,a1,a2,b1
    END TYPE params
    TYPE fvars
      INTEGER::isat
      REAL(FRK)::h,phi,phiS,K,KS
    END TYPE fvars
    TYPE(params),DIMENSION(:),TARGET,ALLOCATABLE::par
    TYPE rapointer
      REAL(FRK),DIMENSION(:),POINTER::p
    END TYPE rapointer
    TYPE(rapointer),ALLOCATABLE::isopar(:,:)
    TYPE qa ! for quadratic approx to K(h) and phi(h)
      REAL(RKQA)::dh1,dh2
      INTEGER(IKQA)::mh1,mh2
      INTEGER(IKQA),DIMENSION(:),POINTER::nh
      REAL(RKQA),DIMENSION(:),POINTER::K,phi
    END TYPE qa
    TYPE qparams ! quadratic params
    ! quad approx for K is Kv(i)+a1*z+a2*z**2 where Kv is array of stored K values
    ! and z=(h-x0)*rdx (0<=z<=1)
      INTEGER::i
      REAL(FRK)::x0,rdx,a1,a2
    END TYPE qparams
    TYPE(qa),DIMENSION(:),TARGET,ALLOCATABLE::fun
    TYPE dryparams ! for near-dry water retention and condy
      REAL(FRK)::Sjoin,Kjoin,phijoin
    END TYPE dryparams
    TYPE dryparamsptr
      TYPE(dryparams),POINTER::p ! to a set of params
    END TYPE dryparamsptr
    TYPE(dryparamsptr),DIMENSION(:),ALLOCATABLE::drypar ! ptr for each prop type
    
    INTEGER,PARAMETER::drypwr=2
    REAL(FRK),PARAMETER::hjoin=-1.5e4_frk,hdry=-1.0e7_frk
    REAL(FRK),PARAMETER::drycS=-6.502290170873973_frk ! ln(hjoin/hdry)
    REAL(FRK),PARAMETER::dryc1=(-hjoin)**(-drypwr),dryc2=(-hdry)**(-drypwr)
    REAL(FRK),PARAMETER::drycK=one/(dryc1-dryc2),dryc3=one/(one-drypwr)
    REAL(FRK),PARAMETER::dryc4=hjoin*(dryc3*dryc1-dryc2)
    REAL(FRK),PARAMETER::S1=0.99_frk
    REAL(FRK),PARAMETER::Rmac1=0.25_frk,hmac1=-4.0_frk,hmac2=-40.0_frk
    REAL(FRK),PARAMETER::rerrK=0.01_frk
    INTEGER::nKv=0
    REAL(FRK)::gf=1.0_frk,hmin=-1.0e12_frk
    CHARACTER(LEN=2),DIMENSION(:,:),ALLOCATABLE::isotype
    REAL(FRK),ALLOCATABLE::fbd(:),dis(:)

    REAL(FRK),PARAMETER::dSfac=1.25_frk,dpmaxr=0.5_frk,h0min=-0.02_frk, &
      Smax=1.001_frk,dh0max=0.01_frk
    CHARACTER(LEN=20)::botbc="free drainage"
    INTEGER::nwsteps=10
    REAL(FRK)::h0max=1.0e10_frk,qprecmax=1.0e10_frk,hbot=0.0_frk
    REAL(FRK)::dSmax=0.05_frk,dSmaxr=0.5_frk,dtmax=1.0e10_frk,dtmin=0.0_frk,dsmmax=1.0_frk
    REAL(FRK),PARAMETER::third=1.0_frk/3.0_frk,sixth=1.0_frk/6.0_frk
    INTEGER::nless,nitsi
    !
    ! Definitions of public entities and private parameters (see above for default
    ! values):
    ! botbc    - bottom boundary condn for water; "constant head", "free drainage",
    !            "seepage", or "zero flux". Constant head means that matric head h
    !            is specified. Free drainage means zero gradient of matric head,
    !            i.e. unit hydraulic gradient. Seepage means zero flux when the
    !            matric head is below zero and an upper limit of zero for the head.
    ! h0max    - max pond depth allowed before runoff.
    ! qprecmax - max precipitation (or water input) rate (cm/h) for use with ponded
    !            constant head infiltration. If qprec > qprecmax then actual input
    !            rate is taken to be equal to infiltration plus evaporation rates.
    ! hbot     - matric head at bottom of profile when botbc set to "constant head".
    ! dSmax    - max change in S (the "effective saturation") of any unsaturated
    !            layer to aim for each time step; controls time step size.
    ! dSmaxr   - maximum negative relative change in S each time step. This
    !            parameter helps avoid very small or negative S.
    ! dtmax    - max time step allowed.
    ! dsmmax   - max solute change per time step (see dSmax); user should set this
    !            according to solute units used. Units for different solutes can be
    !            scaled by the user (e.g. to an expected max of around 1.0).
    ! nwsteps  - the solute routine is called every nwsteps of the RE solution.
    ! dSfac    - a change in S of up to dSfac*dSmax is accepted.
    ! dpmaxr   - relative change in matric flux potential (MFP) phi that is
    !            accepted for convergence when finding head h at soil interfaces.
    ! h0min    - min (negative) value for surface pond when it empties.
    ! Smax     - max value for layer saturation to allow some overshoot.
    ! dh0max   - allowable overshoot when pond reaches max allowed depth.
    ! solve    - sub to call to solve RE and ADE.
    !
    ! params      - type for water parameters. Params the, thre (=the-thr), hg, m,
    !               n, Kx and eta are for the vG and BC functions. he, Ke, KSe,
    !               phie and phiSe are values of variables h, K, KS, phi and phiS at
    !               saturation (denoted by the "e" for "air entry"), needed by
    !               module flow (MF). hr1 and hrS1 are values of h and hS at the
    !               saturation S1, above which a quadratic h(S) is used. bab is a
    !               parameter for the imcomplete beta function and hpwr is a head
    !               below which the vG retention is approximately a power law, for
    !               use in getting the gravity flow conductivity weighting w.
    ! fvars        - type for water variables used by MF and returned by subroutine
    !               hyofS (except for isat, which is 0 for unsaturated layers and 1
    !               for saturated layers).
    ! gf          - gravity factor for flow direction (usually 1 for straight down).
    ! hmin        - minimum matric head h (used by MF).
    ! hdry        - head at zero water content.
    ! hjoin       - head at which optional dry soil logarithmic water retention
    !               function joins vG function.
    ! drypwr      - power in conductivity function used with logarithmic retention.
    ! S1          - saturation above which a quadratic h(S) is used.
    ! par(:)      - hydraulic property params for soil types (also used by MF).
    ! nKv         - no. of stored K values for quadratic approx of K and the MFP
    !               phi. There are also approx nKv/2 stored phi values.
    ! allo        - subroutine to allocate parameter storage.
    ! hyofS       - subroutine to get water variable from saturation S (where S<1).
    ! hyofh       - subroutine to get some water variables from h.
    ! Sofh        - subroutine to get S from h.
    ! hypar       - subroutine to set soil hydraulic params.
    ! weight      - subroutine to get gravity flow conductivity weight w.
    ! fbd(:)       - bulk densities for soil types (used by MF).
    ! dis(:)      - dispersivities for soil types (used by MF).
    ! isotype(:)  - adsorption isotherm code for soil types (used by MF).
    ! isopar(:,:) - adsorption isotherm params for soil types (used by MF).
    ! solpar      - subroutine to set soil solute params.
    ! setiso      - subroutine to set soil solute isotherm type and params.
    ! Definitions of private parameters:
    ! dryc%  - parameters for dry soil wrc and K.
    ! Rmac   - value of SvG R factor at hmac1.
    ! hmac1  - value where SvG exponential "macropore" K starts.
    ! hmac2  - value where SvG non-exp "macropore" K starts.
    ! rerrK  - relative error for quadratic approx of K and MFP phi from hmac1 to
    !          min(hmac2,hg). This can be increased to reduce storage but accuracy
    !          of K and phi will be affected. Reducing it will increase accuracy.

    CONTAINS
    ! The subroutines solve, flux, and tri solve the Richards eqn (RE) and the
    ! advection-dispersion eqn (ADE) for water in a soil
    ! profile for a specified period. Subroutine solve solves the RE in any 
    ! surface pond to solve the ADE. The basic reference for the methods is:
    ! Ross, P.J. 2003. Modeling soil water and solute transport - fast, simplified
    ! numerical solutions. Agron. J. 95:1352-1361.

    !***********************************************************************************
    !***********************************************************************************
    
    SUBROUTINE solve(ts,tfin,qprec,qevap,n,nsol,dx,jt,h0,S,evap,runoff,infil,drn, &
      prk,nsteps,SUCCESS,heads,wex,cin,c0,sm,soff,sinfil,sdrn,nssteps)
    !***********************************************************************************
    !***********************************************************************************

    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n,nsol,jt(n)
    REAL(FRK),INTENT(IN)::ts,tfin,qprec,qevap,dx(n)!,ET(n)
    INTEGER,INTENT(INOUT)::nsteps
    LOGICAL,INTENT(OUT),OPTIONAL::success !temp cdj
    REAL(FRK),INTENT(INOUT)::h0,S(n),evap,runoff,infil,drn
    REAL(FRK),INTENT(IN),OPTIONAL::cin(nsol)
    INTEGER,INTENT(INOUT),OPTIONAL::nssteps(nsol)
    REAL(FRK),INTENT(INOUT),OPTIONAL::wex(n),c0(nsol),sm(n,nsol),soff(nsol), &
      sinfil(nsol),sdrn(nsol)
    REAL(FRK),INTENT(OUT),OPTIONAL::heads(n)
    REAL(FRK),INTENT(OUT)::prk(1:n) !temp cdj
    REAL(FRK)::calcPrk(0:n) !temp cdj
    ! Solves the RE and, optionally, the ADE from time ts to tfin.
    ! Definitions of arguments:
    ! Required args:
    ! ts      - start time (h).
    ! tfin    - finish time.
    ! qprec   - precipitation (or water input) rate (fluxes are in cm/h).
    ! qevap   - potl evaporation rate from soil surface.
    ! n       - no. of soil layers.
    ! nsol    - no. of solutes.
    ! dx(1:n) - layer thicknesses.
    ! jt(1:n) - layer soil type nos.
    ! h0      - surface head, equal to depth of surface pond.
    ! S(1:n)  - degree of saturation ("effective satn") of layers.
    ! evap    - cumulative evaporation from soil surface (cm, not initialised).
    ! runoff  - cumulative runoff.
    ! infil   - cumulative net infiltration (time integral of flux across surface).
    ! drn     - cumulative net drainage (time integral of flux across bottom).
    ! nsteps  - cumulative no. of time steps for RE soln.
    ! Optional args:
    ! heads(1:n)      - matric heads h of layers at finish.
    ! qexsub          - subroutine to get layer water extraction rates (cm/h) by
    !                   plants. Note that there is no solute extraction and osmotic
    !                   effects due to solute are ignored. Arguments:
    !                   jt(1:n) - layer soil type nos; h(1:n) - layer matric heads;
    !                   qex(1:n) - layer extraction rates; qexh(1:n) - partial
    !                   derivs of qex wrt h.
    ! wex(1:n)        - cumulative water extraction from layers.
    ! cin(1:nsol)     - solute concns in water input (user's units/cc).
    ! c0(1:nsol)      - solute concns in surface pond.
    ! sm(1:n,1:nsol)  - solute (mass) concns in layers.
    ! soff(1:nsol)    - cumulative solute runoff (user's units).
    ! sinfil(1:nsol)  - cumulative solute infiltration.
    ! sdrn(1:nsol)    - cumulative solute drainage.
    ! nssteps(1:nsol) - cumulative no. of time steps for ADE soln.
    LOGICAL again,getq0,getqn,init,initpond,maxpond
    INTEGER::i,iflux,ih0,iok,itmp,j,ns,nsat,nsatlast,nsteps0
    REAL(FRK)::accel,dmax,dt,dwinfil,dwoff,fac,infili,Khmin1,Kmin1,phimin1,phip, &
      qpme,qprec1,rsig,rsigdt,sig,t,ti,win
    REAL(FRK),DIMENSION(1)::Sbot
    REAL(FRK),DIMENSION(n-1)::dz
    REAL(FRK),DIMENSION(n)::hint,phimin,qex,qexd
    REAL(FRK),DIMENSION(n)::thi,thf
    REAL(FRK),DIMENSION(0:n)::aa,bb,cc,dd,dy,ee,q,qya,qyb
    REAL(FRK),DIMENSION(nsol)::cav,sinfili
    REAL(FRK),DIMENSION(n,nsol)::c
    TYPE(fvars)::vtop,vbot
    TYPE(fvars)::vcall(1)
    TYPE(fvars),DIMENSION(n),TARGET::var
    TYPE(fvars),POINTER::v
    TYPE(params),POINTER::p
    ! The derived types params and fvars hold soil water parameters and variables.
    ! Parameter names often end in e, which loosely denotes "air entry", i.e.,
    ! values at h=he. While values of water content th and hydraulic conductivity K
    ! at h=he are equal to those at saturation, the derivs wrt S are nonzero. The
    ! MFP phi for h>he is given by phi=phie+Ke*(h-he). The saturation status of a
    ! layer is stored as 0 or 1 in isat since S may be >1 (because of previous
    ! overshoot) when a layer desaturates. Fluxes at the beginning of a time step
    ! and their partial derivs wrt S or phi of upper and lower layers or boundaries
    ! are stored in q, qya and qyb.
    calcPrk=0. !temp cdj
    j=jt(1); p=>par(j)
    phip=max(p%phie-p%he*p%Ke,1.00001_frk*p%phie) ! phi at h=0
    ! get K, Kh and phi at hmin (hmin is smallest h, stored in hyprops)
    call hyofh(hmin,j,Kmin1,Khmin1,phimin1)
    dz=half*(dx(1:n-1)+dx(2:n)) ! flow paths
    !----- set up for boundary conditions
      getq0=.true.
      getqn=.false.
      if (botbc=="constant head") then ! h at bottom bdry specified
        getqn=.true.
        j=jt(n); p=>par(j)
        if (hbot<p%he) then
          Sbot(1)=Sofh(hbot,j)
          call hyofS(Sbot,1,(/j/),vcall)
          vbot=vcall(1)
          vbot%isat=0
        else
          vbot=fvars(1,hbot,(hbot-p%he)*p%Ke+p%phie,zero,p%Ke,zero)
        end if
      end if
    !----- end set up for boundary conditions
    !----- initialise
      t=ts; nsteps0=nsteps; nsat=0
      ! initialise saturated regions
      var%isat=0
      where (S>=one)
        var%phi=par(jt)%phie; var%K=par(jt)%Ke
        var%isat=1
      end where
      if (nsol>0) then
        ! set solute info
        thi=par(jt)%the-par(jt)%thre*(one-S) ! initial th (note: thre=the-thr)
        ti=t; infili=infil; sinfili=sinfil
        if (h0>zero.and.count(c0/=cin)>0) then
          initpond=.true. ! initial pond with different solute concn
        else
          initpond=.false.
        end if
        c=zero ! temp storage for soln concns
      end if
    !----- end initialise
    !----- solve until tfin
      do while (t<tfin)
        !----- take next time step
          do iflux=1,2 ! sometimes need twice to adjust phi at satn
            if (nsteps==nsteps0.and.iflux==1) then
              init=.true. ! flag to initialise h at soil interfaces
            else
              init=.false.
            end if
            nsatlast=nsat ! for detecting onset of profile saturation
            nsat=sum(var%isat) ! no. of sat layers
            sig=half; IF (nsat/=0) sig=one ! time weighting sigma
            rsig=one/sig
            ! update variables
            if (iflux==1) call hyofS(S,n,jt,var) ! for layers where S<1
            ! phi is solution var at satn, so h calc from phi where S>=1
            where (S>=one) var%h=par(jt)%he+(var%phi-par(jt)%phie)/par(jt)%Ke
    !        do i=1,8
    !          if(s(i)>=0.98) print*,s(i) !temp cdj delme !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        end do
            !----- get fluxes and derivs
              ! get surface condition
              p=>par(jt(1))
              if (var(1)%phi<=phip.and.h0<=zero.and.nsat<n) then ! no ponding
                ns=1 ! start index for eqns
                vtop=fvars(0,hmin,phimin1,zero,Kmin1,zero) ! fvars at soil surface
              else ! ponding
                ns=0
                vtop=fvars(1,h0,(h0-p%he)*p%Ke+p%phie,zero,p%Ke,zero)
              end if
              ! get bottom bdry condn
              if (botbc=="seepage") then
                if (var(n)%h>-half*gf*dx(n)) then
                  getqn=.true.
                  p=>par(jt(n))
                  vbot=fvars(1,zero,(zero-p%he)*p%Ke+p%phie,zero,p%Ke,zero)
                else
                  getqn=.false.
                end if
              end if
              ! get fluxes
              call getfluxes(n,jt,dx,dz,vtop,vbot,var,hint,phimin,q,qya,qyb, &
                iflux,init,getq0,getqn,SUCCESS)
              !IF(success==.FALSE.) RETURN !temp cdj - exit solve subroutine
              ! adjust for top and bottom bdry condns
              qpme=qprec-qevap ! input rate at saturation
              qprec1=qprec ! may change qprec1 to maintain pond if required
              if (ns==1) then
                if (q(0)<qpme) then
                  q(0)=qpme; qyb(0)=zero
                end if
                ! correction 10/2/2010
                maxpond=.false.
              else
                if (h0>=h0max.and.qpme>q(0)) then
                  maxpond=.true.
                  ns=1
                else
                  maxpond=.false.
                  ! change qya(0) from dq/dphi (returned by getfluxes) to dq/dh
                  qya(0)=par(jt(1))%Ke*qya(0)
                end if
              end if
              if (botbc/="constant head") then
                select case (botbc)
                  case ("zero flux")
                    q(n)=zero
                    qya(n)=zero
                  case ("free drainage")
                    v=>var(n)
                    q(n)=gf*v%K
                    if (v%isat==0) then
                      qya(n)=gf*v%KS
                    else
                      qya(n)=zero
                    end if
                  case ("seepage")
                    if (var(n)%h<=-half*gf*dx(n)) then
                      q(n)=zero
                      qya(n)=zero
                    end if
                  case default
                    write (*,*) "solve: illegal bottom boundary condn"
                    PAUSE
                end select
              end if
              again=.false. ! flag for recalcn of fluxes
            !----- end get fluxes and derivs
            !----- estimate time step dt
              dmax=zero
              thf=zero ! use thf as temp storage
              where (var%isat==0) thf=abs(q(1:n)-q(0:n-1))/(par(jt)%thre*dx)
              dmax=maxval(thf) ! max derivative |dS/dt|
              if (dmax>zero) then
                dt=dSmax/dmax
                ! if pond going adjust dt
                if (h0>zero.and.(q(0)-qpme)*dt>h0) dt=(h0-half*h0min)/(q(0)-qpme)
              else ! steady state flow
                if (qpme>=q(n)) then
                  ! step to finish
                  dt=tfin-t
                else
                  dt=-(h0-half*h0min)/(qpme-q(n)) ! pond going so adjust dt
                end if
              end if
              if (dt>dtmax) dt=dtmax ! user's limit
              ! if initial step, improve phi where S>=1
              if (nsteps==nsteps0.and.nsat>0.and.iflux==1) then
                again=.true.
                dt=1.0e-20*(tfin-ts)
              end if
              if (nsat==n.and.nsatlast<n.and.iflux==1) then
                ! profile has just become saturated so adjust phi values
                again=.true.
                dt=1.0e-20*(tfin-ts)
              end if
              if (t+1.1*dt>tfin) then ! step to finish
                dt=tfin-t
                t=tfin
              else
                t=t+dt ! tentative update
              end if
            !----- end estimate time step dt
            !----- get and solve eqns
              rsigdt=one/(sig*dt)
              ! aa, bb, cc and dd hold coeffs and rhs of tridiag eqn set
              aa(ns+1:n)=qya(ns:n-1); cc(ns:n-1)=-qyb(ns:n-1)
              dd(1:n)=-(q(0:n-1)-q(1:n))*rsig
              iok=0 ! flag for time step test
              itmp=0 ! counter to abort if not getting solution
              do while (iok==0) ! keep reducing time step until all ok
                itmp=itmp+1
                accel=one-0.05_frk*min(10,max(0,itmp-4)) ! acceleration
                if (itmp>50) then
                  write (*,*) "solve: too many iterations of equation solution"
                  STOP
                end if
                if (ns<1) then
                  bb(0)=-qya(0)-rsigdt
                  dd(0)=-(qpme-q(0))*rsig
                end if
                where (var%isat==0) bb(1:n)=qyb(0:n-1)-qya(1:n)- &
                par(jt)%thre*dx*rsigdt
                where (var%isat/=0) bb(1:n)=qyb(0:n-1)-qya(1:n)
                call tri(ns,n,aa,bb,cc,dd,ee,dy)
                ! dy contains dS or, for sat layers, dphi values
                iok=1
                if (.not.again) then
                  ! check if time step ok, if not then set fac to make it less
                  iok=1
                  do i=1,n
                    if (var(i)%isat==0) then ! check change in S
                      if (abs(dy(i))>dSfac*dSmax) then
                        fac=max(half,accel*abs(dSmax/dy(i))); iok=0; exit
                      end if
                      if (-dy(i)>dSmaxr*S(i)) then
                        fac=max(half,accel*dSmaxr*S(i)/(-dSfac*dy(i))); iok=0; exit
                      end if
                      if (S(i)<one.and.S(i)+dy(i)>Smax) then
                        fac=accel*(half*(one+Smax)-S(i))/dy(i); iok=0; exit
                      end if
                      if (S(i)>=one.and.dy(i)>half*(Smax-one)) then
                        fac=0.25*(Smax-one)/dy(i); iok=0; exit
                      end if
                    end if
                  end do
                  if (iok==1.and.ns<1.and.h0<h0max.and.h0+dy(0)>h0max+dh0max) then
                    ! start of runoff
                    fac=(h0max+half*dh0max-h0)/dy(0); iok=0
                  end if
                  if (iok==1.and.ns<1.and.h0>zero.and.h0+dy(0)<h0min) then
                    ! pond going
                    fac=-(h0-half*h0min)/dy(0); iok=0
                  end if
                  if (iok==0) then ! reduce time step
                    t=t-dt; dt=fac*dt; t=t+dt; rsigdt=1./(sig*dt)
                    nless=nless+1 ! count step size reductions
                  end if
                  v=>var(1)
                  if (v%isat/=0.and.iflux==1.and.v%phi<phip.and. &
                    v%phi+dy(1)>phip) then
                    ! incipient ponding - adjust state of saturated regions
                    t=t-dt; dt=1.0e-20*(tfin-ts); rsigdt=1./(sig*dt)
                    again=.true.; iok=0
                  end if
                end if
              end do
            !----- end get and solve eqns
            !----- update unknowns
              ih0=0
              if (.not.again) then
                dwoff=zero
                if (ns<1) then
                  h0=h0+dy(0)
                  if (h0<zero.and.dy(0)<zero) ih0=1 ! pond gone
                  evap=evap+qevap*dt
                  ! note that fluxes required are q at sigma of time step
                  dwinfil=(q(0)+sig*(qya(0)*dy(0)+qyb(0)*dy(1)))*dt
                else
                  dwinfil=(q(0)+sig*qyb(0)*dy(1))*dt
                  if (maxpond) then
                    evap=evap+qevap*dt
                    if (qprec>qprecmax) then ! set input to maintain pond
                      qpme=(q(0)+sig*qyb(0)*dy(1))
                      qprec1=qpme+qevap
                      dwoff=zero
                    else
                      dwoff=qpme*dt-dwinfil
                    end if
                    runoff=runoff+dwoff
                  else
                    evap=evap+qprec1*dt-dwinfil
                  end if
                end if
                infil=infil+dwinfil
                if (nsol>0) then ! get surface solute balance
                  if (initpond) then ! pond concn /= cin
                    if (h0>zero) then
                      if (ns==1) dy(0)=zero ! if max pond depth
                      cav=((two*h0-dy(0))*c0+qprec1*dt*cin)/(two*h0+dwoff+dwinfil)
                      c0=two*cav-c0
                    else
                      cav=((h0-dy(0))*c0+qprec1*dt*cin)/(dwoff+dwinfil)
                      initpond=.false. ! pond gone
                      c0=cin ! for output if any pond at end
                    end if
                    soff=soff+dwoff*cav
                    sinfil=sinfil+dwinfil*cav
                  else
                    soff=soff+dwoff*cin
                    sinfil=sinfil+(qprec1*dt-dwoff)*cin
                  end if
                end if
                if (botbc=="constant head") then
                  drn=drn+(q(n)+sig*qya(n)*dy(n))*dt
                else
                  drn=drn+(q(n)+sig*qya(n)*dy(n))*dt
                end if
                !begin temp cdj
                calcPrk=calcPrk+(q+sig*qya*dy)*dt
                !end temp cdj
                if (present(wex)) then
                  where (var%isat==0) wex=wex+(qex+sig*qexd*dy(1:n))*dt
                end if
              end if
              do i=1,n
                j=jt(i); p=>par(j); v=>var(i)
                if (v%isat==0) then
                  if (.not.again) then
                    S(i)=S(i)+dy(i)
                    if (S(i)>one.and.dy(i)>zero) then ! saturation of layer
                      v%isat=1; v%K=p%Ke; v%phi=p%phie
                    end if
                  end if
                else
                  v%phi=v%phi+dy(i)
                  if (i==1.and.ih0/=0.and.v%phi>=p%phie) v%phi=0. ! pond gone
                  if (v%phi<p%phie) then ! desaturation of layer
                    v%isat=0; v%K=p%Ke; v%phi=p%phie
                    v%KS=p%KSe; v%phiS=p%phiSe
                  end if
                end if
              end do
            !----- end update unknowns
            if (.not.again) exit
          end do
          if (dt<=dtmin) then
            write (2,*) "solve: time step = ",dt
            PAUSE
          end if
        !----- end take next time step
        ! remove negative h0 (optional)
        if (h0<zero.and.var(1)%isat==0) then
          infil=infil+h0
          S(1)=S(1)+h0/(par(jt(1))%thre*dx(1)); h0=zero
        end if
        nsteps=nsteps+1
      end do
    !----- end solve until tfin
    ! get heads if required
    if (present(heads)) then
      call hyofS(S,n,jt,var)
      heads=var%h
      where (S>=one) heads=par(jt)%he+(var%phi-par(jt)%phie)/par(jt)%Ke
    end if
    !begin temp cdj
    DO j=1,n
      prk(j)=calcPrk(j)
    END DO
    !end temp cdj
    !thOut=thf
    !CONTAINS
    END SUBROUTINE solve
    SUBROUTINE getfluxes(n,jt,dx,dz,vtop,vbot,var,hint,phimin,q,qya,qyb, &
      iflux,init,getq0,getqn,SUCCESS)
    IMPLICIT NONE
    LOGICAL,INTENT(IN)::init,getq0,getqn
    INTEGER,INTENT(IN)::n,jt(n),iflux
    REAL(FRK),INTENT(IN)::dx(n),dz(n-1)
    TYPE(fvars),INTENT(IN)::vtop,vbot
    TYPE(fvars),TARGET,INTENT(IN)::var(n)
    REAL(FRK),INTENT(INOUT)::hint(n),phimin(n)
    REAL(FRK),INTENT(OUT)::q(0:n),qya(0:n),qyb(0:n)
    LOGICAL,INTENT(OUT)::success !temp cdj
    ! Gets fluxes q and partial derivs qya, qyb wrt S (if unsat) or phi (if sat).
    ! Fluxes at top and bottom of profile, and fluxes due to plant extraction of
    ! water are included.
    ! Definitions of arguments:
    ! n           - no. of soil layers.
    ! jt(1:n)     - layer soil type nos.
    ! dx(1:n)     - layer thicknesses.
    ! dz(1:n-1)   - distances between layer centres.
    ! vtop        - water fvars at soil surface.
    ! vbot        - water fvars at bottom of profile.
    ! var(1:n)    - water fvars at layer centres.
    ! hint(1:n)   - values of h at interfaces are stored sequentially in hint.
    ! phimin(1:n) - similarly for phi at hmin in layers above interfaces.
    ! q(0:n)      - fluxes; q(i), i=1,...,n-1 is flux from layer i to layer i+1.
    !               q(0) is surface flux and q(n) is flux at bottom of profile.
    ! qya(0:n)    - partial deriv of q(i), i=0,...,n, wrt the variable to be solved
    !               for (S, phi or h) at upper end of flow path.
    ! qyb(0:n)    - ditto for var at lower end.
    ! iflux       - if iflux/=1, get only fluxes involving sat layers.
    ! init        - true if hint and phimin to be initialised.
    ! getq0       - true if q(0) required.
    ! getqn       - true if q(n) required.
    LOGICAL flag,limit
    INTEGER::i,itmp,j,l,m
    REAL(FRK)::dphii1,dhi,h1,h2,hi,Khi1,Khi2,phii1,q2,qya2,qyb2,y,y1,y2
    TYPE(params),POINTER::p,pm
    TYPE(fvars)::vi1,vi2
    TYPE(fvars),POINTER::v,vp
    v=>var(1)
    if (iflux==1.or.v%isat/=0) then ! get top flux if required
      if (getq0) then
        call flux(jt(1),vtop,v,half*dx(1),q(0),qya(0),qyb(0))
      end if
    end if
    ! get other fluxes
    l=0

    do i=1,n-1
      j=jt(i); p=>par(j)
      v=>var(i); vp=>var(i+1)
      if (iflux==1.or.v%isat/=0.or.vp%isat/=0) then ! get flux
        if (j==jt(i+1)) then ! same soil type, no interface
          call flux(j,v,vp,dz(i),q(i),qya(i),qyb(i))
        else ! interface
          l=l+1; m=jt(i+1); pm=>par(m)
          if (init) then ! initialise
            call hyofh(hmin,j,vi1%K,Khi1,phimin(l)) ! get phi at hmin
            h1=v%h; h2=vp%h
            y1=v%K*dx(i+1); y2=vp%K*dx(i)
            ! equate fluxes (K constant) to get initial estimate of h at interface
            hint(l)=(y1*h1+y2*h2+half*gf*(v%K-vp%K)*dx(i)*dx(i+1))/(y1+y2)
          end if
          hi=hint(l)
          flag=.true.; itmp=0
          ! iterate to get hi at interface for equal fluxes using Newton's method
          ! get dphii1 at interface in upper layer, because of better linearity,
          ! then convert to dhi
          do while (flag)
            itmp=itmp+1
            if (itmp>100) then
              write (*,*) "getfluxes: too many iterations finding interface h"
              success=.FALSE. !temp cdj
              STOP
            end if
            if (hi<p%he) then
              vi1%isat=0
              call hyofh(hi,j,vi1%K,Khi1,phii1)
              vi1%KS=Khi1/vi1%K ! use dK/dphi, not dK/dS
            else
              vi1%isat=1
              vi1%K=p%Ke; phii1=p%phie+(hi-p%he)*p%Ke; vi1%KS=zero
            end if
            vi1%h=hi; vi1%phi=phii1; vi1%phiS=one ! use dphi/dphi not dphi/dS
            call flux(j,v,vi1,half*dx(i),q(i),qya(i),qyb(i))
            if (hi<pm%he) then
              vi2%isat=0
              call hyofh(hi,m,vi2%K,Khi2,vi2%phi)
              vi2%KS=Khi2/vi2%K ! dK/dphi
            else
              vi2%isat=1; vi2%K=pm%Ke; vi2%phi=pm%phie+(hi-pm%he)*pm%Ke
            end if
            vi2%h=hi; vi2%phiS=one ! dphi/dphi
            call flux(m,vi2,vp,half*dx(i+1),q2,qya2,qyb2)
            qya2=qya2*vi2%K/vi1%K ! partial deriv wrt phii1
            ! adjust for equal fluxes
            dphii1=-(q(i)-q2)/(qyb(i)-qya2)
            limit=.false.
            if (phii1+dphii1<=phimin(l)) then ! out of range
              limit=.true.; dphii1=-half*(phii1-phimin(l))
            end if
            phii1=phii1+dphii1
            dhi=dphii1/(vi1%K+half*vi1%KS*dphii1) ! 2nd order Pade approx
            if (-vi1%KS*dphii1>1.5*vi1%K) then ! use 1st order approx for dhi
              dhi=dphii1/vi1%K
            end if
            hi=hi+dhi
            IF(HI<-999990.)FLAG=.FALSE.
            ! check for convergence - dphi/(mean phi)<=dpmaxr
            if (limit.or.abs(dphii1/(phii1-half*dphii1))>dpmaxr) then
              nitsi=nitsi+1 ! accumulate no. of interface its
            else
              flag=.false.
            end if
          end do
          q(i)=q(i)+qyb(i)*dphii1
          hint(l)=hi
          ! adjust derivs
          y=1./(qya2-qyb(i))
          qya(i)=qya(i)*qya2*y; qyb(i)=-qyb2*qyb(i)*y
        end if
      end if
    end do
    v=>var(n)
    if (iflux==1.or.v%isat/=0) then ! get bottom flux if required
      if (getqn) then
        call flux(jt(n),v,vbot,half*dx(n),q(n),qya(n),qyb(n))
      end if
    end if
      END SUBROUTINE getfluxes

    SUBROUTINE flux(j,v1,v2,dz,q,qya,qyb)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::j
    REAL(FRK),INTENT(IN)::dz
    REAL(FRK),INTENT(OUT)::q,qya,qyb
    TYPE(fvars),INTENT(IN)::v1,v2
    ! Gets flux and partial derivs for specified flow path.
    ! Definitions of arguments:
    ! j   - soil type no.
    ! v1  - water fvars at upper end of path.
    ! v2  - ditto at lower end.
    ! dz  - length of path.
    ! q   - flux.
    ! qya - partial deriv of flux wrt S (if unsat) or phi (if sat) at upper end.
    ! qyb - ditto at lower end.
    REAL(FRK)::w,rdz
    TYPE(params),POINTER::p
    ! gf is gravity factor (0 to 1) assumed available in module
    p=>par(j)
      if (gf<zero) then
        if ((v1%isat/=0.and.v2%isat/=0).or.v1%h-gf*(-dz)>=p%he) then
          ! correction 21/5/07
          !w=zero
          w=one
        else
          w=weight(j,v1%h,v1%K,v1%phi,-dz)
          w=one-w
        end if
      else
        if ((v1%isat/=0.and.v2%isat/=0).or.v2%h-gf*dz>=p%he) then
          w=zero
        else
          w=weight(j,v2%h,v2%K,v2%phi,dz)
        end if
      end if
    rdz=one/dz
    q=(v1%phi-v2%phi)*rdz+gf*(w*v1%K+(one-w)*v2%K)
    if (v1%isat==0) then
      qya=v1%phiS*rdz+gf*w*v1%KS
    else
      qya=rdz
    end if
    if (v2%isat==0) then
      qyb=-v2%phiS*rdz+gf*(1.-w)*v2%KS
    else
      qyb=-rdz
    end if
    END SUBROUTINE flux

    SUBROUTINE tri(ns,n,aa,bb,cc,dd,ee,dy)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ns,n
    REAL(FRK),DIMENSION(0:n),INTENT(IN)::aa,cc,dd
    REAL(FRK),DIMENSION(0:n),INTENT(INOUT)::bb,ee,dy
    ! Solves tridiag set of linear eqns. Coeff arrays aa and cc left intact.
    ! Definitions of arguments:
    ! ns      - start index for eqns.
    ! n       - end index.
    ! aa(0:n) - coeffs below diagonal; ns+1:n used.
    ! bb(0:n) - coeffs on diagonal; ns:n used.
    ! cc(0:n) - coeffs above diagonal; ns:n-1 used.
    ! dd(0:n) - rhs coeffs; ns:n used.
    ! ee(0:n) - work space.
    ! dy(0:n) - solution in ns:n.
    INTEGER::i
    dy(ns)=dd(ns) ! decomposition and forward substitution
    do i=ns,n-1
      ee(i)=cc(i)/bb(i)
      dy(i)=dy(i)/bb(i)
      bb(i+1)=bb(i+1)-aa(i+1)*ee(i)
      dy(i+1)=dd(i+1)-aa(i+1)*dy(i)
    end do
    dy(n)=dy(n)/bb(n) ! back substitution
    do i=n-1,ns,-1
      dy(i)=dy(i)-ee(i)*dy(i+1)
    end do
    END SUBROUTINE tri

! The following subroutines implement Schaap and van Genuchten (SvG) soil water retention and
! conductivity functions. It includes "macropore" conductivity components and an
! optional dry soil water retention component. Note that m=1-1/n for the vG wrc.

    SUBROUTINE allo(nt,ns)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nt,ns
    ! Allocate storage for soil water and solute parameters. This cannot be
    ! deallocated but it can be reused.
    ! Definitions of arguments:
    ! nt - no. of soil hydraulic property types.
    ! ns - no. of solutes.
    INTEGER::i,j
    allocate(par(nt),fun(nt),drypar(nt))
    allocate(isotype(nt,ns),fbd(nt),dis(nt),isopar(nt,ns))
    do i=1,nt
      nullify(fun(i)%nh,drypar(i)%p) ! so association can be tested
      do j=1,ns
        nullify(isopar(i,j)%p)
      end do
    end do
    END SUBROUTINE allo

    SUBROUTINE deallo()
    IMPLICIT NONE
    ! Allocate storage for soil water and solute parameters. This cannot be
    ! deallocated but it can be reused.
    ! Definitions of arguments:
    ! nt - no. of soil hydraulic property types.
    ! ns - no. of solutes.
    deallocate(par,fun,drypar)
    deallocate(isotype,fbd,dis,isopar)
    !do i=1,nt
    !  nullify(fun(i)%nh,drypar(i)%p) ! so association can be tested
    !  do j=1,ns
    !    nullify(isopar(i,j)%p)
    !  end do
    !end do
    END SUBROUTINE deallo

    SUBROUTINE hyofS(Sv,nl,jt,var)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nl,jt(:)
    REAL(FRK),INTENT(IN)::Sv(:)
    TYPE(fvars),DIMENSION(:),TARGET,INTENT(OUT)::var
    ! Get soil water variables from S.
    ! Definitions of arguments:
    ! Sv(1:nl)  - degree of saturation ("effective satn") of layers.
    ! nl        - no. of soil layers.
    ! jt(1:nl)  - layer soil type nos.
    ! var(1:nl) - other water fvars of layers.
    LOGICAL::SgeS1
    INTEGER::i,j
    REAL(FRK)::a1,a2,hr,hrS,lnS,r,rn,S,v1,v2,x,z
    TYPE(qparams)::qp
    TYPE(params),POINTER::p
    TYPE(dryparams),POINTER::dp
    TYPE(fvars),POINTER::v
    TYPE(qa),POINTER::f
    do i=1,nl
      if (Sv(i)<one) then
        j=jt(i); v=>var(i)
        S=Sv(i)
        IF(S<0.)CYCLE
        if (associated(drypar(j)%p)) then ! check for dry soil wrc
          dp=>drypar(j)%p
          if (S<dp%Sjoin) then ! use logarithmic dry soil wrc, K and phi
            a1=drycS/dp%Sjoin
            v%h=hdry*exp(a1*S)
            v1=(-v%h)**(-drypwr)
            v%K=drycK*dp%Kjoin*(v1-dryc2)
            v%KS=-drycK*dp%Kjoin*drypwr*v1*a1
            v%phi=dp%phijoin+drycK*dp%Kjoin*(v%h*(dryc3*v1-dryc2)-dryc4)
            v%phiS=v%K*v%h*a1
            cycle ! do next layer
          end if
        end if
        p=>par(j); f=>fun(j)
        ! get hr
        if (S>=S1) then ! quadratic h(S)
          SgeS1=.true.
          r=one/(S1-one)
          z=(S-one)*r
          a1=two*p%hr1-(S1-one)*p%hrS1
          a2=p%hr1-a1
          hr=z*(a1+z*a2)
          hrS=(a1+two*z*a2)*r
        else ! vG h(S)
          SgeS1=.false.
          rn=one/p%n
          lnS=log(S)
          x=exp(lnS/p%m)
          hr=exp(rn*log(one/x-one))
          hrS=-hr/(p%m*p%n*(one-x)*S)
        end if
        v%h=p%hg*hr
        ! get K and phi
        if (v%h>=hmac1) then ! SvG exponential "macropore" K and phi
          v%K=p%Ke*exp(-v%h/p%hp)
          v%KS=-v%K*p%hg*hrS/p%hp
          v%phi=f%phi(1)-p%hp*(v%K-f%K(1))
        elseif (v%h>hmac2+f%mh2*f%dh2) then ! SvG non-exp "macropore" K and phi
          call getqa(v%h,f,qp)
          z=(v%h-qp%x0)*qp%rdx
          v%K=f%K(qp%i)+z*(qp%a1+z*qp%a2)
          v%KS=(qp%a1+two*z*qp%a2)*qp%rdx*p%hg*hrS
          v%phi=f%phi(qp%i/2+1)+z*(f%K(qp%i)+z*(half*qp%a1+third*z*qp%a2))/qp%rdx
        else ! SvG modified vG K and phi
          if (SgeS1) then ! very unlikely, but get required functions of h
            x=one/(one+exp(p%n*log(hr)))
            lnS=p%m*log(x) ! actually ln(x**m), not ln(S)
            S=exp(lnS) ! actually x**m, not S
          end if
          v1=p%m*x*(two-x)/(two+(p%m-two)*x*(one+sixth*(p%m-one)*x))
          ! v1 is Pade approx for 1-(1-x)**m - saves ** and errors for very small x
          v2=exp(p%p*lnS)
          v%K=p%Kxm*v2*v1**2
          v%KS=v%K*(p%p+two*x*(one-v1)/((one-x)*v1))/S
          v1=(p%a0+x*(p%a1+x*p%a2))/(one+x*p%b1)
          ! v1 is Pade approx to series for getting phi
          v%phi=p%Kxm*p%hg*v2*S*x*v1
          v%phiS=v%K*p%hg*hrS
        end if
        v%phiS=v%K*p%hg*hrS
      end if
    end do
    END SUBROUTINE hyofS
    SUBROUTINE hyofh(h,j,K,Kh,phi)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::j
    REAL(FRK),INTENT(IN)::h
    REAL(FRK),INTENT(OUT)::K,Kh,phi
    ! Get soil water variables from h.
    ! Definitions of arguments:
    ! h   - matric head.
    ! j   - soil type no.
    ! K   - hydraulic conductivity.
    ! Kh  - derivative dK/dh.
    ! phi - matric flux potential (MFP).
    REAL(FRK)::hr,lnx,v1,v2,x,z
    TYPE(qparams)::qp
    TYPE(params),POINTER::p
    TYPE(qa),POINTER::f
    TYPE(dryparams),POINTER::dp
    if (associated(drypar(j)%p)) then
      dp=>drypar(j)%p
      if (h<hjoin) then ! use logarithmic dry soil K and phi
        x=(-h)**(-drypwr)
        K=drycK*dp%Kjoin*(x-dryc2)
        Kh=-drycK*dp%Kjoin*drypwr*x/h
        phi=dp%phijoin+drycK*dp%Kjoin*(h*(dryc3*x-dryc2)-dryc4)
        return ! finished
      end if
    end if
    p=>par(j); f=>fun(j)
    if (h>=hmac1) then ! SvG exponential "macropore" K and phi
      K=p%Ke*exp(-h/p%hp)
      Kh=-K/p%hp
      phi=f%phi(1)-p%hp*(K-f%K(1))
    elseif (h>hmac2+f%mh2*f%dh2) then ! SvG non-exp SvG "macropore" K and phi
      call getqa(h,f,qp)
      z=(h-qp%x0)*qp%rdx
      K=f%K(qp%i)+z*(qp%a1+z*qp%a2)
      Kh=(qp%a1+two*z*qp%a2)*qp%rdx
      phi=f%phi(qp%i/2+1)+z*(f%K(qp%i)+z*(half*qp%a1+third*z*qp%a2))/qp%rdx
    else ! SvG modified vG K and phi
      hr=h/p%hg
      x=one/(one+exp(p%n*log(hr)))
      v1=p%m*x*(two-x)/(two+(p%m-two)*x*(one+sixth*(p%m-one)*x))
      ! v1 is Pade approx for 1-(1-x)**m - saves ** and errors for very small x
      lnx=log(x)
      v2=exp(p%p*p%m*lnx)
      K=p%Kxm*v2*v1**2
      Kh=-p%m*p%n*K*(p%p*(one-x)+two*x*(one/v1-one))/(p%hg*hr)
      v1=(p%a0+x*(p%a1+x*p%a2))/(one+x*p%b1)
      ! v1 is Pade approx to series for getting phi
      phi=p%Kxm*p%hg*exp((p%m*(p%p+one)+one)*lnx)*v1
    end if
    END SUBROUTINE hyofh
    SUBROUTINE hypar(j,thr,the,hg,Kx,mn,pp,Ks,drywrc)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::j
    REAL(FRK),INTENT(IN)::thr,the,hg,Kx,mn,pp,Ks
    LOGICAL,INTENT(IN),OPTIONAL::drywrc
    ! Set soil hydraulic property parameters.
    ! Note that a call without drywrc deallocates dry wrc memory.
    ! Definitions of arguments:
    ! j      - soil type no.
    ! thr    - "residual" water content.
    ! the    - saturated water content.
    ! Kx     - saturated hydraulic conductivity of soil "matrix".
    ! mn     - vG retention function shape parameter.
    ! pp     - pore parameter for conductivity.
    ! Ks     - saturated hydraulic conductivity.
    ! drywrc - logarithmic wrc used in dry soil (h<hjoin) if drywrc true.
    INTEGER,PARAMETER::maxpts=9999999 ! for K(h) approx by quadratics - can increase
    INTEGER::i,jj,jt(1),kk,m,n,nh(maxpts)
    REAL(FRK),PARAMETER::rerrKw=0.01_frk ! rel error for power law for w approx
    ! can increase rerrKw for speed, but may affect accuracy
    REAL(FRK)::a,dh,K,Kh,phi,S(1),v1,v2,x
    REAL(FRK),DIMENSION(maxpts)::h,Kv,y
    TYPE(qparams)::qp
    TYPE(params),POINTER::p
    TYPE(fvars)::v(1)
    TYPE(qa),POINTER::f
    if (.not.allocated(par)) then
      write (*,*) "hypar: must call allo before hypar to allocate storage"
      stop
    end if
    if (size(par)<j) then
      write (*,*) "hypar: not enough space allocated for soil types"
      stop
    end if
    jt(1)=j; p=>par(j); f=>fun(j)
    p%n=one+mn; p%m=mn/p%n ! note that kvg must be one
    p%the=the; p%thre=the-thr ! note that the-thr is stored rather than thr
    p%hg=hg; p%he=-.001_frk !temp cdj changed from zero to -one
    p%p=pp; p%Ke=Ks
    p%hp=hmac1/((one-Rmac1)*log(Ks/Kx)) ! K=Ks*exp(-h/p%hp)
    x=exp(log(S1)/p%m)
    p%hr1=exp(log(one/x-one)/p%n) ! hr at S1
    p%hrS1=-p%hr1/(mn*S1*(one-x)) ! dhr/dS at S1
    x=one/(one+exp(p%n*log(hmac1/hg)))
    v1=exp(p%m*pp*log(x))
    v2=one-exp(p%m*log(one-x))
    p%Kxm=Kx/(v1*v2**2) ! v1*v2**2 is modified vG K function at hmac1
    p%hpwr=min(hmac2,hg*exp(-log(rerrKw/(p%m*(p%p+one)+one))/p%n)) ! for w approx
    ! get params for phi approx
    a=p%m*(pp+one)+one
    if (a<=0.001_frk) then
      write (*,*) "hypar: pp for K too small"
      stop
    end if
    v1=(p%m**2+11.0_frk)/(a+two)
    v2=(p%m-one)*p%m**2
    p%b1=-two*(p%m**2+5.0_frk)/((a+3.0_frk)*v1)
    p%a0=v2/a
    p%a1=v2*(one/(a+one)+p%b1/a)
    p%a2=v2*(p%b1/(a+one)+v1/12.0_frk)
    ! get quadratic approx for K and phi
    f%mh1=3+(hmac2-hmac1)/hg ! no. of intervals from hmac1 to hmac2
    f%dh1=(hmac2-hmac1)/f%mh1
    m=f%mh1+1 ! no. of h values
    do i=1,m-1
      h(i)=hmac1+(i-1)*f%dh1 ! h values
    end do
    h(m)=hmac2
    x=one/(one+exp(p%n*log(hmac2/hg)))
    if (x>0.55_frk) then ! get extra quads
      f%mh2=one+(hg-hmac2)/hg
      f%dh2=(hg-h(m))/f%mh2
      do i=m+1,m+f%mh2
        h(i)=h(m)+(i-m)*f%dh2
      end do
      m=m+f%mh2 ! new no. of h values
    else
      f%mh2=0 ! no extra quads
    end if
    jj=1; nh(1)=1
    do i=1,m-1
      call fitquad(h(i),h(i+1),rerrK,p,y,n)
    !  write (2,*) "n ",n
      if (jj+n-1>maxpts)EXIT
        !write (*,*) "hypar: too many points for K(h)"
        !stop
      !end if
      Kv(jj:jj+n-1)=y(1:n)
      jj=jj+n-1
      nh(i+1)=jj
    end do
    if (associated(f%nh)) then
      nKv=nKv-f%nh(size(f%nh)) ! adjust counter
      deallocate(f%nh,f%K,f%phi)
    end if
    nKv=nKv+jj ! cumulative no. of K values
    allocate(f%nh(m),f%K(jj),f%phi(jj/2+1))
    f%nh=nh(1:m); f%K=Kv(1:jj)
    !write (2,"(a5,10i5)") "nh ",f%nh
    !write (2,*) h(1:m)
    !write (2,"(8f10.4)") f%K
    ! get phi values from K by integrating quads
    if (f%mh2>0) x=0.5_frk
    v1=(p%a0+x*(p%a1+x*p%a2))/(one+x*p%b1)
    f%phi(jj/2+1)=p%Kxm*p%hg*exp(a*log(x))*v1
    do kk=m,2,-1 ! step through major intervals from last to first
      if (kk>f%mh1+1) then
        dh=f%dh2/(f%nh(kk)-f%nh(kk-1)) ! half quadratic interval
      else
        dh=f%dh1/(f%nh(kk)-f%nh(kk-1))
      end if
      jj=1
      do i=f%nh(kk)-2,f%nh(kk-1),-2 ! step through quadratics
        call getqa(h(kk)-jj*dh,f,qp)
        ! add integral of K over quad from lower to higher h
        f%phi(i/2+1)=f%phi(i/2+2)-(f%K(i)+half*qp%a1+third*qp%a2)/qp%rdx
        jj=jj+2
      end do
    end do
    !write (2,"(8f10.4)") f%phi
    if (present(drywrc)) then
      if (drywrc) then
        if (abs(thr)>epsilon(one)) then
          write (*,*) "hypar: thr must be zero with dry wrc"
          stop
        end if
        if (.not.associated(drypar(j)%p)) then
          allocate(drypar(j)%p)
        end if
        ! get params for join of wet S and dry S
        drypar(j)%p%Sjoin=Sofh(hjoin,j)
        call hyofh(hjoin,j,K,Kh,phi)
        drypar(j)%p%phijoin=phi
        drypar(j)%p%Kjoin=K
      end if
    else
      if (associated(drypar(j)%p)) then ! no param, deallocate memory
        deallocate(drypar(j)%p)
      end if
    end if
    !WRITE(51,*)ONE,EPSILON(ONE)
    S(1)=one-epsilon(one)
    call hyofS(S,1,jt,v)
    p%KSe=v(1)%KS ! dK/dS at he
    p%phie=v(1)%phi ! MFP at he
    p%phiSe=v(1)%phiS ! dphi/dS at he
    END SUBROUTINE hypar
    SUBROUTINE getqa(x,f,qp)
    IMPLICIT NONE
    REAL(FRK),INTENT(IN)::x
    TYPE(qa),INTENT(IN)::f
    TYPE(qparams),INTENT(OUT)::qp
    ! Get quadratic approx params.
    ! Definitions of arguments:
    ! x  - head h.
    ! f  - quad approx data.
    ! qp - quad params.
    INTEGER::i,j,n
    REAL(FRK)::dx0,dxj,x0,xj
    ! find position in array
    if (x<hmac2) then
      x0=hmac2; dx0=f%dh2; n=f%mh1
    else
      x0=hmac1; dx0=f%dh1; n=0
    end if
    j=one+(x-x0)/dx0
    xj=x0+(j-1)*dx0
    j=j+n
    dxj=dx0/((f%nh(j+1)-f%nh(j))/2)
    qp%rdx=one/dxj
    i=(x-xj)*qp%rdx
    qp%x0=xj+i*dxj
    i=f%nh(j)+2*i
    ! get quadratic params
    qp%a1=-3.0*f%K(i)+4.0*f%K(i+1)-f%K(i+2)
    qp%a2=f%K(i+2)-f%K(i)-qp%a1
    qp%i=i
    END SUBROUTINE getqa
    SUBROUTINE fitquad(x0,x1,e,p,y,n)
    IMPLICIT NONE
    REAL(FRK),INTENT(IN)::x0,x1,e
    TYPE(params),INTENT(IN)::p
    REAL(FRK),INTENT(OUT)::y(:)
    INTEGER,INTENT(OUT)::n
    ! Fit quadratics to K in head interval.
    ! Definitions of arguments:
    ! x0,x1 - head interval from x0 to x1.
    ! e     - rel error.
    ! p     - soil type params.
    ! y     - K values defining quads.
    ! n     - no. of values.
    LOGICAL::ok
    INTEGER::i,j,nmax
    REAL(FRK)::a1,a2,dx,e1,e2
    nmax=size(y)
    ! initialise y
    n=5; dx=(x1-x0)/(n-1)
    do i=1,n
      y(i)=Kf(x0+(i-1)*dx,p)
    end do
    ! loop until error ok
    do
      ! check quadratic fits
      ok=.true.
      do i=1,n-4,4
        a1=-3.0_frk*y(i)+4.0_frk*y(i+2)-y(i+4)
        a2=y(i+4)-y(i)-a1
        e1=abs((y(i)+0.25_frk*(a1+0.25_frk*a2))/y(i+1)-one)
        e2=abs((y(i)+0.75_frk*(a1+0.75_frk*a2))/y(i+3)-one)
        ! require errors at test points and monotonicity to be ok
        if (max(e1,e2)>e.or.a1*(a1+two*a2)<zero) then
          ok=.false.
          exit
        end if
      end do
      if (.not.ok) then
        ! double no. of intervals
        if (2*n-1>nmax) then
          write (*,*) "fitquad: too many points"
          stop
        end if
        dx=half*dx
        do i=n,2,-1
          j=2*i-1
          y(j)=y(i)
          y(j-1)=Kf(x0+(j-2)*dx,p)
        end do
        n=2*n-1
      else
        exit
      end if
    end do
    ! remove test points
    j=1
    do i=3,n,2
      j=j+1
      y(j)=y(i)
    end do
    n=j
    CONTAINS
    FUNCTION Kf(h,p)
    IMPLICIT NONE
    REAL(FRK),INTENT(IN)::h
    TYPE(params),INTENT(IN)::p
    REAL(FRK)::Kf
    ! Get SvG K(h) (for h<=hmac1)
    ! Definitions of arguments:
    ! h - head.
    ! p - soil type params.
    REAL(FRK)::R,v1,v2,x
    x=one/(one+exp(p%n*log(h/p%hg)))
    v1=exp(p%m*p%p*log(x))
    v2=one-exp(p%m*log(one-x))
    Kf=p%Kxm*v1*v2**2
    if (h>hmac2) then
      R=Rmac1*(h-hmac2)/(hmac1-hmac2)
      Kf=Kf*exp(R*log(p%Ke/Kf))
    end if
    END FUNCTION Kf
    END SUBROUTINE fitquad
    SUBROUTINE solpar(j,bdj,disj)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::j
    REAL(FRK),INTENT(IN)::bdj,disj
    ! Set soil solute property parameters.
    ! Definitions of arguments:
    ! j    - soil type no.
    ! bdj  - soil bulk density.
    ! disj - dispersivity.
    if (size(fbd)<j) then
      write (*,*) "solpar: not enough space allocated for solutes"
      stop
    end if
    fbd(j)=bdj
    dis(j)=disj
    isotype(j,:)="no" ! will be changed if required in sub setiso
    END SUBROUTINE solpar
    SUBROUTINE setiso(j,isol,isotypeji,isoparji)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::j,isol
    CHARACTER(LEN=2),INTENT(IN)::isotypeji
    REAL(FRK),INTENT(IN)::isoparji(:)
    ! Set soil solute adsorption isotherm and parameters.
    ! Definitions of arguments:
    ! j           - soil type no.
    ! isol        - solute no.
    ! isotypeji   - isotherm code.
    ! isoparji(:) - isotherm params.
    INTEGER::np
    isotype(j,isol)=isotypeji
    np=size(isoparji)
    if (associated(isopar(j,isol)%p)) deallocate(isopar(j,isol)%p)
    if (isotypeji=="Fr") then ! add params to avoid singularity at zero
      allocate(isopar(j,isol)%p(np+2))
      isopar(j,isol)%p=zero
    else
      allocate(isopar(j,isol)%p(np))
    end if
    isopar(j,isol)%p(1:np)=isoparji
    END SUBROUTINE setiso
    FUNCTION Sofh(h,j)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::j
    REAL(FRK),INTENT(IN)::h
    REAL(FRK)::Sofh
    ! Get saturation S from matric head h.
    ! Definitions of arguments:
    ! h   - matric head.
    ! j   - soil type no.
    REAL(FRK)::a1,a2,hr,z
    TYPE(params),POINTER::p
    TYPE(dryparams),POINTER::dp
    if (associated(drypar(j)%p)) then ! check for dry soil wrc
      dp=>drypar(j)%p
      if (h<hjoin) then ! use logarithmic dry soil wrc
        Sofh=dp%Sjoin*log(h/hdry)/drycS
        return ! finished
      end if
    end if
    p=>par(j)
    hr=h/p%hg
    ! get S
    if (hr<p%hr1) then ! quadratic h(S)
      a1=two*p%hr1-(S1-one)*p%hrS1
      a2=p%hr1-a1
      z=two*hr/(a1+sqrt(a1**2+4.0_frk*a2*hr))
      Sofh=one+(S1-one)*z
    else ! vG wrc
      Sofh=exp(-p%m*log(one+exp(p%n*log(hr)))) ! (1+hr**n)**(-m)
    end if
    END FUNCTION Sofh
    FUNCTION weight(j,h,K,phi,dz)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::j
    REAL(FRK),INTENT(IN)::h,K,phi,dz
    REAL(FRK)::weight
    ! Get conductivity weighting for gravity flux calculations.
    ! Definitions of arguments:
    ! j     - soil type no.
    ! h     - matric head.
    ! K     - conductivity.
    ! phi   - MFP.
    ! dz    - flow path length.
    REAL(FRK)::hz,Kz,Khz,phiz,x,a,w
    TYPE(params),POINTER::p
    p=>par(j)
    a=p%n*(p%m*p%p+two) ! vG asymptotic K equiv of BC lambda*eta
    hz=h-gf*dz ! gf is gravity fac in direction of dz
    x=-gf*dz/h
    if (h<p%hpwr.and.(a<=3.0_frk.or.x*(a-3.0_frk)<=4.0_frk)) then
      ! use predetermined approx (ignore dry wrc as w near 0.5 by then)
      w=(60.0_frk+x*(70.0_frk+10.0_frk*a+x*(16.0_frk+a*(5.0_frk+a))))/ &
        (120.0_frk+x*(120.0_frk+x*(22.0_frk+2.0_frk*a**2)))
    else
      call hyofh(hz,j,Kz,Khz,phiz) ! accurate but slower
      w=-((phiz-phi)/(gf*dz)+K)/(Kz-K)
    end if
    weight=min(max(w,zero),one)
    END FUNCTION weight

    END MODULE HPERC2_lib

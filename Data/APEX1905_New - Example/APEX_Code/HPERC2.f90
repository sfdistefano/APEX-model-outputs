    SUBROUTINE HPERC2(phase) !temp cdj fastRE subroutine
    USE PARM
    USE HPERC2_lib, ONLY:dSmax,solve,allo,hypar,Sofh,deallo

    IMPLICIT NONE

    INTEGER,INTENT(IN)::phase
    INTEGER::i,j,tempCount
    REAL::tempTh, tempThr, tempThs, tempS, tempTemp, X1,X2,DZ,Z2,X3,AVW,SUM,ADD,TOT !,tempOr
    REAL::kunsat, tempSST, tempKunsat, tempRatio, tempStore, slopeAngle, Qv !temp cdj richards lateral
    !integrates Peter Ross' fastRE code (FRENUMBERS, FREHYPVGVG, and FREFLOW) with EPIC
    !phase - 1-INITIALIZE, 2-UPDATE, >2-DAILY
    !i, j, ij - indices
    ADD=0.;TOT=0.;SUM=0.
    !------------------------- INITIALIZE -----------------------------
    !initializes the parameters and flux trackers
    IF(phase .EQ. 1)THEN
       IF(ALLOCATED(FRE%ILYR))THEN
            DEALLOCATE(FRE%ILYR)
            CALL DEALLO()
       END IF
       ALLOCATE(FRE%p(MSA),FRE%h0(MSA),FRE%qevap(MSA),FRE%qprec(MSA),FRE%drn(MSA),FRE%evap(MSA),FRE%infil(MSA),FRE%runoff(MSA))
       ALLOCATE(FRE%ILYR(MXLA,MSA))
       FRE%ns=0
       FRE%p(:)=0.5
       FRE%drn(:)=0;FRE%evap(:)=0;FRE%infil(:)=0;FRE%runoff(:)=0
       FRE%qevap(:)=0;FRE%qprec(:)=0
       FRE%ti=-24;FRE%tf=0;FRE%nsteps=0
       dSmax=0.05
       FRE%h0(:)=0.0 !surface head, equal to depth of surface pond
        
       tempCount = 0
       DO j=1,MSA
            DO i=1,NBSL(j)
                ISL=LID(i,j)
                tempCount = tempCount + 1
                FRE%ILYR(i,j)%x=z(ISL,j)*100. !convert from m to cm
                FRE%ILYR(i,j)%Kx=SATC(ISL,j)/(10.) !convert from mm/hr to cm/hr
                FRE%ILYR(i,j)%vgn=VGN(ISL,j)
                FRE%ILYR(i,j)%vgm=1-1/VGN(ISL,j)
                FRE%ILYR(i,j)%hg=-1./VGA(ISL,j)
                FRE%ILYR(i,j)%vgn=VGN(ISL,j)
                FRE%ILYR(i,j)%ths=vgSat(ISL,j)
                FRE%ILYR(i,j)%thr=vgRes(ISL,j)
                FRE%ILYR(i,j)%Kxh=HCL(ISL,j)/(10.) !temp cdj richards lateral convert from mm/hr to cm/hr
                FRE%ILYR(i,j)%jt = tempCount
                FRE%ILYR(i,j)%PRK=0.
                
                FRE%ILYR(i,j)%mn = FRE%ILYR(i,j)%vgm*FRE%ILYR(i,j)%vgn
            END DO
            FRE%ILYR(:,j)%dx=FRE%ILYR(:,j)%x-eoshift(FRE%ILYR(:,j)%x,-1)
            FRE%ILYR(:,j)%Ks=FRE%ILYR(:,j)%Kx+0.1
																
            FRE%ILYR(:,j)%Ksh=FRE%ILYR(:,j)%Kxh+0.1 !temp cdj richards lateral
        END DO
        
        call allo(FRE%ILYR(NBSL(MSA),MSA)%jt,FRE%ns)
        
        DO j=1,MSA
            do i=1,NBSL(j)
                call hypar(FRE%ILYR(i,j)%jt,FRE%ILYR(i,j)%thr,FRE%ILYR(i,j)%ths,FRE%ILYR(i,j)%hg,FRE%ILYR(i,j)%Kx,FRE%ILYR(i,j)%mn,FRE%p(j),FRE%ILYR(i,j)%Ks)
            end do
        END DO
        
    ELSE
        !----------------------------- UPDATE DAILY ---------------------
        !updates the fastRE soil parameters and drivers from the previous day and executes the fastRE for the current day
        DO i=1,NBSL(ISA)
            ISL=LID(i,ISA)
            CPFH(ISL,ISA)=0. 
            
            FRE%ILYR(i,ISA)%S=(SWST(ISL,ISA)/(FRE%ILYR(i,ISA)%dx*10.)-FRE%ILYR(i,ISA)%thr) /(FRE%ILYR(i,ISA)%ths-FRE%ILYR(i,ISA)%thr)
            FRE%ILYR(i,ISA)%S = max(FRE%ILYR(i,ISA)%S, 1.0E-9)
            tempTh=SWST(ISL,ISA)/(FRE%ILYR(i,ISA)%dx*10.)
            tempThs=FRE%ILYR(i,ISA)%ths
            tempThr=FRE%ILYR(i,ISA)%thr
            tempS=(tempTh-tempThr)/(tempThs-tempThr)
            tempTemp = 1.
        END DO
        
        CPVH(IDO)=0. !Pipe flows are NOT considered with the Richards's percolation method. Jaehak 2017
        FRE%ti=FRE%ti+24.; FRE%tf=FRE%tf+24.
        FRE%qprec(ISA)=SEP/(24.*10.) !Infiltration to the first layer. converted from mm/day to cm/hr
        
																											  
																											
        call solve(FRE%ti,FRE%tf,FRE%qprec(ISA),FRE%qevap(ISA),NBSL(ISA),FRE%ns,FRE%ILYR(:,ISA)%dx,FRE%ILYR(:,ISA)%jt,FRE%h0(ISA),FRE%ILYR(:,ISA)%S,FRE%evap(ISA),FRE%runoff(ISA),FRE%infil(ISA),FRE%drn(ISA),FRE%ILYR(:,ISA)%PRK,FRE%nsteps) !temp cdj richards lateral added SSF and IFLO        
        
        FRE%ILYR(:,ISA)%fvwc=FRE%ILYR(:,ISA)%thr + (FRE%ILYR(:,ISA)%ths - FRE%ILYR(:,ISA)%thr) * FRE%ILYR(:,ISA)%s
        SST(ISA)=0.
 
        SW(ISA)=0.
        DO i=1,NBSL(ISA)
            ISL=LID(i,ISA)
            SWST(ISL,ISA)=FRE%ILYR(i,ISA)%fvwc*FRE%ILYR(i,ISA)%dx*10. !cm to mm
            SEP= FRE%ILYR(i,ISA)%PRK*10. !mm
            SW(ISA)=SW(ISA)+SWST(ISL,ISA)
        END DO
        
        DO i=1,NBSL(ISA)
            ISL=LID(i,ISA)
            SWST(ISL,ISA)=FRE%ILYR(i,ISA)%fvwc*FRE%ILYR(i,ISA)%dx*10. !cm to mm
            SEP= FRE%ILYR(i,ISA)%PRK*10. !mm
            tempKunsat = kunsat(HCL(ISL,ISA), FRE%ILYR(i,ISA)%VGM, FRE%ILYR(i,ISA)%S) !effective lateral conductivity mm/hr
            IFLO = 2  !JAEHAK temporary setting
        
            !begin temp cdj richards lateral
            IF(IFLO.EQ.1)THEN
                !Warrick et al., 2008 WR eq. 14 and assuming homogenous isotropic layers
                slopeAngle = ATAN(STP(ISA))*360./(2.*3.14159)
                tempRatio = 1-SIN(slopeAngle*2.*3.14159/360.)
                Qv = FRE%ILYR(i,ISA)%PRK / tempRatio
                tempSst = Qv - FRE%ILYR(i,ISA)%PRK
            ELSE IF(IFLO.EQ.2)THEN
                !approximation of standard EPIC approach
                tempRatio = HCL(ISL,ISA)/SATC(ISL,ISA)
                tempSst = FRE%ILYR(i,ISA)%PRK * tempRatio * 10.
                tempSst = min(tempSst, (SWST(ISL,ISA)-FC(ISL,ISA)) ) 
                tempSst = max(tempSst, 0.)
                DZ=FRE%ILYR(i,ISA)%dx*0.01 !cm to m
                AVW=tempSst !mm
                IF (AVW>0.001) THEN
                    X1=24./(PO(ISL,ISA)-FC(ISL,ISA))
                    X2=X1*tempKunsat
                    IF(ISL/=IDR(ISA))THEN      
                        Z2=PRMT(90)*DZ*X2/SPLG(ISA)
                        X3=AVW*(1.-EXP(-Z2))
                        X1=MIN(X3,.001*X3*SPLG(ISA)/RCHL(ISA))
                        SST(IDO)=X1
                        QRF(IDO)=X3-X1
                    ELSE
                        QRF(IDO)=X3
                        SST(IDO)=0.
                    END IF
                ELSE
                    QRF(IDO)=0.
                    SST(IDO)=0.
                END IF
                
                IF(RSAE(ISA)>0.)QRF(IDO)=0.
                
            ELSE IF(IFLO.EQ.3)THEN
                !Darcian approach
                tempKunsat = 24. * kunsat(HCL(ISL,ISA), FRE%ILYR(i,ISA)%VGM, FRE%ILYR(i,ISA)%S)
                tempSst = tempKunsat * STP(ISA) !UPS -> STP to transfer from EPIC to APEX Jaehak 2017  
            END IF
            SWST(ISL,ISA)=MAX(1.E-5,SWST(ISL,ISA)-SST(IDO)-QRF(IDO))
    
            IF(ISL/=IDR(ISA))THEN
                SUM=SUM+QRF(IDO)
            ELSE
                SMM(17,MO,ISA)=SMM(17,MO,ISA)+QRF(IDO)
                VAR(17,ISA)=QRF(IDO)
                QDR(IDO)=QRF(IDO)
            END IF
            ADD=ADD+SST(IDO)
            TOT=TOT+CPVH(IDO)
            SSF(ISL,ISA)=SST(IDO)
            QSF(ISL,ISA)=QRF(IDO)
            CPFH(ISL,ISA)=CPVH(IDO)
            PKRZ(ISL)=SEP
           
        END DO
        SST(IDO)=ADD
        QRF(IDO)=SUM
        
    END IF

    END SUBROUTINE HPERC2

    !temp cdj richards lateral
    Function kunsat(ksat, mvg, se)
    implicit none
    real kunsat, ksat
    REAL(SELECTED_REAL_KIND(15)) :: se, mvg
    real, parameter :: l=0.5

    if(se > 0.999) then
        kunsat=ksat
    else
        kunsat = ksat*(se**l)*(1.-(1.-Se**(1./mvg))**mvg)**2. !mm/h
        kunsat = Max(0., kunsat) 
        kunsat = Min(ksat, kunsat) 
        if(kunsat < 1.E-10) THEN
            kunsat = 0.0
        endif
    end if

    return
    end function kunsat





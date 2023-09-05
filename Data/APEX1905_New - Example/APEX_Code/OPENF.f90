      SUBROUTINE OPENF
!     APEX1905
!     THIS SUBPROGRAM OPENS FILES.
      USE PARM
      CHARACTER(4)::AXT
      DIMENSION AXT(50)
	  DATA AXT/".OUT",".MAN",".SUS",".ASA",".SWT",".DPS",".MSA",".AWP",&   !1,2,3,4,5,6,7,8
      ".DHY",".WSS",".SAD",".HYC",".DRS","    ",".MWS",".DWS",".AWS",&     !9,10,11,12,13,14,15,16,17
      ".DGZ",".DUX",".DDD",".ACN",".DCN",".SCX",".ACY",".EFR",".EHY",&     !18,19,20,21,22,23,24,25,26
      ".APS",".MSW",".DPW",".SPS",".ACO",".SWN",".AGZ",".SAO",".RCH",&     !27,28,29,30,31,32,33,34,35
      ".ERX",".GZM",".STR",".MRH",".MGR",".DNC",".DHS",".MSX",".DGN",&     !36,37,38,39,40,41,42,43,44
      ".DPD",".ASL",".MS5",".AS5",".CAG",".CMP"/                           !45,46,47,48,49,50
      OPEN(KW(1),FILE=ASTN//AXT(1))
	  DO I=2,48
	      IF(AXT(I)/="    ".AND.KFL(I)>0)OPEN(KW(I),FILE=ASTN//AXT(I))
	  END DO
!      IF(KFL(49)>0)THEN
!          OPEN(KW(49),FILE=ASTN//AXT(49),ACCESS='APPEND')
!      END IF
      RETURN
      END
! Read files KR(), complied by Liwang Ma
! 	  DATA AXT/".SIT","FPARM","FTILL","FCROP",".SUB","FMLRN","","FPEST",&               !1,2,3,4,5,6,7,8
!      "FFERT","FTR55","APEXRUN","APEXFILE","FSOIL","    ","FOPSC",".OPC","FPRNT",&     !9,10,11,12,13,14,15,16,17
!      "FWPM","FWIND","APEXCONT","FHERD"," ","FSITE","FSUBA","FWLST","APEXDIM",&        !18,19,20,21,22,23,24,25,26
!      "FPSOD","SOIL36K","FRFDT"," ","GRZP",".SWN",".AGZ",".SAO",".RCH",&               !27,28,29,30,31,32,33,34,35
!      ".ERX",".GZM",".STR",".MRH",".MGR",".DNC",".DHS",".MSX",".DGN",&                 !36,37,38,39,40,41,42,43,44
!      ".DPD",".ASL",".MS5",".AS5"/                                                     !45,46,47,48

      
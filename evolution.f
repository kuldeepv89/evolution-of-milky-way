C     1) INFALL TIME_SCALE t1; 2) INFALL TIME_SCALE t2; 3) DELAY TMAX; 4) TIME (AGE=13.7-TIME); 5) [FE/H]; 6) [alpha/Fe]; 7) [M/H]    
              
      SUBROUTINE MW(at1,at2,atauin,aTTspi,aAFEHv, aaspiv,aAMFv)
      PARAMETER (NMAX=37) 






      
C
C C	X(I,1)  H              
c	X(I,2)  D              
C       X(I,3)  He3            
C       X(I,4)  He4            
C       X(I,5)  C              
C       X(I,6)  O              
C       X(I,7)  N14            
C       X(I,8)  C13            
C       X(I,9)  neutron rich   
C       X(I,10) Ne             
C       X(I,11) Mg             
C       X(I,12) Si             
C       X(I,13) S              
C       X(I,14) Ca             
C       X(I,15) Fe             
C       X(I,16) Cu             
C       X(I,17) Zn             
C       X(I,18) Ni             
C       X(I,19) Kr             
C       X(I,20) Li             
C       X(I,21) N15            
C       X(I,22) K              
C       X(I,23) Sc             
C       X(I,24) Ti             
C       X(I,25) V              
C       X(I,26) Cr             
C       X(I,27) Mn             
C       X(I,28) Co             
C       X(I,29) Ba slow
C       X(I,30) Ba rapid
C       X(I,31) Sr rapid
C       X(I,32) Eu 
C       X(I,33) La
C       X(I,34) Y
C       X(I,35) Zr
C       X(I,36) Al
C       X(I,37) Na
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Best model (CMR2001; RMVD2001) corretto
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C This model computes the chemical evolution of the galaxy of the 
C elements ordered as : H,D, He3,He4,C,O,N14,C13,neutron rich, Ne, 
C Mg,Si,S,CA,Fe,Cu,Zn,Ni,Kr,Li,N15,K,Sc,Ti,V,Cr,Mn,Co
C +Eu+Ba
C
C Prescriptions for nuclesynthesis :
C
C     low mass stars
C
C van den Hoek & Groenewegen (1997) for 0.9 - 8 Msun,
C______________________________________________________
C
C       massive stars (SNeII)
C
C    Thielemann, Nomoto & Hashimoto (1996) 
C              or
C     Woosley and Weaver (1995)
C_______________________________________________________
C
C           Type Ia supernovae
C
C    Thielemann et al.(1993) per SNeIa,
C                or
C             Iwamoto 
C________________________________________________________
C
C                   Novae
C
C          Jose` ed Hernanz (1998)
C
C__________________________________________________________
C


CCC   variabili common per il bario!!!
      include 'lettura.inc'
      include 'elementi.inc'
CCC         
      real at1,at2,atauin
      real t1,t2,tauin
      REAL TTspi,AFEHv,aspiv,AMFv
      REAL MI,MS,MU
      integer radius,gr,threshold,insideout,halo
      character*80 titre
      character*20 filename,filename2
      character*5  modello
      character*35  filename3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Variabili condivise con function Q (Francesco)
      REAL M1,M2
      DIMENSION XSSTAR(NMAX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei bursts di nova.
      DIMENSION RNOV(NMAX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei raggi cosmici.
      DIMENSION COSM(NMAX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DIMENSION RR3(NMAX),RR4(NMAX),PRZW1(NMAX),PRZW2(NMAX),PRZW3(NMAX)
      DIMENSION VN(NMAX),TW(NMAX),XCAL(NMAX),XTIN(NMAX),XSTAR(NMAX)
      DIMENSION XSTIN(NMAX)
      DIMENSION BO(3000),RR1(NMAX),RR2(NMAX)
      DIMENSION FALL(3000),SRT(3000),SMS(3000),APUNT(3000),RATE(3000)
      DIMENSION R1(NMAX),SM(200),AMU(200),PARZW(NMAX),XINF(NMAX),T(3000)
      DIMENSION GO(3000),G(3000,NMAX),PAMU(200),XP(NMAX),IND(NMAX)
      DIMENSION WI(NMAX)
      DIMENSION QM(NMAX,NMAX),GP(NMAX),FA(NMAX),X(3000,NMAX),WI1(NMAX)
      DIMENSION WI2(NMAX)
      DIMENSION GRAT(3000)
      DIMENSION AMLI(50),ALI0(50),ALI1(50),ALI2(50),ALI3(50),ALI4(50)
      DIMENSION PMASS(50),YLI0(50),YLI1(50),YLI2(50),YLI3(50),YLI4(50),
     $     YLI5(50),YLI(50)
      DIMENSION COSN1(50),COSN2(50)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei bursts di nova.
      DIMENSION RMAS(3000),RNUM(3000),RNOVAE(3000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       DIMENSION aTTspi(3000),aAFEHv(3000),aaspiv(3000),aAMFv(3000)
      DIMENSION TTspi(3000),AFEHv(3000),aspiv(3000),AMFv(3000)
      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/FR/XINF,FALL,FA,SG
      COMMON/TRE/AMU,GP,XP,WI,PARZW,WI1,WI2,PRZW1,PRZW2,PRZW3
      COMMON/A/QM
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     $     GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     $     KMAX,I,ITER,IND,KDEMAX
      COMMON/CUATRO/COSN1,COSN2
      COMMON/COM/MAXO
      COMMON/DIE/A,B,PAMU,K2,K1,KB,KMAX2,K3,KB3
      COMMON/DIS/EXPO,T1,T2,ETA,RO,RS,AR,BR,RI,RD,SD,SN,AO,
     $R,SMO,SMC,R1,VMR,RN,VNSOL,SRT,ARO,BRO,SMS,GAS,APUNT,
     $RATE,RIEST,RSCA,AC,TW,VN,THALO
      COMMON/MAO/IVO,IPER
      COMMON/STAR/BO,TT1,ISTAR
      COMMON/YIE/ KYIEL
      COMMON/DATA/LET
c     !COMMON/LUCKY/IKU
      COMMON/M/NMAT,NNMAT
      COMMON/WOO/AMLI,ALI0,ALI1,ALI2,ALI3,ALI4
      COMMON/ZOO/PMASS,YLI0,YLI1,YLI2,YLI3,YLI4,YLI5,YLI
      COMMON/Z/ZETA
      COMMON/ELDS/IELD
      COMMON/SN/RR1,RR2,RR3,RR4,COST
      COMMON/NEW/ASNII,GAMMA
      COMMON/PIP/GRAT,DTGRAT
      COMMON/RAF/TAUIN,TNEW,gr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Blocco COMMON per la nuova formulazione di ASP
      COMMON/HEL/ASP
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei bursts di nova.
      COMMON/NOVA/RMAS,RNUM,RNOVAE,WDM,WD,RANOV
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei raggi cosmici.
      COMMON/RAYS/COSM
      COMMON/mod/threshold,modello
     
       COMMON/MCMC/TTspi,AFEHv, aspiv,AMFv
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cf2py intent(in) at1, at2, atauin
Cf2py intent(out) aTTspi, aAFEHv, aaspiv, aAMFv     
      

      gr=0

    
     

    


           
        FILENAME='pasadena4.dat'
        FILENAME2='WW95B.dat'
        
C        write(*,*) 'which is the name of the model?'
c        read(*,*) modello
        modello='prova'

c        write(*,*) 'Inside out yes or no?(1/0)'
c        read(*,*) insideout
        insideout=0

c        write(*,*) 'halo constant (0) or 1/r (1)?'
c        read(*,*) halo
        halo=0

c        write(*,*) 'halo threshold yes or no?(1/0)'
c        read(*,*) threshold
        threshold=1


        



        
        OPEN(5,FILE=FILENAME,STATUS='old')

C
C
C       File 17 contains the yields for SNII
C       It complements the input for file 5  
C

       OPEN(17,FILE=FILENAME2,STATUS='OLD')


C      subroutine di lettura dati del bario messi poi nel blocco common!
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       
       call leggi                                             
       call leggi3
c       call leggi4
       call leggila
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     +-------------------------------------------------+
C     |                                                 |
C     |                 Output files.                   |
C     |                                                 |
C     +-------------------------------------------------+


C     +-------------------------------------------------+
C     |                                                 |
C     |  output files from Cescutti ( a bit             |
C     |                          different from old)    |
C     |                                                 |
C     +-------------------------------------------------+

c        open(88,FILE='infall_2I.gr'           ,STATUS='unknown')

c      open(11,FILE='cesc.model5'           ,STATUS='unknown')
c      OPEN(60,FILE='gratton8h.model5'      ,STATUS='unknown')
c      open(62,FILE='SOL08_best.gr' ,STATUS='unknown')
c      OPEN(63,FILE='rate8h.model5'         ,STATUS='unknown')
c      OPEN(66,FILE='distr8h_best.model5'        ,STATUS='unknown')
c      OPEN(67,FILE='GR08_best.gr' ,STATUS='unknown')
c      OPEN(69,FILE='lpuves.model5'            ,STATUS='unknown')
c      OPEN(80,FILE='fis8h08_2I_best.model5',STATUS='unknown')


c      write(11,*) 't Fe/H Ferro Barios  Barior Stronzior Europio'
c      write(62,*) '#SOL08.gr', radius
c      write(62,7779)  '#T','Z','SFR', '[Fe/H]','[He/Fe]','[Li/Fe]',
c     & '[C/Fe]','[C13/Fe]',
c     & '[O/Fe]', '[N/Fe]', '[Ne/Fe]','[Mg/Fe]','[Si/Fe]','[S/Fe]', 
c     & '[Ca/Fe]','[Cu/Fe]','[Zn/Fe]','[Ni/Fe]','[Kr/Fe]','[K/Fe]',
c     & '[Sc/Fe]','[Ti/Fe]','[V/Fe]','[Cr/Fe]','[Mn/Fe]','[Co/Fe]',
c     & '[Ba/Fe]','[Eu/Fe]','[La/Fe]','[Sr/Fe]','[Y/Fe]','[Zr/Fe]'
c      write(67,*) 'GR08.gr', radius
c      write(67,7779)'T','Z','SFR','log(Fe/H)','log(He/H)','log(Li/H)',
c     & 'log(C/H)','log(C13/H)','log(O/H)', 'log(N/H)' ,'log(Ne/H)',
c     & 'log(Mg/H)','log(Si/H)','log(S/H)' ,'log(Ca/H)','log(Cu/H)', 
c     & 'log(Zn/H)','log(Ni/H)','log(Kr/H)','log(K/H)','log(Sc/H)',
c     & 'log(Ti/H)','log(V/H)','log(Cr/H)','log(Mn/H)','log(Co/H)',
c     & 'log(Ba/H)','log(Eu/H)','log(La/H)','[Sr/Fe]','[Y/Fe]','[Zr/Fe]'


c      write(69,*)'time [Fe/H] [O/Fe] [Mg/Fe] [Si/Fe] [Ca/Fe] [Zn/Fe] Ni'
c      write(80,*)'#time  gaz density infall rate  densite tot SNII SNI' 


 7779 FORMAT(32A13)  



        


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      OPEN(110,FILE='fesolar.dat')
C      OPEN(111,FILE='nesisolar.dat')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      LET=0
      INDEX=0

      READ(5,10) abT2,ETA,A0,AC,SMO,SMC,RO,RS,GAS,TOL,
     $EPSI,CK,UM,M,NNMAX,KMAX,IMAX,KDEMAX
 10   FORMAT(2E12.5/4E12.5/2E12.5/5E12.5/I3/4I10)

      READ(5,15) (R1(K),K=1,12)
 15   FORMAT(6E12.5)

      READ(5,30) (XINF(K),K=1,NMAX)
 30   FORMAT(5E12.5)

      READ(5,20) (AMU(K),K=1,KMAX)
 20   FORMAT(5E12.5)

      READ(5,16) RSCA,VDTEM
 16   FORMAT(2E12.5)

      READ(5,83) (TW(KK),KK=1,12)
      READ(5,83) (VN(K),K=1,12)
 83   FORMAT(4E12.5)

      VVDTEM=VDTEM
      KYIEL=KMAX


c  Some other parameters!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
      SB=1.287E+05
      RB=11.347
      RS=2.
      RF=R1(1)
         eta=13.7
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RD=3.5
      SD=64.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      BULGE=6.28E6*SB*(ALOG(1.+RS/RB)-RS/(RB+RS))

      DISCO=6.28E6*SD*RD**2*((1.+RS/RD)*EXP(-RS/RD)-(1.+RF/RD)
     $*EXP(-RF/RD))

      KHH=KMAX-1
      VNSOL=1.3

C     Here the program starts a loop over different 
C           galactocentric zones 
C     (8 is the index for the radius of solar neigbourhood)

      DO K=8,8

   
C
C     Clean and safe start . Initialisation
C
            DO JK=1,3000
            
            TTspi(JK)=0.
            AFEHv(JK)=0.
            Aspiv(JK)=0.
            AMFv(JK)=0.
         enddo   


               DO JK=1,3000
            
            aTTspi(JK)=0.
            aAFEHv(JK)=0.
            aAspiv(JK)=0.
            aAMFv(JK)=0.
         enddo   


         
         DO JK=1,NMAX
            
            GP(JK)=0.
            XP(JK)=0.
            G(I,JK)=0.
            X(I,JK)=0.
            XSTIN(JK)=0.
            XCAL(JK)=0
            XTIN(JK)=0.
            XSTAR(JK)=0.
            WI1(JK)=0.
            WI2(JK)=0.
            PRZW1(JK)=0.
            PRZW2(JK)=0.
            PRZW3(JK)=0.
            RR3(JK)=0.
            RR4(JK)=0.
            FA(JK)=0.
            WI(JK)=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei bursts di nova.
            RNOV(JK)=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei raggi cosmici.
            COSM(JK)=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            RR1(JK)=0.
            RR2(JK)=0.
            
         ENDDO
         
         
         GAMMA=0.
         ASNII=0.
         
C     
C     End initialisation 
C     
C     T1=TW(K)

         R=R1(K)

c         WRITE(60,131) R
 131     FORMAT(10X,'EVOLUZIONE CHIMICA AL RAGGIO',5X,'R=',F10.5//)

c         WRITE(60,9999) CK
 9999    FORMAT(1X,'CK=',1E12.5)

         VMR=SD
C     Valido solo nei dintorni solari.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Infall ritardato.
C     delayed infall
C     
         t1=at1
         t2=at2
         tauin=atauin
c         T1=.1
        
c            t2=8.
       
c         TAUIN=4.3
         THALO=1.0
c      write(*,*) t1,t2,tauin, 'BOH'
C     Nel solar ring, entro 1 Kpc di distanza dal piano galattico:
C        SIGMAH=17.0    dens. sup. di massa dell'alone (+ thick disk)
C        SIGMAD=54.0    dens. sup. di massa del disco (thin disk)

         SH=136.
        
            sigmah=8
         

         AR=(SIGMAH-GAS)/T1*(1.-EXP(-THALO/T1))

         ARO=AR
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     WRITE(60,110) AR
C     110     FORMAT(10X,'AR=',1E12.5)

c         WRITE(60,115) VMR
 115     FORMAT(10X,'DENSITA SIGMA R,T=',1E12.5//)
      
         
         CALL CHEM
         
         VDTEM=VVDTEM
      
      
         
      ENDDO
      
      
      
      CLOSE(60)
      CLOSE(70)
      CLOSE(80)
      

      
            DO JK=1,3000
            
            aTTspi(JK)=TTspi(JK)
            aAFEHv(JK)=AFEHv(JK)
            aAspiv(JK)=Aspiv(JK)
             aAMFv(JK)= AMFv(JK)
         enddo   


      
      return
      END






      SUBROUTINE CHEM
       
      PARAMETER (NMAX=37)
      include 'lettura.inc'
      character*5  modello
      REAL TTspi,AFEHv,aspiv,AMFv
      REAL MI,MS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei bursts di nova.
      DIMENSION RMAS(3000),RNUM(3000),RNOVAE(3000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei raggi cosmici.
      DIMENSION COSM(NMAX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DIMENSION PRZW1(NMAX),PRZW2(NMAX),PRZW3(NMAX)
      DIMENSION SM(200),GRAT(3000)
      DIMENSION TW(NMAX)
      DIMENSION VN(NMAX)
      DIMENSION BO(3000),STA(3000),XSTIN(NMAX)
      DIMENSION R1(NMAX),AMU(200),PARZW(NMAX),XINF(NMAX),T(3000)
      DIMENSION GO(3000),G(3000,NMAX),PAMU(200),XP(NMAX),IND(NMAX)
      DIMENSION WI(NMAX)
      DIMENSION FALL(3000),SRT(3000),SMS(3000),APUNT(3000),RATE(3000)
      DIMENSION QM(NMAX,NMAX),FA(NMAX),X(3000,NMAX)
      DIMENSION GP(NMAX),WI1(NMAX),WI2(NMAX),XCAL(NMAX),XTIN(NMAX)
      DIMENSION XSTAR(NMAX)
      DIMENSION COSN1(50),COSN2(50)
      DIMENSION Q1(NMAX),QJK(NMAX,NMAX),P(NMAX),ENZO(200),SSM(200)
      integer threshold



       DIMENSION TTspi(3000),AFEHv(3000),aspiv(3000),AMFv(3000)
      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/FR/XINF,FALL,FA,SG
      COMMON/TRE/AMU,GP,XP,WI,PARZW,WI1,WI2,PRZW1,PRZW2,PRZW3
      COMMON/A/QM
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     $GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     $KMAX,I,ITER,IND,KDEMAX
      COMMON/COM/MAXO
      COMMON/DIE/A,B,PAMU,K2,K1,KB,KMAX2,K3,KB3
      COMMON/DIS/EXPO,T1,T2,ETA,RO,RS,AR,BR,RI,RD,SD,SN,AO,
     $R,SMO,SMC,R1,VMR,RN,VNSOL,SRT,ARO,BRO,SMS,GAS,APUNT,
     $RATE,RIEST,RSCA,AC,TW,VN,THALO
      COMMON/MAO/IVO,IPER
      COMMON/STAR/BO,TT1,ISTAR
      COMMON/CUATRO/COSN1,COSN2
c     !COMMON/LUCKY/IKU
      COMMON/DATA/LET
      COMMON/Z/ZETA
      COMMON/M/NMAT,NNMAT
      COMMON/CLAUDIA/AFH
      COMMON/ELDS/IELD
      COMMON/NEW/ASNII,GAMMA
      COMMON/PIP/GRAT,DTGRAT
      COMMON/RAF/TAUIN,TNEW,gr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Blocco COMMON per la nuova formulazione di ASP
      COMMON/HEL/ASP
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei bursts di nova.
      COMMON/NOVA/RMAS,RNUM,RNOVAE,WDM,WD,RANOV
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei raggi cosmici.
      COMMON/RAYS/COSM
      COMMON/mod/threshold,modello
       COMMON/MCMC/TTspi,AFEHv, aspiv,AMFv
      
c     !IKU=0
      THR=0.01
      THALO=1.0
c      TAUIN=4.3
      CK=1.5
      IGRAT=0
      IAB=0
      EXPO=2.*(CK-1.)

      TP=0.
      MS=AMU(KMAX)
      MAXO=NMAX
      TAUINF=TAU(MS)

      I=1
      T(1)=0.
      GRAT(1)=2.0
      BO(1)=0.

C      APUNT(1)=AR+BR


      APUNT(1)=AR
      SRT(1)=GAS
      SMS(1)=GAS
      SG=VMR


C Densita` superficiale di massa totale alla Gilmore.

      GO(1)=SRT(1)/SG
     
c      IF ((GO(1)*SG.LE.4.0).and.(threshold.eq.1)) THEN
c         GRAT(1)=0.
c      ELSE
c         GRAT(1)=2.0
c      END IF
     
      
      RATE(1)=GRAT(1)*(SRT(1)/SMS(1))**EXPO*(SG/SRT(1))**(CK-1.)


      DO L=1,NMAX
         X(1,L)=XINF(L)
      ENDDO


      FALL(1)=APUNT(1)/SG
      DO  L=1,NMAX
         G(1,L)=X(1,L)*GO(1)
      ENDDO

      ZETA=0.

      DO LZ=5,NMAX
         ZETA=ZETA+X(1,LZ)
      ENDDO


 300  CONTINUE
    
      CALL TEMPO

      KRIP=0

 301  CONTINUE
      

      ITER=0
      IELD=0
      TN=T(I)+DTEM
      TT1=TN
      IF(TN.GT.ETA) GO TO 600

      IF(I-1)40,40,50

 40   CONTINUE

      TZ=DTEM/2.
C      SRP=AR*T1*(1.-EXP(-TZ/T1))+BR*T2*(1.-EXP(-TZ/T2))+GAS
C      FAP=(AR*EXP(-TZ/T1)+BR*EXP(-TZ/T2))/SG
      SRP=AR*T1*(1.-EXP(-TZ/T1))+GAS
      FAP=(AR*EXP(-TZ/T1))/SG
C      SMSP=ARO*T1*(1.-EXP(-TZ/T1))+BRO*T2*(1.-EXP(-TZ/T2))+GAS
      SMSP=ARO*T1*(1.-EXP(-TZ/T1))+GAS
      RAP=GRAT(1)*(SRP/SMSP)**EXPO*(SG/SRP)**(CK-1.)
      AN=RATE(1)*GO(1)**(CK-1.)
      
      IF (AN.EQ.0.) GO TO 401
      GASP=0.
      
      DO L=1,NMAX

         FA(L)=(FALL(1)/AN)*(1.-EXP(-AN*DTEM))*XINF(L)
         GP(L)=G(1,L)*EXP(-AN*DTEM)+FA(L)
         GASP=GASP+GP(L)

      ENDDO

      GO TO 70


 401  CONTINUE
      GASP=0.
      DO  L=1,NMAX

         GP(L)=SRP*XINF(L)/SG
         GASP=GASP+GP(L)
         WI(L)=FAP*XINF(L)
         IND(L)=1

      ENDDO

      GO TO 402


 50   CONTINUE                                                         

C LE SCHEDE CHE SEGUONO VANNO BENE NEL CASO KNON KAPPA NONO UNO

c      TAUIN=4.3
      
      IF (T(I).LT.TAUIN) THEN
         BR=0.
         BRO=0.
      ELSE IF (T(I).GE.TAUIN) THEN
         BR=(VMR-GAS)/T2*(1.-EXP(-(ETA-TAUIN)/T2))
         BRO=BR
      END IF

      IF (T(I).LT.TAUIN) THEN
         TTAU=0.
      ELSE
         TTAU=T(I)-TAUIN
      END IF
      
      APUNT(I)=AR/EXP(T(I)/T1)+BR/EXP(TTAU/T2)

c      WRITE(6,*) APUNT(I),AR,BR,TAUIN

C      APUNT(I)=AR/EXP(T(I)/T1)+BR/EXP(T(I)/T2)
C      SRT(I)=AR*T1*(1.-1./EXP(T(I)/T1))+BR*T2*(1.-1./EXP(T(I)/T2))
C     $+GAS
C      SMS(I)=ARO*T1*(1.-1./EXP(T(I)/T1))+BRO*T2*(1.-1./EXP(T(I)/T2))
C     $+GAS
      SRT(I)=AR*T1*(1.-1./EXP(T(I)/T1))+BR*T2*(1.-
     $1./EXP(TTAU/T2))+GAS
      SMS(I)=ARO*T1*(1.-1./EXP(T(I)/T1))+BRO*T2*(1.-
     $1./EXP(TTAU/T2))+GAS
      GRAT(I)=VNSOL

C Pepe's SFR.
C      RATE(I)=VNSOL*SRT(I)**EXPO
C      RATE(I)=VNSOL
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c aggiunta iatus.
c      if(igrat.eq.0) then
c      if(aoh.ge.-0.3) then
c -------------------------------------
c se t < 2 Gyr (halo) entao eficiencia das SFRhalo = valor dado para GRAT(I)
c se t >= 2 Gyr (disco) eficiencia da SFRdisco = VNSOL (valor dado no inicio)
c note que eficiencia disco < eficiencia halo - evolucao halo e' mais rapida

      if (t(i).le.tauin-0.3) then
          vnsol=1.3
          grat(i)=vnsol
       else
          vnsol=1.3
          grat(i)=vnsol
       end if



      
c      IF(T(I).LE.THALO) THEN
c         if((GO(I)*sg.le.4.0).and.(threshold.eq.1)) then
c            grat(i)=0.
c         else
c            grat(i)=2. 
c         end if
c      ELSE
c         if((GO(I)*sg.le.7.0).and.(threshold.eq.1)) then
c            grat(i)=0.
c         else
c            grat(i)=vnsol
c         end if
c      END IF
c --------------------------------------

c      write(6,*) grat(i),ck,rate(i)

      RATE(I)=GRAT(I)*(SRT(I)/SMS(I))**EXPO*(SG/SRT(I))**(CK-1.)
      RAP=GRAT(I)*(SRP/SMSP)**EXPO*(SG/SRP)**(CK-1.)



C      WRITE(6,*) T(I),GRAT(I),CK
c fine aggiunta.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FALL(I)=APUNT(I)/SG
      TDIF=T(I)-T(I-1)
      GASP=0.
C      IF (T(I).GE.2.0.AND.T(I).LE.2.05) GO TO 9364
      DO     L=1,NMAX
         GP(L)=(G(I,L)-G(I-1,L))/TDIF*DTEM+G(I,L)
         GASP=GASP+GP(L)
      ENDDO

      TM=T(I)+DTEM/2.
      SRP=AR*T1*(1.-EXP(-TM/T1))+BR*T2*(1.-EXP(-(TM-TAUIN)/T2))+GAS
      SMSP=ARO*T1*(1.-EXP(-TM/T1))+BRO*T2*(1.-EXP(-(TM-TAUIN)/T2))+GAS
      AP=AR*EXP(-TM/T1)+BR*EXP(-(TM-TAUIN)/T2)
      FAP=AP/SG
      DO    L=1,NMAX
         XP(L)=GP(L)/GASP
      ENDDO

c     write(6,*) (xp(l),l=1,nmax)
 70   CONTINUE
      DO    K=1,NMAX
         IND(K)=0
      ENDDO
C     Inizio iterazione.

C      MAMMA=0

      KDELTA=0
      GPRE=GO(I)
      IF(IAB.EQ.1.) GO TO 786
      If(RATE(I).EQ.0.) go to 785
      IAB=1
      TNEW=T(I)
      
 786  CONTINUE

      CALL CALWI
      
 785  CONTINUE

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Parte per tener conto dei bursts di nova.

      RMAS(I)=WDM
      RNUM(I)=WD
      RNOVAE(I)=RANOV

c      WRITE(63,9861) T(I),RMAS(I),RNUM(I),RNOVAE(I)

 9861 FORMAT(1X,'T=',1E12.5,1X,'RMAS=',1E12.5,1X,'RNUM=',1E12.5,1X,
     $'RNOVAE=',1E12.5)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      wpad=0.
      
      do k1=1,nmax
         wpad=wpad+wi1(k1)
      enddo
      
c      write(6,3002)t(i),wpad
 3002 format(1x,1e12.5,1x,'wpad=',1e12.5)
      
      ITER=1

 200  CONTINUE
      GG1=GPRE**(CK-1.)
      GG2=GASP**(CK-1.)
      GG=(GG1+GG2)*0.5
      RA=RAP*GG
      DERLOG=CK-1.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(ra.eq.0.) go to 658
      PO=RA*DTEM
      P1=1./EXP(PO)
      IVO=0
      IPER=0
      DO 100 J=1,NMAX
      P2=WI(J)/RA
      AI=P1*(G(I,J)-P2)+P2
      BI=0.5*DERLOG*((P1-1.)*P2-PO/EXP(PO)*(G(I,J)-P2))
      DELTA=(GP(J)-AI)/(BI-GP(J))
      GP(J)=GP(J)*(1.+DELTA)
      IF(ABS(DELTA)-EPSI) 120,120,110
 120  IND(J)=1
 110  CONTINUE
 100  CONTINUE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      go to 111
 658  continue
      gasp=0.
      do   j=1,nmax
         gp(j)=g(i,j)+wi1(j)*dtem+ar*t1*exp(-t(i)/T1)*
     $        (1.-exp(-dtem/T1))*
     $        xinf(j)/sg+br*t2*exp(-t(i)/t2)*
     $        (1.-exp(-dtem/t2))*xinf(j)/sg
         gasp=gasp+gp(j)
         ind(j)=1
      enddo
 111  continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 402  continue
      HGA=0.
      NID=0
      DO  L=1,NMAX
         HGA=HGA+GP(L)
         NID=NID+IND(L)
      enddo

      GASP=HGA
      DO  L=1,NMAX
         XP(L)=GP(L)/GASP
      enddo

      IF(NID.EQ.NMAX) GO TO 333
      KDELTA=KDELTA+1
      IF(KDELTA-KDEMAX)103,103,250
 103  GO TO 200
 250  DTEM=DTEM/2.5
      KRIP=KRIP+1
      IF(KRIP.GT.14) GO TO 500
      GO TO 301
 333  CONTINUE
C Trovata una soluzione.

      ISTAR=I
      I=I+1
      IF(I+1.GT.IMAX) GO TO 400

      GO(I)=GASP
      DO  L=1,NMAX
         G(I,L)=GP(L)
         X(I,L)=XP(L)
      enddo

      T(I)=TN
      rate(i)=rap
      FALL(I)=FAP
      RSF=RATE(I)*GO(I)**CK
      SRT(I)=SRP
      ZETA=0.
      L1=5
      DO 320 LZ=L1,NMAX
 320  ZETA=ZETA+X(I,LZ)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C G1 e G2 sono la SFR al tempo I ed al tempo I-1.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      G1=GO(I)**CK*RATE(I)
      G2=GO(I-1)**CK*RATE(I-1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Con SS viene indicata la frazione di massa totale in stelle forma-
C tesi nell' intervallo DTEM.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SS=0.5*(G1+G2)*DTEM
      BO(I)=BO(I-1)+SS
      STARS=0.
      DWARF=0.
      SMSS=0.
C Area solar neighbourhood.
      AA=1.00E+5
C Area bulge.
C      AA=1.256E+4
C      MI=-0.9*ALOG10(ZETA)-1.7
C      DO 77 K=1,KMAX
C      C=AMU(K)-MI
C      IF(C.GT.0.) GO TO 79
C 77   CONTINUE
C 79   K1=K
C      DO 999 LK=K1,KMAX
      DO 999 LK=2,KMAX
      A=AMU(LK-1)
      B=AMU(LK)
      CALL SARA(A,B,SS,BROC,SOM)
      SM(LK-1)=SOM
      SSM(LK-1)=SM(LK-1)*SG*AA
      SMSS= SMSS+SSM(LK-1)
 999  STARS=STARS+SOM
      DO 888 LK=2,KMAX
      A=AMU(LK-1)
      B=AMU(LK)
      CALL SARA(A,B,SS,BROC,SOM)
     
      ENZO(LK-1)=BROC
      ENZO(LK-1)=ENZO(LK-1)*SG*AA

 888  DWARF=DWARF+BROC




C La normalizzazione dei rapporti [el/Fe] e [el/H] e' effettuata
C sui valori osservativi (Anders & Grevesse 1989).
C
C      Normalisation to observed solar abundances
C      is taken from Anders and Grevesse (1989)
C    


       SOF=0.878
       SFH=-2.743               

c     ASPLUND 2005
        sfh=-2.8018
       SMH=-3.135
c asplund 2005       
       smh=-3.089
       
       SCH=-2.365
       SOH=-1.865
       SNH=-2.801
       SSH=-2.994
       ssh=-3.0428
       SNEH=-2.637
       SCU=-3.179
       SZN=-2.784
       SNIF=-0.907
       SKRF=-5.584

c --------------------------------------

       SNEF=-0.106
       SMF=-0.392
       SSF=-0.252
       SCF=0.378
       SNF=-0.058
       SZF=-0.483
       SHEF=2.33
       SC13F=-1.54
       SCA=-1.3114
       STIF=-2.626

         smf=-0.2872
         ssf=-0.2412

c----Presi da Asplund!  2004--

       SLIF = -5.106 
       SKF  = -2.546
       SSCF = -4.505
       SVF  = -3.521
       SCRF = -1.852
       SMNF = -1.988
       SCOF = -2.568
       SBAF = -4.900
       SEUF = -6.526
       SLAF = -5.904
       SSrF = -4.374
       SYF  = -5.067
       SZRF = -4.667
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SNR = -0.95 
       SNAH = -4.204
       SALH = -3.975
       SNAF = SNAH-SFH
       SALF = SALH-SFH
c-------------------------
C      Ratio X/H
c-------------------------      

       AFH  =ALOG10(X(I,15)/X(I,1))-SFH
       AMH  =ALOG10(X(I,11)/X(I,1))-SMH
       ACH  =ALOG10(X(I,5)/X(I,1))-SCH
       AOH  =ALOG10(X(I,6)/X(I,1))-SOH
       ANH  =ALOG10(X(I,7)/X(I,1))-SNH
       ASH  =ALOG10(X(I,12)/X(I,1))-SSH
       ANEH =ALOG10(X(I,10)/X(I,1))-SNEH


c-------------------------
C      Ratio X/H
c-------------------------      

       AHEF =ALOG10(X(I,4)/X(I,15))-SHEF
       ACF  =ALOG10(X(I,5)/X(I,15))-SCF
       AOF  =ALOG10(X(I,6)/X(I,15))-SOF
       ANF  =ALOG10(X(I,7)/X(I,15))-SNF
       AC13F=ALOG10(X(I,8)/X(I,15))-SC13F
       ANR = ALOG10(x(I,9)/X(i,15))-SNR
C      Neutron rich
       ANEF =ALOG10(X(I,10)/X(I,15))-SNEF
       AMF  =ALOG10(X(I,11)/X(I,15))-SMF
       ASIF  =ALOG10(X(I,12)/X(I,15))-SSF
       ASF  =ALOG10(X(I,13)/X(I,15))-SZF
       ACAF =ALOG10(X(I,14)/X(I,15))-SCA
       ACUF =ALOG10(X(I,16)/X(I,15))-SCU
       AZNF =ALOG10(X(I,17)/X(I,15))-SZN
       ANIF =ALOG10(X(I,18)/X(I,15))-SNIF
       AKRF =ALOG10(X(I,19)/X(I,15))-SKRF
       ALIF =ALOG10(X(I,20)/X(I,15))-SLIF 
       AN15F=ALOG10(X(I,21)/X(I,15))-SN15F
       AKF  =ALOG10(X(I,22)/X(I,15))-SKF  
       ASCF =ALOG10(X(I,23)/X(I,15))-SSCF 
       ATIF =ALOG10(X(I,24)/X(I,15))-STIF 
       AVF  =ALOG10(X(I,25)/X(I,15))-SVF  
       ACRF =ALOG10(X(I,26)/X(I,15))-SCRF 
       AMNF =ALOG10(X(I,27)/X(I,15))-SMNF 
       ACOF =ALOG10(X(I,28)/X(I,15))-SCOF 
       ABAF =ALOG10((X(I,29)+X(I,30))/X(I,15))-SBAF
       AEUF =ALOG10(X(I,32)/X(I,15))-SEUF
       ALAF =ALOG10(X(I,33)/X(I,15))-SLAF
       ASrF =ALOG10(X(I,31)/X(I,15))-SSrF
       AYF =ALOG10(X(I,34)/X(I,15))-SYF
       AZrF =ALOG10(X(I,35)/X(I,15))-SZrF
         ANAH=ALOG10(X(I,36)/X(I,1))-SNAH
         ALH =ALOG10(X(I,37)/X(I,1))-SALH
         ANaF=ANaH-AFH
         ALF=ALH-AFH

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C         Non normalized abundance in number of atoms
C
       GRZ=ALOG10(ZETA)+12
       GRD  =ALOG10((X(I,2)/X(I,1))/2)+12
       GRHE3=ALOG10((X(I,3)/X(I,1))/3)+12
       GRHE =ALOG10((X(I,4)/X(I,1))/4)+12
       GRC  =ALOG10((X(I,5)/X(I,1))/12)+12
       GRO  =ALOG10((X(I,6)/X(I,1))/16)+12
       GRN  =ALOG10((X(I,7)/X(I,1))/14)+12
       GRC13=ALOG10((X(I,5)/X(I,1))/13)+12
ccc    Neutron rich
       GRNE =ALOG10((X(I,10)/X(I,1))/20)+12
       GRMG =ALOG10((X(I,11)/X(I,1))/24)+12          
       GRSI =ALOG10((X(I,12)/X(I,1))/28)+12
       GRS  =ALOG10((X(I,13)/X(I,1))/32)+12
       GRCA =ALOG10((X(I,14)/X(I,1))/40)+12  
       GRFE =ALOG10((X(I,15)/X(I,1))/56)+12
       GRCU =ALOG10((X(I,16)/X(I,1))/63.56)+12
       GRZN =ALOG10((X(I,17)/X(I,1))/65)+12
       GRNI =ALOG10((X(I,18)/X(I,1))/58.69)+12
       GRKR =ALOG10((X(I,19)/X(I,1))/83.80)+12
       GRLI =ALOG10((X(I,20)/X(I,1))/6.94)+12   
       GRN15=ALOG10((X(I,21)/X(I,1))/15.00)+12   
       GRK  =ALOG10((X(I,22)/X(I,1))/39.10)+12   
       GRSC =ALOG10((X(I,23)/X(I,1))/44.96)+12   
       GRTI =ALOG10((X(I,24)/X(I,1))/47.87)+12   
       GRV  =ALOG10((X(I,25)/X(I,1))/50.94)+12   
       GRCR =ALOG10((X(I,26)/X(I,1))/52.00)+12   
       GRMN =ALOG10((X(I,27)/X(I,1))/54.94)+12
       GRCO =ALOG10((X(I,28)/X(I,1))/58.93)+12   
       GRBA =ALOG10(( (X(I,29)+(X(I,30))) /X(I,1))/137.33)+12
       GREU =ALOG10((X(I,32)/X(I,1))/151.96)+12
       GRLA =ALOG10((X(I,33)/X(I,1))/138.90)+12    

       GRSR =ALOG10((X(I,31)/X(I,1))/87.62)+12    
       GRY =ALOG10((X(I,34)/X(I,1))/88.90)+12    
       GRZR =ALOG10((X(I,35)/X(I,1))/91.22)+12    
       GRNa =ALOG10((X(I,36)/X(I,1))/22.99)+12    
       GRAl =ALOG10((X(I,37)/X(I,1))/26.98)+12    




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

           
       RSO=ALOG10((X(I,13)/32)/(X(I,6)/16))
       RNO=ALOG10((X(I,7)/14)/(X(I,6)/16))
       RCO=ALOG10((X(I,5)/12)/(X(I,6)/16))
       RNEO=ALOG10((X(I,10)/20)/(X(I,6)/16))
C
C      new definitions for output file 17. 
C

       xofe= ALOG10( (X(I,6)/16)/(X(I,15)/56))
       xmgfe=ALOG10((X(I,11)/24)/(X(I,15)/56))
       xsife=ALOG10((X(I,12)/28)/(X(I,15)/56))
       xcafe=ALOG10((X(I,14)/40)/(X(I,15)/56))
       xznfe=ALOG10((X(I,17)/65)/(X(I,15)/56))
       xnife=ALOG10((X(I,18)/60)/(X(I,15)/56))
       xkpfe=ALOG10((X(I,22)/39)/(X(I,15)/56))
       xscfe=ALOG10((X(I,23)/45)/(X(I,15)/56))
       xtife=ALOG10((X(I,24)/48)/(X(I,15)/56))
       xvafe=ALOG10((X(I,25)/51)/(X(I,15)/56))
       xcrfe=ALOG10((X(I,26)/52)/(X(I,15)/56))
       xmnfe=ALOG10((X(I,27)/55)/(X(I,15)/56))
       xcofe=ALOG10((X(I,28)/59)/(X(I,15)/56))



C
C     +-------------------------------------------+
C     |                                           |
C     | Non normalized Abundances ratios          |
C     | for the LP UVES observed elements         |
C     | this will permit to choose between        |
C     | normalization with the computed values    |
C     | or the solar observed values              |
C     +-------------------------------------------+


c       WRITE(69,6969)t(i),grfe,xofe,xmgfe,xsife,
c     $ xcafe,xznfe,xnife,xkpfe,xscfe,xtife,xvafe,xcrfe,xmnfe,xcofe     
c
 6969  FORMAT(15(1x,1e12.5))


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Aggiunte Gabriele
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 6999 FORMAT(34 (1x,1e12.5))
 6929 FORMAT(35 (1x,1e12.5))

c!!   tolto elio messo neutron rich non normalizzati al solare!

      
                 aspi=log10((x(i,11)+x(i,12))/x(i,15))
     $      -log10(10**smf+10**ssf)

c         write(62,6999) T(i),ZETA,RSF,AFH ,aspi,ALIF,ACF,AC13F,AOF,ANF,
c     $        ANEF,AMF,ASIF,ASF,ACAF,ACUF,AZNF,ANIF,AKRF,AKF,ASCF,ATIF,
c     $        AVF,ACRF,AMNF,ACOF,ABAF,AEUF,ALAF,ASrF,AYF,AZrF, 
c     $        ANaF,AlF
         
c         write(67,6929) T(i),ZETA,RSF,GRFE ,GRHE,GRLI,GRC,GRC13,GRO,GRN,
c     $        GRNE,GRMG,GRSI,GRS,GRCA,GRCU,GRZN,GRNI,GRKR,GRK,GRSC,GRTI,
c     $        GRV,GRCR,GRMN,GRCO,GRBA,GREU,GRLA,GRSR,GRY,GRZR,GRNa,GRAl,
c     $        go(i)

c          write(67,*) T(i),GRFE ,GRMG,GRSI
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C Arquivo output fis.res - dens.gas,taxa infall,dens.total,SNII,SNI

c          write(80,2126)t(i),go(i)*SG,fall(i),srt(i),asnii*SG,gamma*SG
c     $  ,rsf*SG


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 2126 FORMAT(1X,1E12.5,1X,1E12.5,1X,1E12.5,1X,1E12.5,1X,
     $1E12.5,1X,1E12.5,1X,1E12.5)


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Solo per graficare in funzione del tempo che massa  stellare inizia a morire!
c        gigi=t(i)-dtem/2-tnew
c        gigi=assa(gigi)
c        write(57,*) AFH,gigi

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Aggiunte per Bario - Europio!
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c       write(11,6990) t(i),AFH,x(i,15),x(i,29),x(i,30),x(i,31),x(i,32),
c     $  X(i,33)
 6990  FORMAT(11 (1x,1e12.5))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Rame su ferro e zinco su ferro.
C 1    FORMAT(2X,I5)
C      IF (T(I).LT.14.0) GO TO 6662

c          write(88,*) i, t(i),fall(i),R


c      WRITE(60,15) DTEM,KDELTA
c      WRITE(60,2010) I,T(I),GO(I),ZETA,RSF,FALL(I),STARS*SG,SRT(I),
c     $(X(I,L),L=1,NMAX)


 15   FORMAT(20X,'DTEM=',1E15.5,10X,'KDELTA=',1I5)

 2010 FORMAT(2X,I5,2X,'T=',1F10.5,5X,'G=',1E12.5,5X,'Z=',1E12.5/
     $     2X,'RSF=',1E12.5,'FALL=',1E12.5,'STARS=',1E12.5,'SRT=',
     $     1E12.5/
     $     1X,7E12.5/1X,7E12.5/1X,7E12.5/,1X,7E12.5/1X,7E12.5/1X,2E12.5)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      STA(I)=DWARF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C File per ottenere la distribuzione differenziale delle nane G in
C funzione della metallicita`.

c      WRITE(66,2355) T(I),AFH,DWARF*SG,AOH
 2355 FORMAT(7(1X,1E12.5))
        TTspi(i)=T(i)
      aspiv(i)=aspi
      AMFv(i)=AMF
      AFEHv(i)=AFH
     
c      write(65,2356) T(i),ZETA,AFH ,aspi,AMF,ASIF,DWARF*SG,go(i)*SG
c     $,fall(i),srt(i),asnii*SG,gamma*SG
c     $     ,rsf*SG
 2356 FORMAT(13(1X,1E12.5))     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Definizione di ASP
      IF(I.EQ.1) THEN
      ASP=X(I,3)+3./2.*X(I,2)
      ELSE
      ASP=X(I-1,3)+3./2.*X(I-1,2)
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c



 

      VDTEM=DTEM
      TP=TN
 1010 CONTINUE
      ITER=0
      IELD=0

c!!      if (T(i).eq.5.0) then
c!!         open(56,file=modello//'.5.gr1',status='unknown')
c!!         do
c!!            read(56,*,iostat=ios) a
c!!            if (ios.ne.0) exit
c!!         enddo
c!!         
c!!         write(56,6919) T(i),ZETA,RSF,GRFE ,GRHE,GRLI,GRC,GRC13,GRO,
c!!     $        GRN,GRNE,GRMG,GRSI,GRS,GRCA,GRCU,GRZN,GRNI,GRKR,GRK,
c!!     $        GRSC,GRTI,GRV,GRCR,GRMN,GRCO,GRBA,GREU,GRLA,GRSR,GRY,
c!!     $        GRZR, go(i),r
c!!         close(56)
c!!         
c!!         if (r.eq.8.) then
c!!            open(55,file=modello//'.5.8.gr1',status='unknown')
c!!            do j=1,10
c!!                  write(55,6919) T(i),ZETA,RSF,GRFE ,GRHE,GRLI,GRC,
c!!     $              GRC13,GRO,GRN,
c!!     &              GRNE,GRMG,GRSI,GRS,GRCA,GRCU,GRZN,GRNI,GRKR,
c!!     $              GRK,GRSC,GRTI,GRV,   
c!!     &              GRCR,GRMN,GRCO,GRBA,GREU,GRLA,GRSR,GRY,GRZR,go(i),r
c!!               enddo
c!!               close(55)
c!!            endif
c!!            
c!!         endif
c!!

      GO TO 300
 600  CONTINUE

 610  FORMAT(10X,'SUPERATA ETA MASSIMA',1E15.5//)

      GO TO 1000

 400  CONTINUE
c 400  WRITE(60,410) IMAX

 410  FORMAT(10X,'SUPERA LA DIMENSIONE MAX DELLE MATRICI',1I10/)
      GO TO 1000
 500  CONTINUE
c 500  WRITE(60,510) KRIP
 510  FORMAT(10X,'NON TROVA SOLUZIONE DOPO 12 DIV DEL PASSO TEM',
     $1I10/)

 1000 CONTINUE

 6919 FORMAT(34 (1x,1e12.5))

c!!      open(57,file=modello //'.14.gr1',status='unknown')
c!!      do
c!!!      read(57,*,iostat=ios) a
c!!!      if (ios.ne.0) exit
c!!!      enddo
c!!!      
c!!!      write(57,6919) T(i),ZETA,RSF,GRFE ,GRHE,GRLI,GRC,GRC13,GRO,GRN,
c!!!     &     GRNE,GRMG,GRSI,GRS,GRCA,GRCU,GRZN,GRNI,GRKR,GRK,GRSC,GRTI,GRV,   
c!!!     &     GRCR,GRMN,GRCO,GRBA,GREU,GRLA,GRSR,GRY,GRZR,go(i),r
c!!!      close(57)
c!!!
c!!!      if (r.eq.8.) then
c!!!         open(59,file=modello//'14.8.gr1',status='unknown')
c!!!         do j=1,10
c!!!            write(59,6919) T(i),ZETA,RSF,GRFE ,GRHE,GRLI,GRC,GRC13,GRO,
c!!!     $           GRN,
c!!!     &     GRNE,GRMG,GRSI,GRS,GRCA,GRCU,GRZN,GRNI,GRKR,GRK,GRSC,GRTI,GRV,   
c!!!     &     GRCR,GRMN,GRCO,GRBA,GREU,GRLA,GRSR,GRY,GRZR,go(i),r
c!!!         enddo
c!!!         close(59)
c!!!      endif
c!!!


      RETURN
      END




      FUNCTION TINT(HM)
      PARAMETER (NMAX=37) 
      DIMENSION GO(3000),G(3000,NMAX),T(3000),XSTIN(NMAX)
      DIMENSION X(3000,NMAX),BO(3000),XCAL(NMAX),XTIN(NMAX),XSTAR(NMAX)
      DIMENSION COSN1(50),COSN2(50)
      COMMON/CUATRO/COSN1,COSN2           

      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/STAR/BO,TT1,ISTAR
      I=ISTAR+1
      TM=TAU(HM)
      IF((TT1-TM).GT.0.) GO TO 10
      TINT=BO(I)
 10   CONTINUE
      TIN=TT1-TM
      DO 80 J=1,I
      IF(T(J)-TIN) 80,90,100
 80   CONTINUE
      HBO=BO(J+1)-BO(J)
      HTO=T(J+1)-T(J)
      HT=T(J+1)-TIN
      VB=BO(J+1)-HBO/HTO*HT
      GO TO 101
 90   TINT=BO(I)-BO(J)
      RETURN
 100  CONTINUE
      HBO=BO(J)-BO(J-1)
      HTO=T(J)-T(J-1)
      IF(HTO.EQ.0.) GO TO 121
      HT=T(J)-TIN
      VB=BO(J)-HBO/HTO*HT
      GO TO 101
 121  VB=0.
 101  CONTINUE
      TINT=BO(I)-VB
      RETURN
      END


      SUBROUTINE SARA(XL,XU,SS,BROC,SOMMA)

C XL MASSA ESTREMO INFERIORE
C XU MASSA ESTREMO SUPERIORE
      COMMON/SUP/AO1,BO1,CO1
      IF(XU.NE.XL) GO TO 10
      SOMMA=0.
      RETURN
 10   CONTINUE
      HH=XU-XL
      H1=XL
      H2=XU
      F1=PSI(H1)*TINT(H1)
      F2=PSI(H2)*TINT(H2)
      SOMMA=(F1+F2)*HH*0.5
      B1=SS*PNU(H1)
      B2=SS*PNU(H2)
      BROC=(B1+B2)*HH*0.5
      RETURN
      END


            subroutine tempo
        PARAMETER (NMAX=37) 
      REAL MI,MS,MU  
      DIMENSION XINF(nmax),FALL(3000),FA(nmax)
      DIMENSION VN(NMAX),GRAT(3000)
      DIMENSION IND(NMAX)
      DIMENSION TW(NMAX),XCAL(NMAX),XTIN(NMAX),XSTAR(NMAX),XSTIN(NMAX)
      DIMENSION G(3000,NMAX),GO(3000),T(3000),X(3000,NMAX)
      DIMENSION SRT(3000),SMS(3000),APUNT(3000),RATE(3000),R1(NMAX)
      DIMENSION COSN1(50),COSN2(50)
      COMMON/CUATRO/COSN1,COSN2           
      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     $GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     $KMAX,I,ITER,IND,KDEMAX
      COMMON/DIS/EXPO,T1,T2,ETA,RO,RS,AR,BR,RI,RD,SD,SN,AO,
     $R,SMO,SMC,R1,VMR,RN,VNSOL,SRT,ARO,BRO,SMS,GAS,APUNT,
     $RATE,RIEST,RSCA,AC,TW,VN,THALO
      COMMON/MAO/IVO,IPER
      COMMON/PIP/GRAT,DTGRAT
      COMMON/RAF/TAUIN,TNEW,gr
    
     
      if(i-1) 20,20,30
 20   dtem=vdtem
      return
 30   continue
      tdif=t(i)-t(i-1)
      tt1=2.*vdtem
      do 40 l=1,nmax
      dg=g(i,l)-g(i-1,l)
      dg=abs(dg)
      if(dg.lt.1.e-10) go to 50
      t5=tol*tdif*g(i,l)/dg
      go to 60
 50   t5=0.99999*tt1
 60   t6=1.5*tdif
      t5 =amin1(t5,t6)
 40   continue
      if(t(i).lt.4.e-03) go to 110
      tt=t(i)
      to=t(i-1)
      sr1=ar*t1*(1.-exp(-to/t1))+br*t2*(1.-exp(-(to-tauin)/t2))+gas
      sr2=ar*t1*(1.-exp(-tt/t1))+br*t2*(1.-exp(-(tt-tauin)/t2))+gas
c      sr1=grat(i-1)*rate(i-1)*go(i-1)**ck
c      sr2=grat(i)*rate(i)*go(i)**ck
      dsr=sr2-sr1
      if(dsr.eq.0.) go to 120
      if(sr2.eq.0.) go to 120
      ds=dsr/sr2
      ds=abs(ds)
      if(ds.lt.5.e-02) go to 120
      b=1.-exp(-tt/t1)
      ts=-t1*alog(1.-5.e-02*b)
      t5=amin1(ts,t5)
 110  continue
      t5=2.e-05
C QUESTA SCHEDA VALE PER I VECCHI TEMPI DI VITA.
 120  dtem=amin1(t5,0.01)
      return
      end


      SUBROUTINE TEMPO2
      PARAMETER (NMAX=37) 
      REAL MI,MS,MU  
      DIMENSION XINF(nmax),FALL(3000),FA(nmax)
      DIMENSION VN(NMAX),GRAT(3000)
      DIMENSION IND(NMAX)
      DIMENSION TW(NMAX),XCAL(NMAX),XTIN(NMAX),XSTAR(NMAX),XSTIN(NMAX)
      DIMENSION G(3000,NMAX),GO(3000),T(3000),X(3000,NMAX)
      DIMENSION SRT(3000),SMS(3000),APUNT(3000),RATE(3000),R1(NMAX)
      DIMENSION COSN1(50),COSN2(50)
      COMMON/CUATRO/COSN1,COSN2           
      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     $GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     $KMAX,I,ITER,IND,KDEMAX
      COMMON/DIS/EXPO,T1,T2,ETA,RO,RS,AR,BR,RI,RD,SD,SN,AO,
     $R,SMO,SMC,R1,VMR,RN,VNSOL,SRT,ARO,BRO,SMS,GAS,APUNT,
     $RATE,RIEST,RSCA,AC,TW,VN,THALO
      COMMON/MAO/IVO,IPER
      COMMON/PIP/GRAT,DTGRAT
      COMMON/RAF/TAUIN,TNEW,gr
    
 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!                                            ! 
c! first step of   0.05 Gyr                   !
c!                                            !  
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (I.eq.1) then 

         DTEM=0.05
         
         RETURN
         
      endif                                  
                                             
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c    
c!
c! define the temporal step done before
c!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      
      TDIF=T(I)-T(I-1)


C!
c! define a step of 2 times  vdtem (2* 0.2e-4)
c!      

      TT1=2.*VDTEM

   
      SFR=RATE(I-1)*GO(I-1)**CK

c! when SRT grows over the  limit of star formation    (4.0)
c! take a small step (2.e-3) until SRT pass the limit  (4.5)
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!
c!!!!! when born the first stars it takes a smaller step 1.e-4 
c!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if (SRT(I).gt.4.and.SRT(I).lt.4.5) then
         
         t5=2.e-3

         if (SFR.ne.0) then
            
            T5=1.e-4

         endif
      
      else



c! if some metal grows up more then 1.e-10 in a temporal step
c! takes a temporal step of  
c
c!      2*vdtem
c
c! else takes a temporal step of 
c
c!      TOL*TDIF*G(I,L)/DG
c
c! then compare this temporal step with 
c
c!      1.5*TDIF 
c
c! and choose the smaller. 
c
c! It never takes temporal step greater then 0.2 Gyr 
c!
     

         DO  L=1,17
            
            DG=G(I,L)-G(I-1,L)
            
            DG=ABS(DG)
            
            IF(DG.LT.1.E-10) then
               
               T5=0.99999*TT1          
               
            else
               
               T5=TOL*TDIF*G(I,L)/DG
               
            endif
            
            
            T6=1.5*TDIF
            
            if (T6.gt.0.2) then
               
               T6=0.2 
               
            endif
            
            
            T5 =AMIN1(T5,T6)
            
         enddo
         
c!!!         
c!!! Quick Solution to the problem of no star
c!!! formation in the last 2-3 Gyr
c!!!        
         
         IF(T(I).GE.10.0) then
            
            T5=0.1
            
            DTEM=T5
            
            return
               
         endif
         
         
         
         
         TT=T(I)
         
         TO=T(I-1)
      
      
         SR1=AR*T1*(1.-EXP(-TO/T1))+BR*T2*(1.-EXP(-(TO-TAUIN)/T2))+GAS
         SR2=AR*T1*(1.-EXP(-TT/T1))+BR*T2*(1.-EXP(-(TT-TAUIN)/T2))+GAS
      
c         SR1=RATE(I-1)*GO(I-1)**CK
c         SR2=RATE(I)*GO(I)**CK
      
         DSR=(SR2-SR1)
      
         DS=DSR/SR2
         
         DS=ABS(DS)



         IF(DS.lt.5.E-02) then

            DTEM=T5
               
         else 
         
            B=1.-EXP(-TT/T1)
         
            TS=-T1*ALOG(1.-5.E-02*B)
         
            T5=AMIN1(TS,T5)
             
         endif
      
      endif


      IF ((GR.eq.0).and.(T(I).gt.4.80)) then
         GR=1
         T5=5.0-T(I)

      endif


         
      DTEM=T5
         
      RETURN
         
      END


      FUNCTION FII(XL,XU)
      PARAMETER (NMAX=37)  
      DIMENSION RR1(NMAX),RR2(NMAX),RR3(NMAX),RR4(NMAX)
      COMMON/SUP/AO1,BO1,CO1
      COMMON/SN/RR1,RR2,RR3,RR4,COST
      IF(XL.NE.XU) GO TO 10
      FII=0.
      RETURN
 10   CONTINUE
      A=-1.85
      B=-0.6545
      C=-3.7
      D=-2.51
      E=-1.78
      F=-0.86
      BB=1.0
      CC=1.35
      DD=0.77
      EE=0.17
      FF=-0.94
      HH=XU-XL
      H1=ASSA(XL)
      H2=ASSA(XU)
      IF(H1.LE.1.3) THEN
      F1= PNU(H1)*H1/(B*10**(B*ALOG10(H1)+BB))
      ELSE IF(H1.LE.3.) THEN
      F1=PNU(H1)*H1/(C*10**(C*ALOG10(H1)+CC))
      ELSE IF(H1.LE.7.) THEN
      F1= PNU(H1)*H1/(D*10**(D*ALOG10(H1)+DD))
      ELSE IF(H1.LE.15.) THEN
      F1=PNU(H1)*H1/(E*10**(E*ALOG10(H1)+EE))
      ELSE IF(H1.LE.60.) THEN
      F1=PNU(H1)*H1/(F*10**(F*ALOG10(H1)+FF))
      ELSE
      F1=PNU(H1)/(1.2*A*H1**(A-1.))
      END IF
      IF(H2.LE.1.3) THEN
      F2=PNU(H2)*H2/(B*10**(B*ALOG10(H2)+BB))
      ELSE IF(H2.LE.3.) THEN
      F2=PNU(H2)*H2/(C*10**(C*ALOG10(H2)+CC))
      ELSE IF(H2.LE.7.) THEN
      F2=PNU(H2)*H2/(D*10**(D*ALOG10(H2)+DD))
      ELSE IF(H2.LE.15.) THEN
      F2=PNU(H2)*H2/(E*10**(E*ALOG10(H2)+EE))
      ELSE IF(H2.LE.60.) THEN
      F2=PNU(H2)*H2/(F*10**(F*ALOG10(H2)+FF))
      ELSE
      F2=PNU(H2)/(1.2*A*H2**(A-1.))
      END IF
      IF(H1.GE.8.AND.H2.LE.16.) GO TO 20
      FII=(F1+F2)*HH*0.5
      GO TO 15
 20   FII=(F1+F2)*HH*0.5*(1.-COST)
 15   CONTINUE
      RETURN
      END



      FUNCTION FI(XL,XU)
      PARAMETER (NMAX=37) 
      DIMENSION RR1(NMAX),RR2(NMAX),RR3(NMAX),RR4(NMAX)
      COMMON/SN/RR1,RR2,RR3,RR4,COST
      IF(XU.NE.XL) GO TO 10
      FI=0.
      RETURN
 10   CONTINUE
      A=-1.85
      B=-0.6545
      C=-3.7
      D=-2.51
      E=-1.78
      F=-0.86
      BB=1.0
      CC=1.35
      DD=0.77
      EE=0.17
      FF=-0.94
      HH=XU-XL
      H1=ASSA(XL)
      H2=ASSA(XU)
      IF(H1.LE.1.3) THEN
      F1= PSI(H1)*H1/(B*10**(B*ALOG10(H1)+BB))
      ELSE IF(H1.LE.3.) THEN
      F1=PSI(H1)*H1/(C*10**(C*ALOG10(H1)+CC))
      ELSE IF(H1.LE.7.) THEN
      F1= PSI(H1)*H1/(D*10**(D*ALOG10(H1)+DD))
      ELSE IF(H1.LE.15.) THEN
      F1=PSI(H1)*H1/(E*10**(E*ALOG10(H1)+EE))
      ELSE IF(H1.LE.60.) THEN
      F1=PSI(H1)*H1/(F*10**(F*ALOG10(H1)+FF))
      ELSE
      F1=PSI(H1)/(1.2*A*H1**(A-1.))
      END IF
      IF(H2.LE.1.3) THEN
      F2=PSI(H2)*H2/(B*10**(B*ALOG10(H2)+BB))
      ELSE IF(H2.LE.3.) THEN
      F2=PSI(H2)*H2/(C*10**(C*ALOG10(H2)+CC))
      ELSE IF(H2.LE.7.) THEN
      F2=PSI(H2)*H2/(D*10**(D*ALOG10(H2)+DD))
      ELSE IF(H2.LE.15.) THEN
      F2=PSI(H2)*H2/(E*10**(E*ALOG10(H2)+EE))
      ELSE IF(H2.LE.60.) THEN
      F2=PSI(H2)*H2/(F*10**(F*ALOG10(H2)+FF))
      ELSE
      F2=PSI(H2)/(1.2*A*H2**(A-1.))
      END IF
      IF(H1.GE.8.AND.H2.LE.16.) GO TO 20
      FI=(F1+F2)*HH*0.5
      GO TO 15
 20   FI=(F1+F2)*HH*0.5*(1.-COST)
 15   CONTINUE
      RETURN
      END

      SUBROUTINE SUPER
C
C     Subroutine to compute the rate of SNI
C     
C
      PARAMETER (NMAX=37) 
      include 'lettura.inc'
      include 'elementi.inc'
  
      DIMENSION XCAL(NMAX),XTIN(NMAX),XSTAR(NMAX)
      DIMENSION PARZW(NMAX),PRZW1(NMAX),PRZW2(NMAX),PRZW3(NMAX)
      DIMENSION FF(NMAX),FF1(NMAX),GG1(NMAX),FO(NMAX),FFO(NMAX)
      DIMENSION F2(NMAX),FF2(NMAX),F3(NMAX),FF3(NMAX)
      DIMENSION AMU(200),GP(NMAX),WI(NMAX),WI1(NMAX),WI2(NMAX)
      DIMENSION XMB(NMAX),XM2V(NMAX),XMU(NMAX),F(NMAX),G1(NMAX)
      DIMENSION F1(NMAX),RR1(NMAX),RR2(NMAX),RR4(NMAX),RR3(NMAX)
      DIMENSION GO(3000),T(3000),X(3000,NMAX),IND(NMAX),TM(3000)
      DIMENSION XSTIN(NMAX),RATE(3000),APUNT(3000),PAMU(200),XP(NMAX)
      DIMENSION TW(NMAX),VN(NMAX),R1(NMAX),G(3000,NMAX),QM(NMAX,NMAX)
      DIMENSION SMS(3000),SRT(3000),FA(NMAX),FALL(3000),XINF(NMAX)
      DIMENSION GRAT(3000)
      DIMENSION COSN1(50),COSN2(50)
      COMMON/CUATRO/COSN1,COSN2           
      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/DIE/A,B,PAMU,K2,K1,KB,KMAX2,K3,KB3
      COMMON/TRE/AMU,GP,XP,WI,PARZW,WI1,WI2,PRZW1,PRZW2,PRZW3
      COMMON/A/QM
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     1GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     2KMAX,I,ITER,IND,KDEMAX
      COMMON/DIS/EXPO,T1,T2,ETA,RO,RS,AR,BR,RI,RD,SD,SN,AO,
     1R,SMO,SMC,R1,VMR,RN,VNSOL,SRT,ARO,BRO,SMS,GAS,APUNT,
     2RATE,RIEST,RSCA,AC,TW,VN,THALO
      COMMON/FR/XINF,FALL,FA,SG
      COMMON/SN/RR1,RR2,RR3,RR4,COST
      COMMON/SUP/AO1,BO1,CO1
      COMMON/NEW/ASNII,GAMMA
      COMMON/PIP/GRAT,DTGRAT
      COMMON/RAF/TAUIN,TNEW,gr

C RICORDA DI CAMBIARE BO1 QUANDO VUOI SALPETER

      IF(A.GT.8.) GO TO 20

c      IF(T(I).LT.TAUIN) THEN
c      TAUIN=0.
c      ELSE
c      TAUIN=2.
c      END IF
      TMEDIO=T(I)+DTEM/2.
      GAMA=2.
      COST=0.05
c      COST=0.1
      X1=2.*A
      C=A+8.
      XUNF=AMAX1(3.,X1)
      DELT=(16.-XUNF)/12.
      DO 40 J1=1,11
      JJ=J1-1
 40   XMB(J1)=XUNF+JJ*DELT
      XMB(12)=16.
C      WRITE(60,9001) (XMB(KK),KK=1,12)
 9001 FORMAT(1X,5E12.5/)
      IF(B.LT.8.) GO TO 1
      IF(ITER.EQ.1) GO TO 20
      DO 70 JX=1,NMAX
      DO 120 L=1,12
      IF(XMB(L).LT.C) GO TO 777
      D=(XMB(L)-8.)/XMB(L)
      C1=A/XMB(L)
      YMIN=AMAX1(C1,D)
      GO TO 666
 777  YMIN=A/XMB(L)
 666  XLIM=ABS(0.5-YMIN)
      IF(XLIM.LE.1.E-10) GO TO 111
      DELTA=0.5-YMIN
      XMU(1)=YMIN
      XMU(2)=YMIN+0.01*DELTA
      XMU(3)=YMIN+0.02*DELTA
      XMU(4)=YMIN+0.05*DELTA
      XMU(5)=YMIN+0.1*DELTA
      XMU(6)=YMIN+0.2*DELTA
      XMU(7)=YMIN+0.3*DELTA
      XMU(8)=YMIN+0.4*DELTA
      XMU(9)=YMIN+0.6*DELTA
      XMU(10)=YMIN+0.8*DELTA
      XMU(11)=0.5
      DO 100 L2=1,11
      XM2V(L2)=XMU(L2)*XMB(L)
      H22=XM2V(L2)
      TM(L2)=ABS(TMEDIO-TAU(H22))
      IF(I.EQ.1) GO TO 96
      DO 18 KO=1,I
      IF(TM(L2)-T(KO))19,18,18
 18   CONTINUE
 19   CONTINUE
      ISUP=KO
      INF=KO-1
      GO TO 27
 96   ISUP=1
      GINF=GO(ISUP)
      DO 9 LL=1,NMAX
 9    XSTIN(LL)=X(ISUP,LL)
      TIS=T(ISUP)
      RAINF=RATE(ISUP)*GINF**CK
      GO TO 28
 27   CONTINUE
      DTMT=TM(L2)-T(INF)
      TIS=T(ISUP)
      TIN=T(INF)
      DT=TIS-TIN
      DG=GO(ISUP)-GO(INF)
      GINF=GO(INF)+DG/DT*DTMT
      DO 4 KK=1,NMAX
 4    XSTIN(KK)=X(INF,KK)+(X(ISUP,KK)-X(INF,KK))/DT*DTMT
c      SRTM= AR*T1*(1.-EXP(-TM(L2)/T1))+BR*T2*(1.-EXP(-TM(L2)/T2))+GAS
c      SMSM=ARO*T1*(1.-EXP(-TM(L2)/T1))+BRO*T2*(1.-EXP(-TM(L2)/T2))+GAS
      IF(TM(L2).LT.TAUIN) THEN
      br=0.
      bro=0.
      else
      br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
      bro=br
      end if
      SRTM=AR*T1*(1.-EXP(-TM(L2)/T1))+BR*T2*(1.-
     $EXP(-(TM(L2)-TAUIN)/T2))+GAS
      SMSM=ARO*T1*(1.-EXP(-TM(L2)/T1))+BRO*T2*(1.-
     $EXP(-(TM(L2)-TAUIN)/T2))+GAS
      GRATM=(GRAT(INF)+GRAT(ISUP))*0.5
      RAINF=GRATM*(SRTM/SMSM)**EXPO*(SG/SRTM)**(CK-1.)*GINF**CK
 28   CONTINUE
      SIJ=0.0
      H3=(1.-XMU(L2))/XMU(L2)*H22
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Modifica ZETIN
      ZETIN=0.
      DO 57 LLL=5,NMAX
 57   ZETIN=ZETIN+XSTIN(LLL)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 5 JS=1,NMAX
 5    SIJ=SIJ+QBIN(H3,JX,JS,ZETIN)*XSTIN(JS)
      HH=XMB(L)
      F(L2)=(XMU(L2)**GAMA)*SIJ*RAINF*PSI(HH)
      FF(L2)=(XMU(L2)**GAMA)*RAINF*PNU(HH)
 100  CONTINUE
      SAM=0.
      SUM=0.
      DO 110 L3=1,10
      DA=(F(L3)+F(L3+1))*(XMU(L3+1)-XMU(L3))/2.
      DAA=(FF(L3)+FF(L3+1))*(XMU(L3+1)-XMU(L3))/2.
      SAM=SAM+DAA
      SUM=SUM+DA
 110  CONTINUE
      F1(L)=SUM
      FF1(L)=SAM
      GO TO 120
 111  F1(L)=0.
      FF1(L)=0.
 120  CONTINUE
      SAM1=0.
      SUM1=0.
      DO 140 L4=1,11
      DA=(F1(L4)+F1(L4+1))*(XMB(L4+1)-XMB(L4))/2.
      DAA=(FF1(L4)+FF1(L4+1))*(XMB(L4+1)-XMB(L4))/2.
      SAM1=SAM1+DAA
      SUM1=SUM1+DA
 140  CONTINUE
      GAMMA=SAM1*COST*2.**(1.+GAMA)*(1.+GAMA)
      RR1(JX)=SUM1*COST*2.**(1.+GAMA)*(1.+GAMA)
C Si e` calcolato l' integrale per ogni elemento.
 70   CONTINUE
C      WRITE(60,3127) GAMMA
 3127 FORMAT(1X,'SNI=',1E12.5)
      GO TO 20
 1    CONTINUE
      XME=2.*B
      DELT1=(XME-XUNF)/7.
      DO 8790 J1=1,7
      JJ=J1-1
 8790 XMB(J1)=XUNF+JJ*DELT1
      XMB(8)=XME
      DELT2=(16.-XME)/4.
      DO 8690 J1=9,11
      JJ=J1-8
 8690 XMB(J1)=XME+JJ*DELT2
      XMB(12)=16.
C      IF(ITER.EQ.1) GO TO 45
C      WRITE(60,2453) (XMB(K),K=1,12)
 2453 FORMAT(1X,'XMB=',1E12.5)
      DO 71 JX=1,NMAX
      DO 130 L=1,7
      IF(XMB(L).LT.C) GO TO 555
      D=(XMB(L)-8.)/XMB(L)
      C1=A/XMB(L)
      YMIN=AMAX1(C1,D)
      GO TO 444
 555  YMIN=A/XMB(L)
 444  CONTINUE
      XLIM=ABS(0.5-YMIN)
      IF(XLIM.LE.1E-10) GO TO 161
      DELTA=0.5-YMIN
      XMU(1)=YMIN
      XMU(2)=YMIN+0.01*DELTA
      XMU(3)=YMIN+0.02*DELTA
      XMU(4)=YMIN+0.05*DELTA
      XMU(5)=YMIN+0.1*DELTA
      XMU(6)=YMIN+0.2*DELTA
      XMU(7)=YMIN+0.3*DELTA
      XMU(8)=YMIN+0.4*DELTA
      XMU(9)=YMIN+0.6*DELTA
      XMU(10)=YMIN+0.8*DELTA
      XMU(11)=0.5
      DO 150 N3=1,11
      XM2V(N3)=XMU(N3)*XMB(L)
      H22=XM2V(N3)
C      WRITE(60,6666) H22
 6666 FORMAT(1X,'H22=',1E12.5)
      TM(N3)=ABS(TMEDIO-TAU(H22))
      IF(I.EQ.1) GO TO 86
      DO 24 KO=1,I
      IF(TM(N3)-T(KO))25,24,24
 24   CONTINUE
 25   CONTINUE
      ISUP=KO
      INF=KO-1
      GO TO 26
 86   ISUP=1
      GINF=GO(ISUP)
      DO 3 LL=1,NMAX
 3    XSTIN(LL)=X(ISUP,LL)
      TIS=T(ISUP)
      RAINF=RATE(ISUP)*GINF**CK
      GO TO 38
 26   CONTINUE
      DTMT=TM(N3)-T(INF)
      TIS=T(ISUP)
      TIN=T(INF)
      DT=TIS-TIN
      DG=GO(ISUP)-GO(INF)
      GINF=GO(INF)+DG/DT*DTMT
      DO 17 KK=1,NMAX
 17   XSTIN(KK)=X(INF,KK)+(X(ISUP,KK)-X(INF,KK))/DT*DTMT
c      SRTM=AR*T1*(1.-EXP(-TM(N3)/T1))+BR*T2*(1.-EXP(-TM(N3)/T2))+GAS
c      SMSM=ARO*T1*(1.-EXP(-TM(N3)/T1))+BRO*T2*(1.-EXP(-TM(N3)/T2))+GAS
      IF(TM(N3).LT.TAUIN) THEN
      br=0.
      bro=0.
      else
      br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
      bro=br
      end if
      SRTM=AR*T1*(1.-EXP(-TM(N3)/T1))+BR*T2*(1.-
     $EXP(-(TM(N3)-TAUIN)/T2))+GAS
      SMSM=ARO*T1*(1.-EXP(-TM(N3)/T1))+BRO*T2*(1.-
     $EXP(-(TM(N3)-TAUIN)/T2))+GAS
      GRATM=(GRAT(INF)+GRAT(ISUP))*0.5
      RAINF=GRATM*(SRTM/SMSM)**EXPO*(SG/SRTM)**(CK-1.)*GINF**CK
 38   CONTINUE
      SIJ1=0.0
      H3=(1.-XMU(N3))/XMU(N3)*H22
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Modifica ZETIN
      ZETIN=0.
      DO 556 LKL=5,NMAX
 556  ZETIN=ZETIN+XSTIN(LKL)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 55 JS=1,NMAX
 55   SIJ1=SIJ1+QBIN(H3,JX,JS,ZETIN)*XSTIN(JS)
      HH=XMB(L)
      FO(N3)=(XMU(N3)**GAMA)*SIJ1*RAINF*PSI(HH)
      FFO(N3)=(XMU(N3)**GAMA)*RAINF*PNU(HH)
 150  CONTINUE
      SIM=0.
      SOM=0.
      DO 189 L3=1,10
      DA=(FO(L3)+FO(L3+1))*(XMU(L3+1)-XMU(L3))/2.
      DAA=(FFO(L3)+FFO(L3+1))*(XMU(L3+1)-XMU(L3))/2.
      SIM=SIM+DAA
      SOM=SOM+DA
 189  CONTINUE
      F2(L)=SOM
      FF2(L)=SIM
      GO TO 130
 161  F2(L)=0.
      FF2(L)=0.
 130  CONTINUE
      SIM1=0.
      SOM1=0.
      DO 180 L4=1,6
      DA=(F2(L4)+F2(L4+1))*(XMB(L4+1)-XMB(L4))/2.
      DAA=(FF2(L4)+FF2(L4+1))*(XMB(L4+1)-XMB(L4))/2.
      SIM1=SIM1+DAA
      SOM1=SOM1+DA
 180  CONTINUE
      GAM0=SIM1*COST*2.**(1.+GAMA)*(1.+GAMA)
      RR4(JX)=SOM1*COST*2.**(1.+GAMA)*(1.+GAMA)
C Il primo pezzo di integrale e` risolto...
 45   CONTINUE
      CC=B+8.
      DO 170 N=8,12
      IF(XMB(N).LT.CC) GO TO 222
      D=(XMB(N)-8.)/XMB(N)
      C1=B/XMB(N)
      YMIN=AMAX1(C1,D)
      GO TO 3331
 222  YMIN=B/XMB(N)
 3331 CONTINUE
      XLIM=ABS(0.5-YMIN)
      IF(XLIM.LE.1.E-10) GO TO 181
      DELTA=0.5-YMIN
      XMU(1)=YMIN
      XMU(2)=YMIN+0.01*DELTA
      XMU(3)=YMIN+0.02*DELTA
      XMU(4)=YMIN+0.05*DELTA
      XMU(5)=YMIN+0.1*DELTA
      XMU(6)=YMIN+0.2*DELTA
      XMU(7)=YMIN+0.3*DELTA
      XMU(8)=YMIN+0.4*DELTA
      XMU(9)=YMIN+0.6*DELTA
      XMU(10)=YMIN+0.8*DELTA
      XMU(11)=0.5
      DO 190 N4=1,11
      XM2V(N4)=XMU(N4)*XMB(N)
      H22=XM2V(N4)
      TM(N4)=ABS(TMEDIO-TAU(H22))
      TDIZ=T(I)
      TSTAR=TM(N4)
      GSTAR=GO(I)+(GASP-GO(I))/DTEM*(TSTAR-TDIZ)
      DO 6 L=1,NMAX
 6    XSTIN(L)=X(I,L)+(XP(L)-X(I,L))/DTEM*(TSTAR-TDIZ)
      IF(TSTAR.LT.TAUIN) THEN
      br=0.
      bro=0.
      else
      br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
      bro=br
      end if
      SSTAR=AR*T1*(1.-EXP(-TSTAR/T1))+BR*T2*
     $(1.-EXP(-(TSTAR-TAUIN)/T2))
     $+ GAS
      SMTAR=ARO*T1*(1.-EXP(-TSTAR/T1))+BRO*T2*
     $(1.-EXP(-(TSTAR-TAUIN)/T2))
     $+GAS
      RAINF=GRAT(I)*(SSTAR/SMTAR)**EXPO*(SG/SSTAR)**(CK-1.)*GSTAR**CK
      SIJ2=0.
      H3=(1.-XMU(N4))/XMU(N4)*H22
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Modifica ZETIN
      ZETIN=0.
      DO 557 LVL=5,NMAX
 557  ZETIN=ZETIN+XSTIN(LVL)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 7 JS=1,NMAX
 7    SIJ2=SIJ2+QBIN(H3,JX,JS,ZETIN)*XSTIN(JS)
      HH=XMB(N)
      GG1(N4)=(XMU(N4)**GAMA)*RAINF*PNU(HH)
      G1(N4)=(XMU(N4)**GAMA)*SIJ2*RAINF*PSI(HH)
 190  CONTINUE
      SEM1=0.
      SEM=0.
      DO 200 N5=1,10
      DA=(G1(N5)+G1(N5+1))*(XMU(N5+1)-XMU(N5))/2.
      DAA=(GG1(N5)+GG1(N5+1))*(XMU(N5+1)-XMU(N5))/2.
      SEM1=SEM1+DAA
      SEM=SEM+DA
 200  CONTINUE
      F3(N)=SEM
      FF3(N)=SEM1
      GO TO 170
 181  F3(N)=0.
      FF3(N)=0.
 170  CONTINUE
      SAM2=0.
      SUM2=0.
      DO 99 NN=8,11
      DAA=(FF3(NN)+FF3(NN+1))*(XMB(NN+1)-XMB(NN))/2.
      DA=(F3(NN)+F3(NN+1))*(XMB(NN+1)-XMB(NN))/2.
      SAM2=SAM2+DAA
      SUM2=SUM2+DA
 99   CONTINUE
      GAM2=SAM2*COST*2.**(1.+GAMA)*(1.+GAMA)
      GAMMA=GAM0+GAM2
      RR3(JX)=SUM2*COST*2.**(1.+GAMA)*(1.+GAMA)
      RR2(JX)=RR3(JX)+RR4(JX)
 71   CONTINUE
C      WRITE(60,8127) RR4,RR3
 8127 FORMAT(1X,'RR4=',6E12.5/1X,7E12.5/1X,'RR3=',6E12.5/1X,7E12.5)
C      WRITE(60,9994) GAMMA
 9994 FORMAT(1X,'SNI=',1E12.5)
 20   CONTINUE
C      WRITE(60,9111) RR1
 9111 FORMAT(1X,'RR1=',6E12.5/1X,7E12.5/)
      RETURN
      END





      FUNCTION Q(H,J1,J2,ZETAII)
C
C     compute the Q matrix for single stars
C
      PARAMETER (NMAX=37) 
C Si utilizzano gli yields stellari tabulati da van den Hoek & Groenewegen
C (1997) per low and intermediate mass stars; si utilizzano i modelli di
C nucleosintesi di Thielemann, Nomoto & Hashimoto (1996) per massive stars.
      include 'lettura.inc'
      include 'elementi.inc'

      REAL M1,M2
      character*80 titre

      DIMENSION AMLI(50),ALI0(50),ALI1(50),ALI2(50),ALI3(50),ALI4(50)
      DIMENSION PMASS(50),YLI0(50),YLI1(50),YLI2(50),YLI3(50),YLI4(50),
     $YLI5(50),YLI(50)
      REAL BID
      DIMENSION QM(NMAX,NMAX),XSSTAR(NMAX)
      DIMENSION GO(3000),G(3000,NMAX),T(3000),X(3000,NMAX),
     $XCAL(NMAX),XTIN(NMAX),XSTIN(NMAX),XSTAR(NMAX)
      COMMON/A/ QM
      COMMON/UNO/ GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/WOO/AMLI,ALI0,ALI1,ALI2,ALI3,ALI4
      COMMON/ZOO/PMASS,YLI0,YLI1,YLI2,YLI3,YLI4,YLI5,YLI
      COMMON/DATA/LET
c     !COMMON/LUCKY/IKU
      COMMON/Z/ZETA
      COMMON/M/NMAT,NNMAT
      COMMON/ELDS/IELD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Blocco COMMON per la nuova formulazione di ASP
      COMMON/HEL/ASP
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      M1=0.5
      M2=8.


 511  FORMAT(5E12.5)
 52   FORMAT(5E12.5)
 51   FORMAT(1I10)
 53   FORMAT(4E12.5)
 4567 FORMAT(F6.1,2X,5(F6.4,2X)) 
 1919 FORMAT(A80)
 777  FORMAT(A80)  
 9988 FORMAT(6E12.5)
 9989 FORMAT(F6.2,6E10.2)


c --------------------------------------
C Acquisizione delle masse espulse sotto forma dei diversi elementi chi-
C mici considerati (n.b.: per low and intermediate mass stars c'e` dipen-
C denza dalla metallicita` iniziale); acquisizione dei remnants stellari.
c --------------------------------------

      IF(LET.EQ.1) GO TO 501
c      write(*,*) let

      READ(5,51) NMAT
      READ(5,51) NNMAT
      DO 661 J=1,NMAT
         READ(5,511) HE1(J),C121(J),O161(J),AN141(J),C131(J)
         READ(5,52) ANE1(J),AMG1(J),ASI1(J),AFE1(J),S141(J)
         READ(5,53) C13S1(J),S321(J),CA1(J),RM1(J)
 661  CONTINUE
      DO 662 J=1,NMAT
         READ(5,511) HE2(J),C122(J),O162(J),AN142(J),C132(J)
         READ(5,52) ANE2(J),AMG2(J),ASI2(J),AFE2(J),S142(J)
         READ(5,53) C13S2(J),S322(J),CA2(J),RM2(J)
 662  CONTINUE
      DO 663 J=1,NMAT
         READ(5,511) HE3(J),C123(J),O163(J),AN143(J),C133(J)
         READ(5,52) ANE3(J),AMG3(J),ASI3(J),AFE3(J),S143(J)
         READ(5,53) C13S3(J),S323(J),CA3(J),RM3(J)
 663  CONTINUE
      DO 664 J=1,NMAT
         READ(5,511) HE4(J),C124(J),O164(J),AN144(J),C134(J)
         READ(5,52) ANE4(J),AMG4(J),ASI4(J),AFE4(J),S144(J)
         READ(5,53) C13S4(J),S324(J),CA4(J),RM4(J)
 664  CONTINUE
      DO 665 J=1,NMAT
         READ(5,511) HE5(J),C125(J),O165(J),AN145(J),C135(J)
         READ(5,52) ANE5(J),AMG5(J),ASI5(J),AFE5(J),S145(J)
         READ(5,53) C13S5(J),S325(J),CA5(J),RM5(J)
 665  CONTINUE
      DO 655 J=NMAT+1,NMAT+NNMAT
         READ(5,52) HE(J),C12(J),O16(J),AN14(J),C13(J)
         READ(5,52) ANE(J),AMG(J),ASI(J),AFE(J),S14(J)
         READ(5,53) C13S(J),S32(J),CA(J),RM(J)
 655  CONTINUE

 669  CONTINUE
       
       METDEP=1                 
       IF (METDEP.eq.1) then

C _____________________________________________________________________    
C     Metal dependent yields for O16 & Carbon 12
C     from Woosley and Weaver 95 + C12 from Maeder '92
     
          open(unit=19,file='WW95-O16.dat',status='old')
            READ(19,1919)TITRE
c            write(6,1919)TITRE
            DO IK=NMAT+1,NMAT+NNMAT
               READ(19,*)BID,(XXO16(IK,IL),IL=1,5)
     
c               write(6,*)BID,(XXO16(IK,IL),IL=1,5)
            END DO
            CLOSE(19)
c!!!
            open(unit=20,file='WW95-C12MM2.dat',status='old')
            READ(20,1919)TITRE
c            write(6,1919)TITRE
            DO IK=NMAT+1,NMAT+NNMAT
               READ(20,*)BID,(XXC12(IK,IL),IL=1,5)
c               write(6,*)BID,(XXC12(IK,IL),IL=1,5)
            END DO
            CLOSE(20)
         
      
c!!!
c!!!! Yields dipendenti dalla metallicità secondo prescrizioni
c!!!! Meynet 2002 (4 metallicità), variazioni da Hirschi et al.
c!!!! per metallicità 10^-8
c
c            
c!!!            open(unit=20,file='MM-Hirshi-C12.dat',status='old')
c!!!            READ(20,1919)TITRE
c!!!            
c!!!            DO IK=NMAT+1,NMAT+NNMAT
c!!!               READ(20,*)BID,(XXC12(IK,IL),IL=1,5)
c!!!               
c!!!               
c!!!            END DO
      
         end if


C     
C      End reading of Metal dependent yields
C________________________________________________________________________




         read(17,777) titre
c         write(6,777) titre
         read(17,777) titre
c         write(6,777) titre
 
         DO 659  J=NMAT+1,NMAT+NNMAT
            read(17,*)Xmas,ZN(J),KP(J),SC(J),TI(J),VANA(J),CR(J),
     $           Manga(J),CO(J),ani(j)
c            write(6,*)Xmas,ZN(J),KP(J),SC(J),TI(J),VANA(J),CR(J),
c     $           Manga(J), CO(J),ani(j)
 659     CONTINUE
         CLOSE(17)

c     --------------------------------------
C     Acquisizione masse stellari iniziali e (solo per massive stars) acqui-
C     sizione masse dell' helium core.
c     --------------------------------------
         
         DO 150 J=1,NMAT
            READ(5,53) AMM(J)            
 150     CONTINUE
         
         DO 787 J=NMAT+1,NMAT+NNMAT
            READ(5,53) AMM(J),AM(J)
 787     CONTINUE
         
C
C   Output for yields ratios
C___________________________________________________________________-


c!      open(unit=75,file='WW95-alpha.dat', status='unknown')
c
c!      write(75,*)'          WW95 original data '
c!      write(75,*)' mass     O      Mg       Si      Ca      Fe'
c
c!      DO 9876  JL=NMAT+1,NMAT+NNMAT
c!         WRITE(75,4567)AMM(JL),O16(JL),AMG(JL),ASI(JL),CA(JL),AFE(JL)
c! 9876 END DO
     

c --------------------------------------
C Acquisizione delle masse espulse sotto forma di litio7 per diversi va-
C lori della metallicita` iniziale da stelle con massa compresa tra 11.0
C 65 e 100 masse solari.
c --------------------------------------

      DO 1511 J=1,15
         READ(5,9988) AMLI(J),ALI0(J),ALI1(J),
     $        ALI2(J),ALI3(J),ALI4(J)
 1511 CONTINUE


c --------------------------------------
C Acquisizione degli yields di litio da stelle di AGB.
c --------------------------------------
      DO 1512 J=1,10
         READ(5,9989) PMASS(J),YLI0(J),YLI1(J),
     $        YLI2(J),YLI3(J),YLI4(J),YLI5(J)
 1512 CONTINUE
      CLOSE(5)

C-----------------------------------------
C     La rotazione ha l'effetto di aumentare
C     lo yield di C da stelle massicce:
C     TNH96 yields:
C     DO 7761 J=NMAT+NNMAT-2,NMAT+NNMAT
C     WW95 yields:
C     DO 7761 J=NMAT+NNMAT-6,NMAT+NNMAT
C     C12(J)=3.*C12(J)
C     7761    CONTINUE
C-----------------------------------------


      DO J=1,NMAT           
         ANE(J)=0.
         AMG(J)=0.
         ASI(J)=0.
         AFE(J)=0.
         S14(J)=0.
         C13S(J)=0.
         S32(J)=0.
         CA(J)=0.
         KP(J)=0.
         SC(J)=0.
         TI(J)=0.
         VANA(J)=0.
         CR(J)=0.
         Manga(J)=0.
         CO(J)=0.
      ENDDO

      LET=1

 501  CONTINUE

C     
C     Sceglie il corretto yield per l'ossigeno a seconda della metallicita' presente!!!
C
      
      METDEP=1

      IF (Metdep.eq.1) then     

         zetaIII=zetaII/0.012
         IF (ZETAIII.LE.1E-4) THEN
            E=(0.00000120-ZETAII)/(0.00000120-0.0)  
            F=(ZETAII-0.0)/(0.00000120-0.0)
            DO J=NMAT+1,NMAT+NNMAT
               O16(J)=E*XXO16(J,5)+F*xxO16(J,4)
               C12(J)=E*XXC12(J,5)+F*XXC12(J,4)
               ANA23(J)=XXNA23(J,5)
               AL27(J)=XXAL27(J,5)
            END DO
         ELSE IF (ZETAIII.GT.1E-4.AND.ZETAIII.LE.1E-2) THEN
            E=(0.000120-ZETAII)/(0.000120-0.0000012)  
            F=(ZETAII-0.0000012)/(0.000120-0.0000012)
            DO J=NMAT+1,NMAT+NNMAT
               O16(J)=E*XXO16(J,4)+F*xxO16(J,3) 
               C12(J)=E*XXC12(J,4)+F*XXC12(J,3)
               ANA23(J)=XXNA23(J,4)
               AL27(J)=XXAL27(J,4)
            END DO
         ELSE IF (ZETAIII.GT.1E-2.AND.ZETAIII.LE.1E-1) THEN
            E=(0.00120-ZETAII)/(0.00120-0.00012)  
            F=(ZETAII-0.00012)/(0.00120-0.00012)
            DO J=NMAT+1,NMAT+NNMAT
               O16(J)=E*XXO16(J,3)+F*xxO16(J,2)
               C12(J)=E*XXC12(J,3)+F*XXC12(J,2)
               ANA23(J)=XXNA23(J,3) 
               AL27(J)=XXAL27(J,3) 
            END DO
         ELSE IF (ZETAIII.GT.1E-1.AND.ZETAIII.LE.1) THEN
            E=(0.0120-ZETAII)/(0.0120-0.0012)  
            F=(ZETAII-0.0012)/(0.0120-0.0012)
            DO J=NMAT+1,NMAT+NNMAT
               O16(J)=E*XXO16(J,2)+F*xxO16(J,1)
               C12(J)=E*XXC12(J,2)+F*XXC12(J,1)
               ANA23(J)=XXNA23(J,2)
               AL27(J)=XXAL27(J,2)
            END DO
         ELSE IF (ZETAIII.GT.1.) THEN
            DO J=NMAT+1,NMAT+NNMAT
               O16(J)=XXO16(J,1)
               C12(J)=XXC12(J,1)
               ANA23(J)=XXNA23(J,1)
               AL27(J)=XXAL27(J,1)
            END DO
         END IF



c!!!!!!! Se voglio carbonio alla chiappinì/MM2002 
c!!!
c!!         IF (ZETAII.LE.1E-5) THEN
c!!            DO J=NMAT+1,NMAT+NNMAT
c!!               C12(J)=XXC12(J,1)
c!!            END DO
c!!         ELSE IF (ZETAII.LE.4E-3) THEN
c!!            DO J=NMAT+1,NMAT+NNMAT
c!!               C12(J)=XXC12(J,2)
c!!            END DO
c!!         ELSE IF (ZETAII.LE.1.2E-2) THEN
c!!            DO J=NMAT+1,NMAT+NNMAT
c!!               C12(J)=XXC12(J,3)
c!!            END DO
c!!         ELSE 
c!!            DO J=NMAT+1,NMAT+NNMAT
c!!               C12(J)=XXC12(J,5)
c!!            END DO
c!!         END IF
c
c!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
c!!!
      END IF

c --------------------------------------
C Definizione delle masse espulse sotto forma dei diversi elementi chi-
C mici considerati per valori della metallicita` iniziale intermedi ri-
C spetto quelli considerati da van den Hoek & Groenewegen (n.b.: si in-
C terpola linearmente sugli intervalli).
c --------------------------------------

 


      IF (ZETAII.LE.0.001) THEN
         DO  J=1,NMAT
            HE(J)=HE1(J)
            C12(J)=C121(J)
            O16(J)=O161(J)
            AN14(J)=AN141(J)
            C13(J)=C131(J)
            RM(J)=RM1(J)
         ENDDO


      ELSE IF (ZETAII.GT.0.001.AND.ZETAII.LE.0.004) THEN
         A=(0.004-ZETAII)/(0.004-0.001)
         B=(ZETAII-0.001)/(0.004-0.001)
         DO  J=1,NMAT
            HE(J)=A*HE1(J)+B*HE2(J)
            C12(J)=A*C121(J)+B*C122(J)
            O16(J)=A*O161(J)+B*O162(J)
            AN14(J)=A*AN141(J)+B*AN142(J)
            C13(J)=A*C131(J)+B*C132(J)
            RM(J)=A*RM1(J)+B*RM2(J)
            
         ENDDO
         

      ELSE IF (ZETAII.GT.0.004.AND.ZETAII.LE.0.008) THEN
         C=(0.008-ZETAII)/(0.008-0.004)
         D=(ZETAII-0.004)/(0.008-0.004)
         DO J=1,NMAT
            HE(J)=C*HE2(J)+D*HE3(J)
            C12(J)=C*C122(J)+D*C123(J)
            O16(J)=C*O162(J)+D*O163(J)
            AN14(J)=C*AN142(J)+D*AN143(J)
            C13(J)=C*C132(J)+D*C133(J)
            RM(J)=C*RM2(J)+D*RM3(J)
         ENDDO
         
      ELSE IF (ZETAII.GT.0.008.AND.ZETAII.LE.0.020) THEN
         E=(0.020-ZETAII)/(0.020-0.008)
         F=(ZETAII-0.008)/(0.020-0.008)
         DO  J=1,NMAT
            HE(J)=E*HE3(J)+F*HE4(J)
            C12(J)=E*C123(J)+F*C124(J)
            O16(J)=E*O163(J)+F*O164(J)
            AN14(J)=E*AN143(J)+F*AN144(J)
            C13(J)=E*C133(J)+F*C134(J)
            RM(J)=E*RM3(J)+F*RM4(J)
         ENDDO
      ELSE IF (ZETAII.GT.0.020.AND.ZETAII.LE.0.040) THEN
         U=(0.040-ZETAII)/(0.040-0.020)
         V=(ZETAII-0.020)/(0.040-0.020)
         DO J=1,NMAT
            HE(J)=U*HE4(J)+V*HE5(J)
            C12(J)=U*C124(J)+V*C125(J)
            O16(J)=U*O164(J)+V*O165(J)
            AN14(J)=U*AN144(J)+V*AN145(J)
            C13(J)=U*C134(J)+V*C135(J)
            RM(J)=U*RM4(J)+V*RM5(J)
         ENDDO
      ELSE IF (ZETAII.GT.0.040) THEN
         DO J=1,NMAT
            HE(J)=HE5(J)
            C12(J)=C125(J)
            O16(J)=O165(J)
            AN14(J)=AN145(J)
            C13(J)=C135(J)
            RM(J)=RM5(J)
         ENDDO
      END IF

      IF (J2.NE.1) GO TO 50
         
            CALL INTERP(H,ZETAII,QMA,QMR,QMHE,QMC,QMO,QMN,QM3,
     $WNS,WC3S,WCUS,WZNS,WKRS,QNE,QMG,QSI,QZO,QCA,
     $QFE,QCU,QZN,QNI,QKR,QLI,
     $WKPS,WSCS,WTIS,WVANAS,WCRS,WMangaS,WCOS,WNaS,WAlS,
     $QKP,QSC,QTI,QVANA,QCR,QManga,QCO,QBAs,QBAr,QSrr,QEU,QLa,
     $QY,QZr,QNa,QAl)
 


C     Nel caso di massa inferiore di 0.5 Msun ==> niente massa espulsa!!!!
C     
      IF (H.GE.M1) GO TO 5
      
c      !!! per stelle di massa minore  di  0.5 masse solari
c      !!! pone i termini a zero non avendo massa espulsa!!!!!!!!!!!!
 
      D=1.
      QC=D
      Q4=D
      WC=0.0

 17   CCWC=0.
      COWC=0.
      CNWC=0.
      C3WC=0.
      CNEWC=0.
      CMGWC=0.
      CSWC=0.
      CZWC=0.0
      CAWC=0.0
      CFWC=0.
      CUWC=0.
      CZNWC=0.
      CNIWC=0.
      CKRWC=0.
C     
C     New  isotopes
C     
      CKPWC=0.
      CSCWC=0.
      CTIWC=0.
      CVANAWC=0.
      CCRWC=0.
      CMangaWC=0.
      CCOWC=0.
      
c Aggiunte GabC
cccccccccccccccc
      CBAsWC=0.
      CBArWC=0.
      CSrrWC=0.
      CEUWC=0.
      CLaWC=0.
      CYWC =0.
      CZrWC=0.
      CNaWC=0.
      CAlWC=0.

ccccccccccccccccc

C
C     end new isotopes
C
      GO TO 200

 5    CONTINUE


c Aggiunte GabC
cccccccccccccccc
      CBAsWC=0.
      CBArWC=0.
      CSrrWC=0.
      CEUWC=0.
      CLaWC=0.
ccccccccccccccccc



      D=QMR
      Q4=QMA
      QC=Q4-QMHE
      WC=QC-D
C____________________________________________________________
C
C     WC=QMC+QMO+QMN+QM3+QNE+QMG+QSI+QZO+QCA+QFE+QCU+QZN+QNI+QKR.
C____________________________________________________________
C

      if((h.gt.8.).and.(h.lt.11.)) then
         
         wc=0.0665

      endif
C
C     forse qualche problema? comunque definisce quasi tutto torna  come
C     zero e per una stella di massa inferiore a 0.5 Msun
C     ????????????????????????????

      IF(WC.LE.0.) THEN
      WC=0.
      GO TO 17
      ENDIF
c     ???????????????????????????

  
C      WRITE(60,2000) WC
C 2000 FORMAT(5X,1E12.5)

      CCWC=QMC/WC
      COWC=QMO/WC
      CNWC=QMN/WC
      C3WC=QM3/WC
      CNEWC=QNE/WC
      CMGWC=QMG/WC
      CSWC=QSI/WC
      CZWC=QZO/WC
      CAWC=QCA/WC
      CFWC=QFE/WC
      CUWC=QCU/WC
      CZNWC=QZN/WC
      CNIWC=QNI/WC
      CKRWC=QKR/WC
C
C     new isotopes
C

      CKPWC=QKP/WC

      CSCWC=QSC/WC

      CTIWC=QTI/WC

      CVANAWC=QVANA/WC

      CCRWC=QCR/WC

      CMangaWC=QManga/WC

      CCOWC=QCO/WC
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CEUWC=QEU/WC
    
      CBAsWC=QBAs

      CBArWC=QBAr/WC
 
      CSrrWC=QSrr/WC

      CLaWC=QLa/WC

      CYWC=QY/WC

      CZrWC=QZr/WC

      CNaWC=QNa/WC

      CAlWC=QAl/WC




c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
 200  CONTINUE


C Calcolo di Q3, valido per tutto il range di masse.


      DD1=D
      ZUT=0.6
      Q3=AMAX1(DD1,ZUT)

C Calcolo di W3 nuovo (Tosi 1995).
c Correzione per l'He3 fatta nel Novembre 1999:
c prima il valore X3(MS) rimaneva fisso (2.1E-4), MA questo non va bene 
c perche' se X3(MS) = X3i +3/2 X2i dove X3i e X2i sono le abbondanze di 
c 3He e D alla formazione della stella, X3(MS) varia durante l'evoluzione 
c del disco!
c Quindi nelle equazioni contenenti X3(MS), X3(MS) deve essere fatto 
c variare. X3(MS) va ricalcolato ad ogni passo temporale:
c X3_MS(I)=X3(I-1)+3/2*X2(I-1) 
c utilizando le abbondanze al momento T(I-1).
c Si usa il risultato di Charbonnel: il 93% delle stelle fra 1-2 Msun
c presenta extra-mixing e non produce 3He.
C      ASP=2.1E-4
      
      IF(H.GT.0.65.AND.H.LE.2.5) THEN
         ATOS=(0.00007+0.00135*H**(-2.2))
         BTOS=(0.55-0.3*ALOG10(H))
         W3=((ATOS+BTOS*(ASP-2.1E-4))/(H*0.77))*0.07
      ELSE IF(H.GT.2.5.AND.H.LE.25) THEN
         ATOS=0.000111+0.00085/(H**2)-0.5E-6*H**1.5
         BTOS=0.485-0.022*EXP(H/10)
         CTOS=1.7E-5*EXP(-(ASP/7.5E-6)**2)
         W3=(ATOS+BTOS*(ASP-2.1E-4)-CTOS)/(H*0.77)
      ELSE IF(H.GT.25) THEN
         ATOS=0.00011-4E-7*H
         BTOS=0.332
         W3=(ATOS+BTOS*(ASP-2.1E-4))/(H*0.77)
      END IF
 127  CONTINUE

      W2=0.


c --------------------------------------
C Calcolo di W7.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Si puo` considerare Mup variabile o Mup costante, commentando
C opportunamente le righe che seguono.
C      EMUP=333.3*(XTIN(6)-1.E-03)+5.
C      IF(EMUP.GT.8.0) EMUP=8.0

      EMUP=8.0


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      IF(H.GT.10.) THEN
C      W7=QLI
C      W7=0.
C      ELSE IF(H.GE.5.0.AND.H.LE.EMUP) THEN
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Vengono date diverse prescrizioni per la produzione di litio
C da stelle con massa compresa tra 5 Msol e Mup, in dipendenza
C dall' abbondanza in numero di litio nel loro inviluppo.
C N(Li) = 5.0 dex:
C      W7=(4.0E-7*H-5.50E-7)/H
C N(Li) = 4.15 dex:
C      W7=(5.6E-8*H-7.7E-8)/H
C N(Li) = 3.5 dex:
C      W7=(1.2E-8*H-1.65E-8)/H
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      ELSE IF(H.GE.2.0.AND.H.LT.5.0) THEN
C      W7=5.E-9/H
C      ELSE
C      W7=0.
C      END IF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Nuove prescrizioni 
C Interpolazione sulla griglia di metallicita`:
C determina semplicemente YLI 




      IF(ZETAII.LE.0.0002) THEN
         DO 707 J=1,10
            YLI(J)=YLI0(J)
 707     CONTINUE


      ELSE IF(ZETAII.GT.0.0002.AND.ZETAII.LE.0.0006) THEN
         DO 708 J=1,10
            YLI(J)=(.0006-ZETAII)/.0004*YLI0(J)
     $           +(ZETAII-.0002)/.0004*YLI1(J)


 708     CONTINUE
      ELSE IF(ZETAII.GT.0.0006.AND.ZETAII.LE.0.001) THEN
         DO 709 J=1,10
            YLI(J)=(.001-ZETAII)/.0004*YLI1(J)
     $           +(ZETAII-.0006)/.0004*YLI2(J)
 709     CONTINUE


      ELSE IF(ZETAII.GT.0.001.AND.ZETAII.LE.0.004) THEN
         DO 710 J=1,10
            YLI(J)=(.004-ZETAII)/.003*YLI2(J)
     $           +(ZETAII-.001)/.003*YLI3(J)
 710     CONTINUE


      ELSE IF(ZETAII.GT.0.004.AND.ZETAII.LE.0.01) THEN
         DO 711 J=1,10
            YLI(J)=(.01-ZETAII)/.006*YLI3(J)
     $           +(ZETAII-.004)/.006*YLI4(J)
 711     CONTINUE


      ELSE IF(ZETAII.GT.0.01.AND.ZETAII.LE.0.02) THEN
         DO 712 J=1,10
            YLI(J)=(.02-ZETAII)/.01*YLI4(J)
     $           +(ZETAII-.01)/.01*YLI5(J)
 712     CONTINUE


      ELSE IF(ZETAII.GT.0.02) THEN
         DO 713 J=1,10
            YLI(J)=YLI5(J)
 713     CONTINUE
      END IF
      
C Interpolazione sulla griglia di massa:


      IF(H.GE.2.2.AND.H.LE.2.5) THEN
         J=1
         K=J
         KK=J+1
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.2.5.AND.H.LE.3.0) THEN
         J=2
         K=J
         KK=J+1
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.3.0.AND.H.LE.3.3) THEN
      J=3
      K=J
      KK=J+1
      W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $     +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.3.3.AND.H.LE.3.5) THEN
         J=4
         K=J
         KK=J+1
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.3.5.AND.H.LE.4.0) THEN
         J=5
         K=J
         KK=J+1
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.4.0.AND.H.LE.4.5) THEN
         J=6
         K=J   
         KK=J+1
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.4.5.AND.H.LE.5.0) THEN
         J=7
         K=J
         KK=J+1
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.5.0.AND.H.LE.5.5) THEN
         J=8
         K=J
         KK=J+1   
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.5.5.AND.H.LE.6.0) THEN
         J=9
         K=J
         KK=J+1
         W7=(PMASS(KK)-H)/(PMASS(KK)-PMASS(K))*YLI(K)
     $        +(H-PMASS(K))/(PMASS(KK)-PMASS(K))*YLI(KK)
      ELSE IF(H.GT.10.) THEN
         W7=QLI
      ELSE IF(H.GE.1.0.AND.H.LE.2.0) THEN
         W7=4.9E-8/H
      ELSE
         W7=0.
      END IF
c     --------------------------------------
      
      W4=1.
      
      
C     Calcolo di Q15.

      IF(H.GT.2.) THEN
         Q15=(1.-0.68+0.01*ALOG10(H))*(1.-Q4)
      ELSE
         Q15=(1.-0.705+0.516*ALOG10(H))*(1.-Q4)
      END IF


C Definizione della matrice QM(I,J).

      DO 160 J=1,NMAX
         DO 140 I=1,NMAX
            QM(I,J)=0.
 140     CONTINUE
 160  CONTINUE


C Quanto segue e` valido per NMAX=21.
c --------------------------------------
c Q(1,1), Q(2,2), Q(3,3), Q(3,2), Q(3,1) foram modificados para
c usar os novos yields da Tosi para He3.
c --------------------------------------


      QM(1,1)=1.-Q4-W7

      QM(2,2)=0.

      QM(3,1)=W3

      QM(3,2)=0.

      QM(3,3)=0.

      QM(4,1)=Q4-QC

      QM(4,2)=Q3-QC

      QM(4,3)=Q3-QC

      QM(4,4)=1.-QC-W7

      QM(4,20)=W4


      DO L=1,4

         QM(20,L)=W7

         QM(5,L)=CCWC*WC
  
         QM(6,L)=COWC*WC
  
         QM(7,L)=CNWC*WC
 
         QM(8,L)=C3WC*WC
         
         QM(10,L)=CNEWC*WC
      
         QM(11,L)=CMGWC*WC

         QM(12,L)=CSWC*WC

         QM(13,L)=CZWC*WC
 
         QM(14,L)=CAWC*WC
  
         QM(15,L)=CFWC*WC
 
         QM(16,L)=CUWC*WC
 
         QM(17,L)=CZNWC*WC
 
         QM(18,L)=CNIWC*WC
      
         QM(19,L)=CKRWC*WC

C
C     new isotopes
C
    
         QM(22,L)=CKPWC*WC
         
         QM(23,L)=CSCWC*WC
         
         QM(24,L)=CTIWC*WC
         
         QM(25,L)=CVANAWC*WC
         
         QM(26,L)=CCRWC*WC
         
         QM(27,L)=CMangaWC*WC
         
         QM(28,L)=CCOWC*WC


C      aggiunte GabC
CCCCCCCCC
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         QM(29,L)=CBAsWC

         QM(30,L)=CBArWC*WC

         QM(31,L)=CSrrWC*WC
          
         QM(32,L)=CEUWC*WC

         QM(33,L)=CLAWC*WC

         QM(34,L)=CYWC*WC

         QM(35,L)=CZrWC*WC

         QM(36,L)=CNaWC*WC

         QM(37,L)=CAlWC*WC
     

      ENDDO
CCCCCCCCC

      QM(7,21)=Q15

      QM(8,21)=Q15

      QM(21,7)=0.

      QM(21,8)=0.

      QM(21,21)=0.


      IF(H.GT.M2) GO TO 1936

c --------------------------------------
C I 6 elementi di matrice seguenti sono stati modificati rispetto quelli
C utilizzati da Chiappini et al. (1997).
c --------------------------------------
      QM(5,5)=1.-Q4-WNAS-WALS
      QM(6,6)=1.-Q4-WNAS-WALS
      QM(7,5)=0.
      QM(7,6)=0.
      QM(8,5)=0.
      QM(8,6)=0.
      GO TO 1937

 1936 CONTINUE
      QM(5,5)=1.-Q4-WNS-WC3S-WNAS-WALS
      QM(6,6)=1.-Q4-WNS-WNAS-WALS
      QM(7,5)=WNS+Q4-QC
      QM(7,6)=WNS+Q4-QC
      QM(8,5)=WC3S
      QM(8,6)=0.

 1937 QM(7,7)=1.-QC-WNAS-WALS
      QM(8,8)=1.-QC

      DO  L=5,8

         QM(9,L)=WC

      ENDDO

 319  CONTINUE

      QM(16,15)=WCUS
      QM(17,15)=WZNS
      QM(19,15)=WKRS
      QM(22,15)=WKPS
      QM(23,15)=WSCS
C     QM(24,15)=WTIS
      QM(25,15)=WVANAS
      QM(26,15)=WCRS
      QM(27,15)=WMangaS
      QM(28,15)=WCOS



C Aggiunta Cescutti
CCCCCCCCCCCCCCCCCCCCCCCCCCC

      QM(31,15)=0.

CCCCCCCCCCCCCCCCCCCCCCCCC


      QM(9,9)=  1.-D
      QM(10,10)=1.-D
      QM(11,11)=1.-D
      QM(12,12)=1.-D
      QM(13,13)=1.-D
      QM(14,14)=1.-D
      QM(15,15)=1.-D
      QM(16,16)=1.-D
      QM(17,17)=1.-D
      QM(18,18)=1.-D
      QM(19,19)=1.-D
      QM(22,22)=1.-D
      QM(23,23)=1.-D
      QM(24,24)=1.-D
      QM(25,25)=1.-D
      QM(26,26)=1.-D
      QM(27,27)=1.-D
      QM(28,28)=1.-D




C     Cescutti   

CCCCCCCCCCCCCCCCCCCCC
      QM(29,29)=1.-D
      QM(30,30)=1.-D
      QM(31,31)=1.-D
      QM(32,32)=1.-D
      QM(33,33)=1.-D
      QM(34,34)=1.-D
      QM(35,35)=1.-D
      QM(36,36)=1.-D
      QM(37,37)=1.-D

CCCCCCCCCCCCCCCCCCCCCc
C     KUNA


      QM(36,5)=WNaS
      QM(36,6)=WNaS
      QM(36,7)=WNaS
      QM(37,5)=WAlS
      QM(37,6)=WAlS
      QM(37,7)=WAlS



      QM(20,20)=1.-W4

  50  CONTINUE

      Q=QM(J1,J2)

      RETURN

      END


      FUNCTION QBIN(H,J1,J2,ZETAII)
      PARAMETER (NMAX=37)  
      include 'lettura.inc'
      include 'elementi.inc'
      REAL M4
      
      DIMENSION QM(NMAX,NMAX),XSTIN(NMAX)
      DIMENSION COSN1(50),COSN2(50) 
      COMMON/CUATRO/COSN1,COSN2
      COMMON/A/ QM
      COMMON/DATA/ LET
c     !COMMON/LUCKY/ IKU
      COMMON/Z/ ZETA
C      COMMON/ZII/ZETAIN,ZETIN
      COMMON/M/NMAT,NNMAT
      COMMON/ELDS/IELD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Blocco COMMON per la nuova formulazione di ASP
      COMMON/HEL/ASP
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      M4=1.5

c!      IF(IKU.EQ.1) GO TO 1002
      IF(J2.NE.1) GO TO 60
c! 1002 CONTINUE


      CALL INTERP(H,ZETAII,QMA,QMR,QMHE,QMC,QMO,QMN,QM3,
     $     WNS,WC3S,WCUS,WZNS,WKRS,QNE,QMG,QSI,QZO,QCA,
     $     QFE,QCU,QZN,QNI,QKR,QLI,
     $     WKPS,WSCS,WTIS,WVANAS,WCRS,WMangaS,WCOS,WNaS,WAlS,
     $     QKP,QSC,QTI,QVANA,QCR,QManga,QCO,QBAs,QBAr,QSrr,QEU,QLa,
     $     QY,QZr,QNa,QAl)





      IF(H.GE.M4) GO TO 1
      D=1.
      QC=D
      Q4=D
17    CCWC=0.
      COWC=0.
      CNWC=0.
      C3WC=0.
      CNEWC=0.
      CMGWC=0.
      CSWC=0.
      CZWC=0.0
      CAWC=0.0
      CFWC=0.
      CUWC=0.
      CZNWC=0.
      CNIWC=0.
      CKRWC=0.

C     New  isotopes
C       
      CKPWC=0.
      CSCWC=0.
      CTIWC=0.
      CVANAWC=0.
      CCRWC=0.
      CMangaWC=0.
      CCOWC=0.

c!! aggiunte GabC
C      Niente produzione in SNIa
C      -----------------------------
      CEUWC=0.
      CBAsWC=0.
      CSrrWC=0.
      CBArWC=0.
      CLaWC =0.
      CYWC=0.
      CZrWC =0.
      CNaWC=0.
      CAlWC =0.

      
C
C     end new isotopes


      GO TO 100
  1   CONTINUE
      D=0.
      Q4=QMA
      QC=QMA-QMHE
C      WC=QMO+QMC+QMN+QM3+(1.4/(0.99*H))


      WC=QC-D

c!! aggiunte GabC
C      Niente produzione in SNIa
C      -----------------------------
      CEUWC=0.
      CBAsWC=0.
      CSrrWC=0.
      CBArWC=0.
      CLaWC=0.
      CYWC=0.
      CZrWC =0.
      CNaWC=0.
      CAlWC =0.
      
C



      if(wc.le.0.) then
         wc=0.
      go to 17
      endif
C------------------------------------------------------
C                                                     |
C                 prescriptions for SNIa              |
C                                                     |
C------------------------------------------------------


      DO JK=1,NMAX
      cosn1(JK)=1.                         
      END DO
C
C     Coefficients
C
       cosn1(6)=1.                       
       cosn1(11)=7.      
       cosn1(12)=1.      
       cosn1(14)=1.      
       cosn1(17)=22.5    
       cosn1(18)=0.15    
       cosn1(22)=1/10.   
       cosn1(23)=58.     
       cosn1(24)=2.0     
       cosn1(26)=1/1.5   
       cosn1(27)=0.7     
       cosn1(28)=15.     
       cosn1(29)=0.      
       cosn1(30)=0.      
       cosn1(31)=0.      
       cosn1(32)=0.      
       cosn1(33)=0.      

c!!!!!!!
       cosn1(34)=0.          
       cosn1(35)=0.         
       cosn1(36)=0.        
       cosn1(37)=0.       
       


       
C
C
C         
      CCWC=0.048/(0.99*WC*H)+QMC/WC            
      COWC=0.143/(0.99*WC*H)+QMO/WC* cosn1(6) 
      CNWC=(QMN+1.16e-6)/WC                  
      C3WC=(QM3+1.40e-6)/WC                  
      CNEWC=0.00202/(0.99*WC*H)             
      CMGWC=(0.0085/(0.99*WC*H))*cosn1(11) 
      CSWC=0.154/(0.99*WC*H)*cosn1(12)    
      CZWC=0.0846/(0.99*WC*H)            
      CAWC=0.0119/(0.99*WC*H)*cosn1(14) 
      CFWC=0.6/(0.99*WC*H)             
      CUWC=2.00E-4/(0.99*WC*H)        
C      CZNWC=6.30E-4/(0.99*WC*H)                 
      CZNWC=2.82E-5/(0.99*WC*H)*cosn1(17)

C      CNIWC=0.0695/(0.99*WC*H))                
      
C       CNIWC=0.0695/(4*(0.99*WC*H))           
        CNIWC=0.122/(0.99*WC*H)*cosn1(18)     
      CKRWC=0.
C      CZNWC=0.
C      CUWC=0.
C     New  isotopes            
       CKPWC=1.19E-2/(0.99*WC*H)*cosn1(22)
       CSCWC=2.21e-7/(0.99*WC*H)*cosn1(23)
       CTIWC=2.08e-4/(0.99*WC*H)*cosn1(24)
       CVANAWC=7.49E-5/(0.99*WC*H)
       CCRWC=8.5E-3/(0.99*WC*H)*COSN1(26)

       CCOWC=1.04E-4/(0.99*WC*H)*cosn1(28)

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!VARIAZIONI NEGLI YIELDS DI SNIa per il Manganese!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       CmangaWC=(8.87E-3/(0.99*WC*H)) 
       CmangaWC=((83*zetaii)**.65)*(8.87E-3/(0.99*WC*H)) 
       CmangaWC=((83*zetaii)**.15)*(8.87E-3/(0.99*WC*H)) 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



C     end new isotopes



 100   CONTINUE
C     Calcolo di Q3, valido per tutto il range di masse.
       DD1=D
       ZUT=0.6
       Q3=AMAX1(DD1,ZUT)
C     Calcolo di W3 nuovo (Tosi 1995).
C     ASP=2.1E-4
       IF(H.GT.0.65.AND.H.LE.2.5) THEN
          ATOS=0.00007+0.00135*H**(-2.2)
          BTOS=0.55-0.3*ALOG10(H)
          W3=(ATOS+BTOS*(ASP-2.1E-4))/(H*0.77)
       ELSE IF(H.GT.2.5.AND.H.LE.25) THEN
          ATOS=0.000111+0.00085/(H**2)-0.5E-6*H**1.5
          BTOS=0.485-0.022*EXP(H/10)
          CTOS=1.7E-5*EXP(-(ASP/7.5E-6)**2)
          W3=(ATOS+BTOS*(ASP-2.1E-4)-CTOS)/(H*0.77)
       ELSE IF(H.GT.25) THEN
          ATOS=0.00011-4E-7*H
          BTOS=0.332
          W3=(ATOS+BTOS*(ASP-2.1E-4))/(H*0.77)
       END IF
C Calcolo di W3.
C      A=0.                                                             
C      W3=4.7E-4*A/(H**2)                                               
C Calcolo di W2. (Tosi 1995: W2 e` posto uguale a zero.)
C      IF(H.GT.2.) GO TO 120
      W2=0.
C      GO TO 130
C 120  IF(H.GT.4.) GO TO 125
C      W2=(1.-Q3)*(H-2.)/2.
C      GO TO 130
C 125  W2=1.-Q3
C 130  CONTINUE
C Definizione della matrice QM(I,J).
      DO 160 J=1,NMAX
      DO 140 I=1,NMAX
      QM(I,J)=0.
  140 CONTINUE
  160 CONTINUE
C Quanto segue e` valido per NMAX=21.
c --------------------------------------
c Mudei QM(1,1), QM(2,2), QM(3,1), QM(3,2) e QM(3,3)
c para usar novos yields do He3 de Tosi.
c -------------------------------------- 
      QM(1,1)=1.-Q4
      
      QM(2,2)=0.
      
      QM(3,1)=W3
      
      QM(3,2)=0.
      
      QM(3,3)=0.
      
       QM(4,1)=Q4-QC
      
      QM(4,2)=Q3-QC
      
      QM(4,3)=Q3-QC
      
      QM(4,4)=1.-QC
      
      DO  L=1,4
  
         QM(5,L)=CCWC*WC
        
         QM(6,L)=COWC*WC
        
         QM(7,L)=CNWC*WC
        
         QM(8,L)=C3WC*WC
        
         QM(10,L)=CNEWC*WC
        
         QM(11,L)=CMGWC*WC
        
         QM(12,L)=CSWC*WC
        
         QM(13,L)=CZWC*WC
        
         QM(14,L)=CAWC*WC
        
         QM(15,L)=CFWC*WC
        
         QM(16,L)=CUWC*WC
        
         QM(17,L)=CZNWC*WC
        
         QM(18,L)=CNIWC*WC
        
         QM(19,L)=CKRWC*WC

C
C     new isotopes
C
         QM(22,L)=CKPWC*WC
  
         QM(23,L)=CSCWC*WC
     
         QM(24,L)=CTIWC*WC
     
         QM(25,L)=CVANAWC*WC
    
         QM(26,L)=CCRWC*WC
   
         QM(27,L)=CMangaWC*WC
     
         QM(28,L)=CCOWC*WC
 
c!! aggiunte GabC
 
         QM(29,L)=CBAsWC*WC
         
         QM(30,L)=CBArWC*WC
         
         QM(31,L)=CSrrWC*WC
      
         QM(32,L)=CEUWC*WC
      
         QM(33,L)=CLaWC*WC
  
         QM(34,L)=CYWC*WC
  
         QM(35,L)=CZrWC*WC
  
         QM(36,L)=CNaWC*WC
  
         QM(37,L)=CAlWC*WC
  

      ENDDO


c --------------------------------------
C I 6 elementi di matrice seguenti sono stati modificati rispetto quelli
C utilizzati da Chiappini et al.
c --------------------------------------

      QM(5,5)=1.-Q4
      
      QM(6,6)=1.-Q4
      
      QM(7,5)=0.
      
      QM(7,6)=0.
      
      QM(8,5)=0.
      
      QM(8,6)=0.
      
      QM(7,7)=1.-QC
      
      QM(8,8)=1.-QC

 
c!! MODIFICHE E AGGIUNTE  GabC

      
      DO  L=5,8
         
         QM(9,L)=WC
      
      ENDDO

      DO L=9,19

         QM(L,L)=1.-D
      
      ENDDO

      DO L=22,37

         QM(L,L)=1.-D

      ENDDO




 
c!! MODIFICHE E AGGIUNTE  GabC


  
 60   CONTINUE
      QBIN=QM(J1,J2)

 7000 FORMAT(1X,6E12.5 /)

      RETURN
      END






      FUNCTION TAU(H)
      IF(H.LE.1.3) THEN
      TAU=10**(-0.6545*ALOG10(H)+1.0)
      ELSE IF(H.LE.3.) THEN
      TAU=10**(-3.7*ALOG10(H)+1.35)
      ELSE IF(H.LE.7.) THEN
      TAU=10**(-2.51*ALOG10(H)+0.77)
      ELSE IF(H.LE.15) THEN
      TAU=10**(-1.78*ALOG10(H)+0.17)
      ELSE IF(H.LE.60) THEN
      TAU=10**(-0.86*ALOG10(H)-0.94)
      ELSE
      TAU=1.2*H**(-1.85)+3.E-3
      END IF
      RETURN
      END






      FUNCTION ASSA(T)
      IF(T.LE.3.36E-3) THEN
      ASSA=8.00E1
      ELSE IF(T.LE.3.39E-3) THEN
      C=(ALOG10((T-3.E-3)/1.2))/(-1.85)
      ASSA=10**C
      ELSE IF(T.LE.0.0119) THEN
      ASSA=10**((ALOG10(T)+0.94)/(-0.86))
      ELSE IF(T.LE.0.0445) THEN
      ASSA=10**((ALOG10(T)-0.17)/(-1.78))
      ELSE IF(T.LE.0.37) THEN
      ASSA=10**((ALOG10(T)-0.77)/(-2.51))
      ELSE IF(T.LE.8.48) THEN
      ASSA=10**((ALOG10(T)-1.35)/(-3.7))
      ELSE
      ASSA=10**((ALOG10(T)-1.0)/(-0.6545))
      END IF
      RETURN
      END





      FUNCTION PSI(H)
      PARAMETER (NMAX=37)  
      REAL MI,MS,M1,M2,M3,M4,M6
      DIMENSION IND(NMAX),GO(3000),G(3000,NMAX),T(3000)
      DIMENSION X(3000,NMAX)
      DIMENSION XCAL(NMAX),XTIN(NMAX),XSTAR(NMAX),XSTIN(NMAX)
C      DIMENSION COSN1(NMAX),COSN2(NMAX)
C      COMMON/CUATRO/COSN1,COSN2           

      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     $GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     $KMAX,I,ITER,IND,KDEMAX
      COMMON/Z/ZETA
      COMMON/SUP/AO1,BO1,CO1
C Calcolo della funzione iniziale di massa (Scalo Zita).
C      IF(T(I).LE.0.0036) GO TO 1953
      UM1=1.35
C      IF (T(I).GT.2.) THEN
      UM2=1.7
C      ELSE
C      UM2=1.35
C      END IF
      M3=1.0
      M6=2.0
      MS=80.0
      MI=0.10
      ZITA=0.30
      A=MS**(1.-UM2)
      B=M6**(1.-UM2)
      C=M6**(1.-UM1)
      D=M3**(1.-UM1)
      CCO1=1.-UM1
      CO2=1.-UM2
      BO1=ZITA/(M6**(UM1-UM2)*(C-D)/CCO1+(A-B)/CO2)
C      WRITE(60,1977) AO1,BO1
C 1977 FORMAT(1X,'AO1=',1E12.5,1X,'BO1=',1E12.5)
      AO1=BO1*M6**(UM1-UM2)
      IF(H.LT.MI) GO TO 10
      IF(H.GE.MI.AND.H.LT.M6) GO TO 9
      IF(H.GE.M6.AND.H.LE.MS) GO TO 8
      IF(H.GT.MS) GO TO 10
 8    PSI=BO1/(H**UM2)
      GO TO 11
 9    PSI=AO1/(H**UM1)
      GO TO 11
C      IF(ZETA.LT.2.E-3) GO TO 20
C      UM1=ALOG10(ZETA)+4.05
C      GO TO 30
C 20   UM1=1.05
C 30   CONTINUE
C      UM1=1.35
C      UM1=0.5*ALOG10(ZETA)+3.00
C      MI=-0.9*ALOG10(ZETA)-1.7
C 1953 CONTINUE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      UM1=1.1
C      MI=0.1
C      MS=80.0
C      M3=10.0
C      CO=1.-UM1
C      ZITA=1.0
C      A=MS**(1.-UM1)
C      B=MI**(1.-UM1)
C      CO1=CO*ZITA/(A-B)
C      IF(H.LT.MI) GO TO 10
C      IF(H.GT.MS) GO TO 10
C      PSI=CO1/(H**UM1)
C      GO TO 11
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Funzione iniziale di massa a tre pendenze.
C      UM1=0.25
C      UM2=1.30
C      UM3=1.7
C      M1=0.1
C      M2=1.0
C      M3=10.0
C      M4=80.0
C      ZITA=0.30
C      C1=UM1-UM2
C      C2=UM3-UM2
C      C3=1.-UM1
C      C4=1.-UM2
C      C5=1.-UM3
C      A=M1**C3
C      B=M2**C3
C      C=M2**C4
C      D=M3**C4
C      E=M3**C5
C      F=M4**C5
C      BO1=1./(((M2**C1)/C3)*(B-A)+(1./C4)*(D-C)+((M3**C2)/C5)*(F-E))
C      AO1=BO1*M2**C1
C      CO1=BO1*M3**C2
C      IF(H.LT.M1) GO TO 10
C      IF(H.GE.M1.AND.H.LE.M2) GO TO 9
C      IF(H.GT.M2.AND.H.LE.M3) GO TO 12
C      IF(H.GT.M3.AND.H.LE.M4) GO TO 13
C      IF(H.GT.M4) GO TO 10
C 9    PSI=AO1/(H**UM1)
C      GO TO 11
C 12   PSI=BO1/(H**UM2)
C      GO TO 11
C 13   PSI=CO1/(H**UM3)
C      GO TO 11
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 10   CONTINUE
      PSI=0.0
 11   CONTINUE
C      WRITE(60,2897) AO1,BO1,CO1
 2897 FORMAT(1X,'A=',1E12.5,1X,'B=',1E12.5,1X,'C=',1E12.5)
      RETURN
      END

      FUNCTION PNU(H)
      PARAMETER (NMAX=37) 
      REAL MI,MS,M1,M2,M3,M4,M6
      DIMENSION IND(NMAX),XCAL(NMAX),XTIN(NMAX),XSTAR(NMAX),XSTIN(NMAX)
      DIMENSION GO(3000),G(3000,NMAX),X(3000,NMAX),T(3000)
C      DIMENSION COSN1(NMAX),COSN2(NMAX)
C      COMMON/CUATRO/COSN1,COSN2           

      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     $GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     $KMAX,I,ITER,IND,KDEMAX
      COMMON/Z/ZETA
      COMMON/SUP/AO1,BO1,CO1
C SCALO IMF CON ZITA
C      IF(T(I).LE.0.0036) GO TO 1953                                   
      UM1=2.35
C      IF (T(I).GT.2.) THEN
      UM2=2.7
C      ELSE
C      UM2=2.35
C      END IF
      M3=1.0
      M6=2.0
      MS=80.0
      MI=0.10
C      ZITA=0.30
C      A=MS**(1.-UM2)
C      B=M6**(1.-UM2)
C      C=M6**(1.-UM1)
C      D=M3**(1.-UM1)
C      CCO1=1.-UM1
C      CO2=1.-UM2
C      BO1=ZITA/(M6**(UM1-UM2)*(C-D)/CCO1+(A-B)/CO2)
C      AO1=BO1*M6**(UM1-UM2)
      IF(H.LT.MI) GO TO 10
      IF(H.GE.MI.AND.H.LT.M6) GO TO 9
      IF(H.GE.M6.AND.H.LE.MS) GO TO 8
      IF(H.GT.MS) GO TO 10
 8    PNU=BO1/(H**UM2)
      GO TO 11
 9    PNU=AO1/(H**UM1)
      GO TO 11
C CALCOLA IL NUMERO DI STELLE
C 1953 CONTINUE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      UM1=2.1
C      MI=0.10
C      MS=80.0
C      M3=1.0
C      IF(H.LT.MI) GO TO 10
C      IF(H.GT.MS) GO TO 10
C      PNU=CO1/(H**UM1)
C      GO TO 11
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Funzione iniziale di massa Scalo (1986) normalizzata ad 1.
C      UM1=1.25
C      UM2=2.35
C      UM3=2.7
C      M1=0.05
C      M2=1.0
C      M3=2.0
C      M4=200.0
C      IF(H.LT.M1) GO TO 10
C      IF(H.GE.M1.AND.H.LE.M2) GO TO 9
C      IF(H.GT.M2.AND.H.LE.M3) GO TO 12
C      IF(H.GT.M3.AND.H.LE.M4) GO TO 13
C      IF(H.GT.M4) GO TO 10
C 9    PNU=AO1/(H**UM1)
C      GO TO 11
C 12   PNU=BO1/(H**UM2)
C      GO TO 11
C 13   PNU=CO1/(H**UM3)
C      GO TO 11
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 10   CONTINUE
      PNU=0.0
 11   CONTINUE
      RETURN
      END

      SUBROUTINE COSMIC
      PARAMETER (NMAX=37) 
      DIMENSION COSM(NMAX)
      COMMON/Z/ZETA
      COMMON/RAYS/COSM
      ZSUN=.20E-1
      ZSUN2=.20E-3
      ZSUN4=.20E-5
      COEFF=15.45264
      DO 72 J=1,NMAX
 72   COSM(J)=0.
      IF(ZETA.LT.ZSUN4) THEN
      GCR=0.
      ELSE IF(ZETA.GE.ZSUN4.AND.ZETA.LT.ZSUN2) THEN
      GCR=(ZSUN2-ZETA)/(ZSUN2-ZSUN4)*COEFF*0.1E-12+
     $(ZETA-ZSUN4)/(ZSUN2-ZSUN4)*COEFF*0.2E-12
      ELSE IF(ZETA.GE.ZSUN2.AND.ZETA.LT.ZSUN) THEN
      GCR=(ZSUN-ZETA)/(ZSUN-ZSUN2)*COEFF*0.2E-12+
     $(ZETA-ZSUN2)/(ZSUN-ZSUN2)*COEFF*0.3E-11
      ELSE
      GCR=COEFF*0.3E-11
      END IF
      COSM(20)=GCR
      RETURN
      END


      SUBROUTINE INTERP(H,ZETAIIn,QMA,QMR,QMHE,QMC,QMO,QMN,QM3,
     $WNS,WC3S,WCUS,WZNS,WKRS,QNE,QMG,QSI,QZO,QCA,
     $QFE,QCU,QZN,QNI,QKR,QLI,
     $WKPS,WSCS,WTIS,WVANAS,WCRS,WMangaS,WCOS,WNaS,WAlS, 
     $QKP,QSC,QTI,QVANA,QCR,QManga,QCO,QBAs,QBAr,QSrr,QEU2,QLa,
     $QY,QZr,QNa,QAl )
      
      include 'lettura.inc'
      include 'elementi.inc'
      include 'interp.inc'
     
      INTEGER NMAX
      PARAMETER (NMAX=37)  
    
      REAL M1,M2,MAN,M4
      DIMENSION AMLI(50),ALI0(50),ALI1(50),ALI2(50),ALI3(50),ALI4(50)
      DIMENSION QM(NMAX,NMAX),X(3000,NMAX),XCAL(NMAX),XSTAR(NMAX)
      DIMENSION XTIN(NMAX)
      DIMENSION COSN1(50),COSN2(50)
      DIMENSION XSTIN(NMAX), GO(3000), G(3000,NMAX), T(3000)
      COMMON/CUATRO/COSN1,COSN2
      COMMON/WOO/AMLI,ALI0,ALI1,ALI2,ALI3,ALI4
      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/Z/ZETA
C      COMMON/ZII/ZETAIN,ZETIN
      COMMON/DATA/LET
      COMMON/M/NMAT,NNMAT
      COMMON/A/QM
      COMMON/CLAUDIA/AFH
      real bar(2),mbar(2)

      
      
      M1=0.5
      M2=8.0
      M4=8.0
      X1=0.7
      X4=0.28
      X12=0.00361
      X16=0.0108
       X14=0.00111

      CL12=X12/10.
      CL14=X14/10.
      CL16=X16/10.
     
     

 

C     multiplying factor for SNII yields
C     che per prima cosa  azzero!

c------------------------------------------
      DO IK=1,NMAX      
      COSN2(IK)=1.
      end do
C-----------------------------------------------

       COSN2(6)=1.     
       COSN2(11)=7.   
       cosn2(12)=1.  
       cosn2(14)=1. 
       cosn2(18)=1.    
       COSN2(22)=2.   
       COSN2(23)=1.15
       COSN2(24)=15.   
       COSN2(26)=3.   
       cosn2(28)=0.3 
       cosn2(29)=0.    
       cosn2(30)=0.   
       cosn2(31)=0.  
       cosn2(32)=0. 
       cosn2(33)=0.    
       cosn2(34)=0.   
       cosn2(35)=0.  
       cosn2(36)=1. 
       cosn2(37)=1.
 



C aggiunte di GabC

    

C canale primario del processo s main del Bario

       QBAs=0.
       Qla=0.
C canale primario del processo r  dello Stronzio


       QSrr=0.

C canale primario del processo r del Bario

       QBAr=0.

C canale primario del processo r dell'Europio

       QEU2=0.

       QY=0.
       QZr=0.
       QNa=0.
       QAl=0.
     








       


c-------------------------------------------------------------------------------------
c OBS: com os novos yields para SNII de WW95
c ja esta sendo levado em conta a parte que cai de novo
c no core e entao nao preciso mais corrigir por BETA
      BETA=1
C------------------------------------------------------------------------------------ 

    
C  Sceglie il corretto intervallo di massa in cui fare i conti che risulta compredo
C  fra il  KappaKappaesimo e il Kappesimo e
      J=1

      
 3    IF(AMM(J)-H) 10,20,20

 10   J=J+1

      GO TO 3

 20   K=J

      KK=J-1


                             
C---------------------------------------------------------------------------
C Interpolazione sulla griglia di masse per massive stars.                  
C---------------------------------------------------------------------------
C***************************************************************************
C---------------------------------------------------------------------------
      IF(H.gt.M4) then          


         DMM=AMM(K)-AMM(KK)     
         DM=AMM(K)-H           
         DMA=AM(K)-AM(KK)     
         MAN=AM(K)-DMA*DM/DMM
         DMR=RM(K)-RM(KK)   
         DMHE=HE(K)-HE(KK) 
         DMC=C12(K)-C12(KK)     
         DMO=O16(K)-O16(KK)    
         DMN=AN14(K)-AN14(KK) 
         DMC3=C13(K)-C13(KK) 
         DMNS=S14(K)-S14(KK)
         DC3S=C13S(K)-C13S(KK)
         DNE=ANE(K)-ANE(KK)
         DMG=AMG(K)-AMG(KK)
         DZO=BETA*S32(K)-BETA*S32(KK) 
         DSI=BETA*ASI(K)-BETA*ASI(KK)
         DCA=BETA*CA(K)-BETA*CA(KK)
         DFE=AFE(K)-AFE(KK)
         DCU=CU(K)-CU(KK)
         DZN=ZN(K)-ZN(KK)
         DNI=ANI(K)-ANI(KK)
         DKR=AKR(K)-AKR(KK)
C     new isotopes
         DKP=KP(K)-KP(KK)
         DSC=SC(K)-SC(KK)
         DTI=TI(K)-TI(KK)
         DVANA=VANA(K)-VANA(KK)
         DCR=CR(K)-CR(KK)
         DManga=Manga(K)-Manga(KK)
         DCO=CO(K)-CO(KK)
         DMAN=AM(K)-MAN
C     KUNA
         DNA23=ANA23(K)-ANA23(KK)
         DAL27=AL27(K)-AL27(KK)
c     SODIO E ALLUMINIO SECONDARI
         DNAS=ANAS(K)-ANAS(KK)
         DALS=ALS(K)-ALS(KK)
      
         
C    end new isotopes


C Parte per calcolare V del litio!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

      JL=1
 4    IF(AMLI(JL)-H) 11,21,21
 11   JL=JL+1
      GO TO 4
 21   KL=JL
      KKK=JL-1
      DMLI=AMLI(KL)-AMLI(KKK)
      DLI=AMLI(KL)-H
      IF(ZETAIIN.LE.1.7E-9) THEN
      DLITIO=ALI0(KL)-ALI0(KKK)
      VLI=ALI0(KL)-DLITIO*DLI/DMLI
      ELSE IF(ZETAIIN.GT.1.7E-9.AND.ZETAIIN.LE.1.9E-6) THEN
      DLITIO=ALI1(KL)-ALI1(KKK)
      VLI=ALI1(KL)-DLITIO*DLI/DMLI
      ELSE IF(ZETAIIN.GT.1.9E-6.AND.ZETAIIN.LE.1.9E-4) THEN
      DLITIO=ALI2(KL)-ALI2(KKK)
      VLI=ALI2(KL)-DLITIO*DLI/DMLI
      ELSE IF(ZETAIIN.GT.1.9E-4.AND.ZETAIIN.LE.1.9E-3) THEN
      DLITIO=ALI3(KL)-ALI3(KKK)
      VLI=ALI3(KL)-DLITIO*DLI/DMLI
      ELSE
      DLITIO=ALI4(KL)-ALI4(KKK)
      VLI=ALI4(KL)-DLITIO*DLI/DMLI
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      VMR=RM(K)-DMR*DMAN/DMA
      VMHE= HE(K)-DMHE*DMAN/DMA
      VMC=C12(K)-DMC*DMAN/DMA
      VMO=(O16(K)-DMO*DMAN/DMA)*COSN2(6)

C instruction for primary Nitrogen

      VMN=AN14(K)-DMN*DMAN/DMA
      VMC3=C13(K)-DMC3*DMAN/DMA
      VMNS=S14(K)-DMNS*DMAN/DMA
      VC3S=C13S(K)-DC3S*DMAN/DMA
      VNE=ANE(K)-DNE*DMAN/DMA
      VMG=(AMG(K)-DMG*DMAN/DMA)*COSN2(11)
      VZO=BETA*S32(K)-DZO*DMAN/DMA
      VSI=(BETA*ASI(K)-DSI*DMAN/DMA)*COSN2(12)
      VCA=(BETA*CA(K)-DCA*DMAN/DMA)*cosn2(14)
      VFE=AFE(K)-DFE*DMAN/DMA
      VNI=(ANI(K)-DNI*DMAN/DMA)*cosn2(18) 
      VZN=ZN(K)-DZN*DMAN/DMA             
C     new isotopes
      VKP=(KP(K)-DKP*DMAN/DMA)*COSN2(22)
      VSC=(SC(K)-DSC*DMAN/DMA)*COSN2(23)
      VTI=(TI(K)-DTI*DMAN/DMA)*COSN2(24)
      VVANA=VANA(K)-DVANA*DMAN/DMA 
      VCR=(CR(K)-DCR*DMAN/DMA)*COSN2(26)
      VMANga=(MaNga(K)-DMANGA*DMAN/DMA)
      VCO=(CO(K)-DCO*DMAN/DMA)*cosn2(28) 
c     KUNA
      VNA23=(ANA23(K)-DNA23*DMAN/DMA)
      VAL27=(AL27(K)-DAL27*DMAN/DMA)
c     SODIO E ALLUMINIO SECONDARI
      VNAS=(ANAS(K)-DNAS*DMAN/DMA)
      VALS=(ALS(K)-DALS*DMAN/DMA)
     
c     end isotopes       

C Componente primaria di rame e zinco (Thielemann).
       VCU=7.0E-6
       VKR=0.

       QMR=VMR/H
       QMA=MAN/H
       QMHE=VMHE/H
       QMC=VMC/H
       QMO=VMO/H
       QMN=VMN/H
       QM3=VMC3/H
       QNE=VNE/H
       QMG=VMG/H
       QSI=VSI/H
       QZO=VZO/H
       QCA=VCA/H
       QFE=VFE/H
       QCU=VCU/H
       QZN=VZN/H
       QNI=VNI/H
       QKR=VKR/H
       QLI=VLI/H
C     new isotopes
       QKP=VKP/H
       QSC=VSC/H
       QTI=VTI/H
       QVANA=VVANA/H
       QCR=VCR/H
       QManga=VManga/H
       QCO=VCO/H
C     KUNA
       QNA=VNA23/H
       QAL=VAL27/H    
C     end new isotopes

  


      
       WNS=VMNS/H/(X12+X16)     
       WC3S=VC3S/H/X12         

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     QUI CI VA NA/AL SECONDARIO     C   (MASSIVE STARS, s-process)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       WNAS=VNAS/H/(CL12+CL14+CL16)
       WALS=VALS/H/(CL12+CL14+CL16)


       WCUS=0.
       WZNS=0.
       WKRS=0.
       WKPS=0.
       WSCS=0.
       WVANAS =0.
       WCRS=0.
       WMaNgaS=0.
       WCOS=0.

c!!!
c!!!
c!!!       if ((h.ge.12).and.(h.le.40)) then
c!!!          call value4(zetaiin,H,qmanga)
c!!!          call value3(zetaiin,H,qfe)
c!!!          call value2(zetaiin,H,qzn)
c!!!       else
c!!!          if (h.lt.12) then
c!!!             call value4(zetaiin,12.,qmanga)
c!!!             qmanga=qmanga/4*(H-8)
c!!!             call value3(zetaiin,12.,qfe)
c!!!             qfe=qfe/4*(H-8)
c!!!             call value2(zetaiin,12.,qzn)
c!!!             qzn=qzn/4*(H-8)             
c!!!          else
c!!!             call value4(zetaiin,40.,qmanga)
c!!!             call value3(zetaiin,40.,qfe)
c!!!             call value2(zetaiin,40.,qzn)
c!!!          endif
c!!!       endif
c!!!       
c!!!
c!!!
c
c                                                                                      
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!                                                                          !!!!!!
c!!!!!!                   BARIO & EUROPIO                                        !!!!!!
c!!!!!!                   '''''''''''''''                                        !!!!!!
c!!!!!!                                                                          !!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   

       IF(H.GE.12.and.H.le.15.) THEN
c          !BARIUM
          bar(1)=9.e-7
          bar(2)=3.e-8
          mbar(1)=12
          mbar(2)=15.
          call polint(mbar,bar,2,H,qbar,dy)
c          !EUROPIUM
          bar(1)=9.e-7*0.05
          bar(2)=3.e-8*0.1
          mbar(1)=12
          mbar(2)=15.
          call polint(mbar,bar,2,H,qeu2,dy)
c          !LANTHANUM
          bar(1)=9.e-7*0.08
          bar(2)=3.e-8*0.08
          mbar(1)=12
          mbar(2)=15.
          call polint(mbar,bar,2,H,qla,dy)
c          !STRONTIUM
          bar(1)=9.e-7*1.8
          bar(2)=3.e-8*1.8
          mbar(1)=12
          mbar(2)=15.
c          !YTTRIUM
          call polint(mbar,bar,2,H,qsrr,dy)
          bar(1)=9.e-7*.4
          bar(2)=3.e-8*.4
          mbar(1)=12
          mbar(2)=15.
          call polint(mbar,bar,2,H,qy,dy)
c          !ZIRCONIUM
          bar(1)=9.e-7*2.
          bar(2)=3.e-8*6.
          mbar(1)=12
          mbar(2)=15.
          call polint(mbar,bar,2,H,qzr,dy)
          

       endif
c      !BARIUM
      if (H.Gt.15.and.H.lt.30) then
         
          bar(1)=3.e-8
          bar(2)=1.e-9
          mbar(1)=15.
          mbar(2)=30.
          call polint(mbar,bar,2,H,qbar,dy)

c      !EUROPIUM
          bar(1)=3.e-8*0.1
          bar(2)=1.e-9*0.5
          mbar(1)=15.
          mbar(2)=30.
          call polint(mbar,bar,2,H,qeu2,dy)
c      !LANTHANUM  
          bar(1)=3.e-8*0.3
          bar(2)=1.e-9*0.1
          mbar(1)=15.
          mbar(2)=30.
          call polint(mbar,bar,2,H,qla,dy)
c      !STRONTIUM
          bar(1)=3.e-8*1.3*2.5
          bar(2)=1.e-9*1.3*2.5
          mbar(1)=15.
          mbar(2)=30.
          call polint(mbar,bar,2,H,qsrr,dy)
c      !YTTRIUM    

          bar(1)=3.e-8*1.
          bar(2)=1.e-9*1.
          mbar(1)=15.
          mbar(2)=30.
          call polint(mbar,bar,2,H,qy,dy)
c!
c      !ZIRCONIUM   

          bar(1)=3.e-8*5
          bar(2)=1.e-9*5
          mbar(1)=15.
          mbar(2)=30.
          call polint(mbar,bar,2,H,qzr,dy)

       endif 




c$$$c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
c$$$c$$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$c$$$C     FINE BARIO EUROPIO MASSIVE STARS
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC









c$$$     

C Componente secondaria (s-process) di rame e zinco.
C in funzione della massa della stella
C --------------------------------------
 

      
C variazioni per masse maggiori di 15 e minori di 20
      

      IF(H.GE.15.AND.H.LE.20.) THEN

         wcus=7.57e-7/h/1.15e-3 
         wzns=5.78e-7/h/1.15e-3 
         WKRS=1.64e-9/h/1.15e-3
        
      endif

      IF(H.GT.20.AND.H.LE.25.) THEN
         wcus=(2.5049e-5*h-5.0022e-4)/h/1.15e-3
         wzns=(1.0428e-4*h-2.0851e-3)/h/1.15e-3
         wkrs=(5.8167e-7*h-1.1632e-5)/h/1.15e-3
      endif
       
      IF(H.GT.25.) THEN
         wcus=(3.64e-5*h-7.84e-4)/h/1.15e-3
         wzns=(2.97e-5*h-2.225e-4)/h/1.15e-3
         wkrs=(1.78e-7*h-1.54e-6)/h/1.15e-3
      endif
         
C strano ma funziona cosi' in best wzns per stelle massive e' definito zero!

      wzns=0





      endif

    
c-----------------------------------------------------------------------------------
C Interpolazione sulla griglia di masse per low and intermediate mass stars.
c ----------------------------------------------------------------------------------
c***********************************************************************************
c ----------------------------------------------------------------------------------      
      if (H.le.m4) then


         DMM=AMM(K)-AMM(KK)
         DM=AMM(K)-H
         DMR=RM(K)-RM(KK)
         DMHE=HE(K)-HE(KK)
         DMC=C12(K)-C12(KK)
         DMO=O16(K)-O16(KK)
         DMN=AN14(K)-AN14(KK)
         DMC3=C13(K)-C13(KK)
         DMNS=S14(K)-S14(KK)
         DC3S=C13S(K)-C13S(KK)
         DNE=ANE(K)-ANE(KK)
         DMG=AMG(K)-AMG(KK)
         DZO=BETA*S32(K)-BETA*S32(KK)
         DSI=BETA*ASI(K)-BETA*ASI(KK)
         DCA=BETA*CA(K)-BETA*CA(KK)
         DFE=AFE(K)-AFE(KK)
         DCU=CU(K)-CU(KK)
         DZN=ZN(K)-ZN(KK)
         DNI=ANI(K)-ANI(KK)
         DKR=AKR(K)-AKR(KK)

C        new isotopes
         DKP=KP(K)-KP(KK)
         DSC=SC(K)-SC(KK)
         DTI=TI(K)-TI(KK)
         DVANA=VANA(K)-VANA(KK)
         DCR=CR(K)-CR(KK)
         DManga=Manga(K)-Manga(KK)
         DCO=CO(K)-CO(KK)
C     KUNA
         DNA23=ANA23(K)-ANA23(KK)
         DAL27=AL27(K)-AL27(KK)
c     SODIO E ALLUMINIO SECONDARI
         DNAS=ANAS(K)-ANAS(KK)
         DALS=ALS(K)-ALS(KK)        
C        end new isotopes



c ----------------------------------------------------------------------
C Le masse espulse sotto forma dei diversi elementi dipendono ora dalla
C metallicita` iniziale (solo per low and intermediate mass stars).
c ----------------------------------------------------------------------


         IF (ZETAIIN.LE.0.001) THEN
            X1=0.756
            X4=0.243
         ELSE IF (ZETAIIN.GT.0.001.AND.ZETAIIN.LE.0.004) THEN
            X1=0.744
            X4=0.252
         ELSE IF (ZETAIIN.GT.0.004.AND.ZETAIIN.LE.0.008) THEN
            X1=0.728
            X4=0.264
         ELSE IF (ZETAIIN.GT.0.008.AND.ZETAIIN.LE.0.020) THEN
            X1=0.68
            X4=0.30
         ELSE IF (ZETAIIN.GT.0.020) THEN
            X1=0.62
            X4=0.34
         END IF
c----------------------------------------------------------------------------
         
         
         IF(RM(K).EQ.0.) then 
            VMR=0.0
         else
            VMR=RM(K)-DMR*DM/DMM
         endif
         

         VMHE=HE(K)-DMHE*DM/DMM      
         VMC=C12(K)-DMC*DM/DMM
         VMO=O16(K)-DMO*DM/DMM
         VMN=AN14(K)-DMN*DM/DMM
         VMC3=C13(K)-DMC3*DM/DMM
         VMNS=S14(K)-DMNS*DM/DMM
         VC3S=C13S(K)-DC3S*DM/DMM
         VNE=ANE(K)-DNE*DM/DMM
         VMG=AMG(K)-DMG*DM/DMM
         VSI=BETA*ASI(K)-DSI*DM/DMM
         VZO=BETA*S32(K)-DZO*DM/DMM
         VCA=BETA*CA(K)-DCA*DM/DMM
         VFE=AFE(K)-DFE*DM/DMM
c     KUNA
         VNA23=ANA23(K)-DNA23*DM/DMM
         VAL27=AL27(K)-DAL27*DM/DMM
c secondari
         VNAS=ANAS(K)-DNAS*DM/DMM
         VALS=ALS(K)-DALS*DM/DMM
C     VCU=CU(K)-DCU*DM/DMM
C     VZN=ZN(K)-DZN*DM/DMM
         VCU=0.
         VZN=0.
         VNI=0.
         VKR=0.
         VKP=0.
         VSC=0.
         VTI=0.
         VVANA =0.
         VCR=0.
         VManga=0.
         VCO=0.

C
C         "Q" funzione della metallicita' tramite X1 & X4 
         
         QMR=VMR/H
         QMA=QMR+VMHE/H/X1+(VMC+VMO+VMN+VMC3)/H/(X1+X4)
         QMHE=VMHE/H/X1        
         QMC=VMC/H/(X1+X4)      
         QMO=VMO/H/(X1+X4)      
         QMN=VMN/H/(X1+X4)     
         QM3=VMC3/H/(X1+X4)
         QNE=VNE/H/(X1+X4)
         QMG=VMG/H/(X1+X4)
         QSI=VSI/H/(X1+X4)
         QZO=VZO/H/(X1+X4)
         QCA=VCA/H/(X1+X4)
         QFE=VFE/H/(X1+X4)
         QCU=VCU/H/(X1+X4)
         QZN=VZN/H/(X1+X4)
         QNI=VNI/H/(X1+X4)
         QKR=VKR/H/(X1+X4)
C     new isotopes        
         QKP=VKP/H/(X1+X4)
         QSC=VSC/H/(X1+X4)
         QTI=VTI/H/(X1+X4)
         QVANA=VVANA/H/(X1+X4)
         QCR=VCR/H/(X1+X4)
         QManga=VManga/H/(X1+X4)
         QCO=VCO/H/(X1+X4)
C     end new isotopes
C     KUNA - parti primarie
         QNA=VNA23/H/(X1+X4)
         QAL=VAL27/H/(X1+X4)

         X12=0.00361
         X16=0.0108
         X14=0.00111
         
         CL12=X12/10.
         CL14=X14/10.
         CL16=X16/10.    

C Rame e zinco secondari nelle masse intermedie.
C---------------------------------------------------------------- 
         WNS=VMNS/H/(X12+X16)   
         WC3S=VC3S/H/X12       
C----------------------------------------------------------------
    
C NA E AL VENGONO PRODOTTI (ma non da S-PROC) DALLE STELLE LIMS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC     QUI CI VA NA/AL SEC     CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       WNAS=VNAS/H/(CL12+CL14+CL16)
       WALS=VALS/H/(CL12+CL14+CL16)

C eventuali altre componenti secondarie definite in generale zero    
C-------------------------------------------------------------------
C     
         WZNS=0.
         WCUS=0.
         WKRS=0.
         WKPS=0.
         WSCS=0.
         WVANAS=0. 
         WCRS=0.
         WMangas=0.
         WCOS=0.
C-------------------------------------------------------------------

 


C yield secondari di Zn-Kr-Ba-Cu per stelle comprese fra 1 massa e 3 masse solari
C------------------------------------------------------------------------------
         IF(H.GE.1.AND.H.LE.3.) then
            WZNS   = 3.6E-7/H/3.538E-4 
            WKRS   = 1.6E-9/H/3.538E-4          
    
C yield secondario di Cu con Fe/H maggiore di -1
C------------------------------------------------
  
            IF(AFH.GE.-1.0) THEN
               WCUS=2.8E-7/H/3.538E-4 
            ELSE
               WCUS=0.
            END IF
C-----------------------------------------------

         endif 
C-------------------------------------------------------------------------------

         QY=0.
         QSRR=0.
         QZR=0.
         QBas=0.
         QLa=0.


c!!!!!! Tentativo di inserire una nucleosintesi per il processo main s 
c!!!!!! "estratta" da Gallino et al. per gli elementi Sr,Y e Zr          


         If (H.ge.1.and.H.le.3) then
            if (zetaiin.le.0.18e-3) then
               QY  = H*0.107e-7*10**(-0.2)
               QSrr= H*0.518e-7*10**(-0.2)
               QZr = H*0.256e-7*10**(-0.2) 
            else
               if (zeta.le.0.42e-2) then
                  QY = H*0.107e-7*3715*  (zetaiin/0.121e-1)**1.71
                  QSrr=H*0.518e-7*3715*  (zetaiin/0.121e-1)**1.71
                  QZr= H*0.256e-7*3715*  (zetaiin/0.121e-1)**1.71 
               else
                  QY = H*0.107e-7*10**2.2
                  QSrr=H*0.518e-7*10**2.2
                  QZr= H*0.256e-7*10**2.2 
               endif
            endif
            QY = QY /(68*1.5)
            QSrr=QSrr/(68*1.30)
            QZr= QZr /(68.*2.3)
         endif





C     inserisco il valore dello yield del bario secondario  dai dati di busso
C     interpolandoli in funzione della metallicita' delle stelle e della massa
         
         if (H.ge.1.5.and.H.le.3) then
            
            call bario(zetaiin,H,QBAs)
        
         endif




CCCCCCCCCCCCCCCCCCCCCCCC

         if (H.gt.1.and.H.lt.1.5) then
            call bario(zetaiin,1.5,Qbas)
            qbas=qbas/1.5*H
         endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

         if (H.gt.1.5.and.H.lt.3) then
    
            call lanta(zetaiin,H,QLA)

         endif

         if (H.gt.1.and.H.lt.1.5) then
            call lanta(zetaiin,1.5,QLA)
            qla=qla/1.5*H
         endif
    

         

C  aggiungo dei valori costanti per stelle di agb fino a 8
c
c$$$
c$$$         if (H.gt.3.and.H.le.8) then
c$$$                        
c$$$            QBas=.4e-7
c$$$        
c$$$         endif
c$$$         
      endif 
c
c



C-------------------------------------------------------------------------------- -
C  Fine ciclo  per stelle low-intermediate 
C---------------------------------------------------------------------------------
C*********************************************************************************
C---------------------------------------------------------------------------------
    
      RETURN
      END





      subroutine value(zetac,massac,y)
      
      include 'lettura.inc'      
      real zetac,massac
      integer i,icc,irc
      real y,ym(2),yz(2),dy,zetav(2),massav(2)
      
      icc=0
      irc=0

      if (zetac.gt.zetael(col)) then
         zetac=zetael(col)
      endif
      if (zetac.lt.zetael(1)) then
         zetac=zetael(1)
      endif
      if (massac.gt.mass(row)) then
          y=0.
         return
      endif

      if (massac.lt.mass(1)) then
         y=0.
 
         return
      endif
      
          

      do i=1,row
         
         if (mass(i).le.massac) then 
            irc=i
         endif
         
      enddo
      
      do i=1,col
         
         if (zetael(i).le.zetac) then 
            icc=i
         endif
         
      enddo
      
       
      massav(1)=mass(irc)
      massav(2)=mass(irc+1)
      
      zetav(1)=zetael(icc)
      zetav(2)=zetael(icc+1)
      
      
      ym(1)=yield(irc,icc)
      ym(2)=yield(irc+1,icc)
        
        
      call polint(massav,ym,2,massac,y,dy)
      
      yz(1)=y
        
      ym(1)=yield(irc,icc+1)
      
      ym(2)=yield(irc+1,icc+1)
        
        
      call polint(massav,ym,2,massac,y,dy)

              
      yz(2)=y
        
      call polint(zetav,yz,2,zetac,y,dy)
        
           
      return
      
      end
   



     
      subroutine leggi
      implicit none
      include 'lettura.inc'
        
      

      real bainiz,bafin,mhelium
      real bainiz2,bafin2,mhelium2
      integer i,ios,count3
      real batemp,batemp2
      
C Dati in cui si è aggiunto un valore per z=0 pari a
c quello di z=0.00001 ps cambiamento piccolo ma apprezzabile
c una volta sostituite la AMU!!
c
      open (19,file='bario1.dat',status='old')
c
c
c
C Dati originali di Gallino senza valori per z=0
C      open (19,file='bario2.dat',status='old')
C


      count=1
      count3=1
      do i=1,200
         zbario(i)=0.
        
         ba1(i)=0.
         
         ba2(i)=0.
      enddo

C leggi i dati relativi agli yield in funzione di z e della massa

      do 
                 
         read(19,*,iostat=ios) 
     $    zbario(count3),bainiz,bafin,Mhelium,
     $            bainiz2,bafin2,Mhelium2
   
         batemp= (bafin-bainiz)*mhelium
         batemp2=(bafin2-bainiz2)*mhelium2
        
         if (count3.gt.1) then
            if (zbario(count3-1).eq.zbario(count3)) then
               zbario(count)=zbario(count3)
            else
               count=count+1
            endif
         endif
            if (ios.ne.0) exit
            count3=count3+1
        
          
         ba1(count)=batemp+ba1(count)
         ba2(count)=batemp2+ba2(count)
         
      enddo
      close(19)
      
      return
      end




      subroutine leggila
   
      implicit none
      include 'lettura.inc'
      real lainiz, lafin, mhelium
      real lafin2,mhelium2,z
      integer i,ios
      character*5 isot
      
C Dati in cui si è aggiunto un valore per z=0 pari a
c quello di z=0.00001 ps cambiamento piccolo ma apprezzabile
c una volta sostituite la AMU!!
c
      open (19,file='Lanta1.dat',status='old')
c
c
c
C Dati originali di Gallino 

C
C leggi i dati relativi agli yield in funzione di z e della massa
      do i=1,200
         la1(i)=0.
         la2(i)=0.
         zlanta(i)=0.

      enddo
      countla=1
      read(19,*) 
      read(19,*)

      do 
                 
         read(19,*,iostat=ios) 
     $        isot,z, lainiz,lafin,Mhelium,lafin2,Mhelium2
         
         if (isot.eq.'la139') then

            zlanta(countla)=z
            la1(countla)= (lafin-lainiz)*mhelium
            la2(countla)= (lafin2-lainiz)*mhelium2
            countla=countla+1

         endif
     


         if (ios.ne.0) exit
      enddo
      close(19)
      
      return
      end



      subroutine Bario (Zcerc,massa,qba) 
      implicit none
      include 'lettura.inc'
      
      real xd(2),yd(2)
      real xd2(2),yd2(2)
      integer count2,i
      real zcerc,y,dy,qba,qbav(2),massa,massav(2)
      
    
      massav(1)=1.5
      massav(2)=3

C azzera le variabili

      do i=1,2
         xd(i)=0.
         yd(i)=0.
         xd2(i)=0.
         yd2(i)=0.
      enddo
   
      xd(1)=-99999.


C trova in che intervallo di zeta bisogna interpolare i punti 
      
      do count2=1,count-1
         
         if (Zcerc.gt.zbario(count2)) then 
      
            xd(1)=zbario(count2)
            xd(2)=zbario(count2+1)
            yd(1)=ba1(count2)
            yd(2)=ba1(count2+1)
            xd2(1)=zbario(count2)
            xd2(2)=zbario(count2+1)
            yd2(1)=ba2(count2)
            yd2(2)=ba2(count2+1)
      
            
         endif
      enddo
      
C Se non siamo all'interno dell'intervallo in zeta definisce
C zero lo yield (che siamo nell'intervallo giuste nelle masse
C viene definito a monte!)

      if (xd(1).le.-9999) then
      
      qba=0.

      else
         if (xd(1).eq.zbario(count-1)) then 
            
            qba=0.

         else

C
C Altrimenti interpola per entrambe le masse definite
C allo zeta corretto          
     
            call polint(xd,yd,2,zcerc,y,dy)
            qbav(1) = y

            call polint(xd2,yd2,2,zcerc,y,dy)
            qbav(2) = y
C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
            call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
           qba=y
      
         endif
      endif
      
      return
      end


      subroutine Lanta (Zcerc,massa,qla) 
      implicit none
      include 'lettura.inc'
      
      real xd(2),yd(2)
      real xd2(2),yd2(2)
      integer count2,i
      real zcerc,y,dy,qla,qlav(2),massa,massav(2)
      
    
      massav(1)=1.5
      massav(2)=3

C azzera le variabili

      do i=1,2
         xd(i)=0.
         yd(i)=0.
         xd2(i)=0.
         yd2(i)=0.
      enddo
   
      xd(1)=-99999.


C trova in che intervallo di zeta bisogna interpolare i punti 
      
      do count2=1,countla-1
         
         if (Zcerc.gt.zlanta(count2)) then 
      
            xd(1)=zlanta(count2)
            xd(2)=zlanta(count2+1)
            yd(1)=la1(count2)
            yd(2)=la1(count2+1)
            xd2(1)=zlanta(count2)
            xd2(2)=zlanta(count2+1)
            yd2(1)=la2(count2)
            yd2(2)=la2(count2+1)
      
            
         endif
      enddo
      
C Se non siamo all'interno dell'intervallo in zeta definisce
C zero lo yield (che siamo nell'intervallo giuste nelle masse
C viene definito a monte!)

      if (xd(1).le.-9999) then
      
      qla=0.

      else
         if (xd(1).eq.zlanta(countla-1)) then 
            
            qla=0.

         else

C
C Altrimenti interpola per entrambe le masse definite
C allo zeta corretto          
     
            call polint(xd,yd,2,zcerc,y,dy)
            qlav(1) = y

            call polint(xd2,yd2,2,zcerc,y,dy)
            qlav(2) = y
C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
            call polint(massav,qlav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
           qla=y
      
         endif
      endif
      
      return
      end

      
      subroutine  leggi3
      implicit none
      include 'lettura.inc'
      
      integer i,j


     
      
      open (21,file='ZnFe.dat',status='old')
         
      read(21,10) rowzn, colzn 
      rowfe=rowzn
      colfe=colzn
      
 10   format(2i2)
          

      do i=1,colzn
         read(21,*) zetaelzn(i)
         zetaelfe(i)= zetaelzn(i)
      enddo
      
      do i=1,rowzn
         read(21,*) masszn(i)
         massfe(i)= masszn(i)
      enddo

      do i=1,rowzn
         
         do j=1,colzn
            
            read(21,*) yieldzn(i,j), yieldfe(i,j), yieldmn(i,j)
            
         enddo
         
      enddo
    
    

      close(21)


      return
      
      end

      

      subroutine value2(zetac,massac,y)
      
      include 'lettura.inc'      
      real zetac,massac
      integer i,icc,irc
      real y,ym(2),yz(2),dy,zetav(2),massav(2)

      if (zetac.gt.zetaelzn(1)) then
         zetac=zetaelzn(1)
      endif

      if (zetac.lt.zetaelzn(colzn)) then
         zetac=zetaelzn(colzn)
      endif


      if (massac.gt.masszn(rowzn)) then
          y=0.
         return
      endif

      if (massac.lt.masszn(1)) then
         y=0.
         return
      endif
     
        

      do i=1,rowzn
         
         if (masszn(i).le.massac) then 
            irc=i
         endif
         if (irc.eq.rowzn) then
            irc=rowzn-1
         endif
      enddo
      


      do i=1,colzn
         
         if (zetac.lt.zetaelzn(i)) then 
            icc=i
         endif
         if (zetac.eq.zetaelzn(1)) then
            icc=1
         endif
      enddo
      
       
      massav(1)=masszn(irc)
      massav(2)=masszn(irc+1)
      

      zetav(1)=zetaelzn(icc)
      zetav(2)=zetaelzn(icc+1)
      
      
      ym(1)=yieldzn(irc,icc)
      ym(2)=yieldzn(irc+1,icc)
        
        
      call polint(massav,ym,2,massac,y,dy)
      
      yz(1)=y
        
      ym(1)=yieldzn(irc,icc+1)
      
      ym(2)=yieldzn(irc+1,icc+1)
        
        
      call polint(massav,ym,2,massac,y,dy)

              
      yz(2)=y
        
      call polint(zetav,yz,2,zetac,y,dy)
        
           
      return
      
      end
   

      subroutine value3(zetac,massac,y)
      
      include 'lettura.inc'      
      real zetac,massac
      integer i,icc,irc
      real y,ym(2),yz(2),dy,zetav(2),massav(2)

      if (zetac.gt.zetaelfe(1)) then
         zetac=zetaelfe(1)
      endif

      if (zetac.lt.zetaelfe(colfe)) then
         zetac=zetaelfe(colfe)
         y=0.
         return
      endif


      if (massac.gt.massfe(rowfe)) then
          y=0.
         return
      endif

      if (massac.lt.massfe(1)) then
         y=0.
         return
      endif
     
         

      do i=1,rowfe
         
         if (massfe(i).le.massac) then 
            irc=i
         endif
         if (irc.eq.rowfe) then
            irc=rowfe-1
         endif
      enddo
      


      do i=1,colfe
         
         if (zetac.lt.zetaelfe(i)) then 
            icc=i
         endif
         if (zetac.eq.zetaelfe(1)) then
            icc=1
         endif
      enddo
      
       
      massav(1)=massfe(irc)
      massav(2)=massfe(irc+1)
      

      zetav(1)=zetaelfe(icc)
      zetav(2)=zetaelfe(icc+1)
      
      
      ym(1)=yieldfe(irc,icc)
      ym(2)=yieldfe(irc+1,icc)
        
        
      call polint(massav,ym,2,massac,y,dy)
      
      yz(1)=y
        
      ym(1)=yieldfe(irc,icc+1)
      
      ym(2)=yieldfe(irc+1,icc+1)        
        
      call polint(massav,ym,2,massac,y,dy)

              
      yz(2)=y
        
      call polint(zetav,yz,2,zetac,y,dy)
        
           
      return
      
      end
   


      subroutine value4(zetac,massac,y)
      
      include 'lettura.inc'      
      real zetac,massac
      integer i,icc,irc
      real y,ym(2),yz(2),dy,zetav(2),massav(2)

      if (zetac.gt.zetaelfe(1)) then
         zetac=zetaelfe(1)
      endif

      if (zetac.lt.zetaelfe(colfe)) then
         zetac=zetaelfe(colfe)
         y=0.
         return
      endif


      if (massac.gt.massfe(rowfe)) then
          y=0.
         return
      endif

      if (massac.lt.massfe(1)) then
         y=0.
         return
      endif
     
         

      do i=1,rowfe
         
         if (massfe(i).le.massac) then 
            irc=i
         endif
         if (irc.eq.rowfe) then
            irc=rowfe-1
         endif
      enddo
      


      do i=1,colfe
         
         if (zetac.lt.zetaelfe(i)) then 
            icc=i
         endif
         if (zetac.eq.zetaelfe(1)) then
            icc=1
         endif
      enddo
      
       
      massav(1)=massfe(irc)
      massav(2)=massfe(irc+1)
      

      zetav(1)=zetaelfe(icc)
      zetav(2)=zetaelfe(icc+1)
      
      
      ym(1)=yieldmn(irc,icc)
      ym(2)=yieldmn(irc+1,icc)
        
        
      call polint(massav,ym,2,massac,y,dy)
      
      yz(1)=y
        
      ym(1)=yieldmn(irc,icc+1)
      
      ym(2)=yieldmn(irc+1,icc+1)
        
        
      call polint(massav,ym,2,massac,y,dy)

              
      yz(2)=y
        
      call polint(zetav,yz,2,zetac,y,dy)
        
           
      return
      
      end
   




 
C
C     routine di interpolazione polinomiale tratta tra numerical r.
C     (usata solo la interpolazione lineare in reltà n=2!)
C
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=1000)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
c
c

      SUBROUTINE CALWI
      
      PARAMETER (NMAX=37)      
      include 'lettura.inc'
      include 'elementi.inc'
      REAL MI,MS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
      DIMENSION RMAS(3000),RNUM(3000),RNOVAE(3000)
      DIMENSION RNOV(NMAX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei raggi cosmici.
      DIMENSION COSM(NMAX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DIMENSION GRAT(3000)
      DIMENSION VN(NMAX),PRZW1(NMAX),PRZW2(NMAX),PRZW3(NMAX)
      DIMENSION TW(NMAX),RR1(NMAX),RR2(NMAX),PARZW(NMAX)
      DIMENSION QM(NMAX,NMAX),RR3(NMAX),RR4(NMAX)
      DIMENSION XP(NMAX),R1(NMAX),IND(NMAX),XSTIN(NMAX)
      DIMENSION FALL(3000),AMU(200),GO(3000),G(3000,NMAX)
      DIMENSION XTIN(NMAX),XCAL(NMAX),XSTAR(NMAX),PAMU(200)
      DIMENSION APUNT(3000),WI1(NMAX),WI2(NMAX)
      DIMENSION T(3000),GP(NMAX),WI(NMAX),X(3000,NMAX),SRT(3000)
      DIMENSION SMS(3000)
      DIMENSION RATE(3000),XINF(NMAX),FA(NMAX)
      DIMENSION COSN1(50),COSN2(50)
      COMMON/CUATRO/COSN1,COSN2           
      COMMON/UNO/GO,G,T,X,XCAL,XSTAR,XTIN,XSTIN
      COMMON/FR/XINF,FALL,FA,SG
      COMMON/TRE/AMU,GP,XP,WI,PARZW,WI1,WI2,PRZW1,PRZW2,PRZW3
      COMMON/A/QM
      COMMON/DUE/CK,UM,TOL,EPSI,TAUINF,VDTEM,DTEM,SRP,
     $     GASP,GPRE,FAP,FPRE,RAP,RPRE,SINT,CNORM,MI,MS,M,IMAX,
     $     KMAX,I,ITER,IND,KDEMAX
      COMMON/COM/MAXO
      COMMON/DIE/A,B,PAMU,K2,K1,KB,KMAX2,K3,KB3
      COMMON/DIS/EXPO,T1,T2,ETA,RO,RS,AR,BR,RI,RD,SD,SN,AO,
     $     R,SMO,SMC,R1,VMR,RN,VNSOL,SRT,ARO,BRO,SMS,GAS,APUNT,
     $     RATE,RIEST,RSCA,AC,TW,VN,THALO
      COMMON/DATA/LET
c     !COMMON/LUCKY/IKU
      COMMON/M/NMAT,NNMAT
      COMMON/SN/RR1,RR2,RR3,RR4,COST
      COMMON/NEW/ASNII,GAMMA
      COMMON/PIP/GRAT,DTGRAT
      COMMON/RAF/TAUIN,TNEW,gr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
      COMMON/NOVA/RMAS,RNUM,RNOVAE,WDM,WD,RANOV
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei raggi cosmici.
      COMMON/RAYS/COSM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real tmedio

CCC Initialization variables
     
      tmedio=0.
      rawd1=0.
      rawdm1=0.
      tmt=0.
      rawd3=0.
      rawdm3=0
      asn1=0.

CCCCC     
 
      IF (ITER.ne.1) then
         
         TMEDIO=T(I-1)+DTEM/2.
         
         TNEW1=TMEDIO-TNEW
         
         A=ASSA(TNEW1)

         
         IF(A.ge.AMU(KMAX)) then
            
            DO IJ=1,NMAX
               WI(IJ)=FAP*XINF(IJ)
            enddo

C     return because no contribution in this case
C     the interval of integration is always zero!

            return
C     
         endif
         
         B=ASSA(DTEM/2.)
         IF(B.ge.AMU(KMAX)) then 
            B=AMU(KMAX)
         endif

         
         DO K=1,KMAX
            C=AMU(K)-A
            IF(C.GT.0.) exit
         enddo
         K1=K
         
         
         DO K=1,KMAX
            C=AMU(K)-B
            IF(C.GT.0.) exit
         enddo
         
         
         IF(K.GT.KMAX) K=KMAX
         K2=K-1
         
         DO K=1,KMAX
            E=AMU(K)-8.
            IF(E.GT.0.) exit
         enddo
         
         K3=K-1
         IF(K2.lt.K3) then
            
            if (K2.lt.K1) then 
               PAMU(1)=A
               PAMU(2)=B
               PAMU(3)=8.
               L4=KMAX-K3+1
               DO LL=1,L4
                  PAMU(3+LL)=AMU(K3-1+LL)
               enddo
               KMAX2=KMAX-K3+4
               KB=2
               KB3=3
            else
               
               PAMU(1)=A
               L1=K2-K1+1
               DO L=1,L1
                  PAMU(L+1)=AMU(K1+L-1)
               enddo
               PAMU(K2-K1+3)=B
               KB=K2-K1+3
               L2=K3-K2
               DO LL=1,L2
                  PAMU(K2-K1+3+LL)=AMU(K2+LL)
               enddo
               PAMU(K3-K1+3)=8.
               KB3=K3-K1+3
               L3=KMAX-K3
               DO LL=1,L3
                  PAMU(K3-K1+3+LL)=AMU(K3+LL)
               enddo
               KMAX2=KMAX-K1+3
            endif
         else
            if (K1.le.K3) then
               
               PAMU(1)=A
               L1=K3-K1+1
               DO L=1,L1
                  PAMU(L+1)=AMU(K1+L-1)
               enddo
               
               KB3=K3-K1+2
               L2=K2-K3
               DO LL=1,L2
                  PAMU(K3-K1+3+LL)=AMU(K3+LL)
               enddo
               PAMU(K2-K1+3)=B
               KB=K2-K1+3
               L2=KMAX-K2
               DO LL=1,L2
                  PAMU(K2-K1+3+LL)=AMU(K2+LL)
               enddo
               KAKA=0
               IF(B.EQ.AMU(KMAX)) KAKA=-1
               KMAX2=KMAX-K1+3+KAKA
            else
               PAMU(1)=A
               L1=K2-K1+1
               DO L=1,L1
                  PAMU(L+1)=AMU(K1+L-1)
               enddo
               PAMU(K2-K1+3)=B
               KB=K2-K1+3
               L2=KMAX-2
               DO LL=1,L2
                  PAMU(K2-K1+3+LL)=AMU(K2+LL)
               enddo
               KAKA=0
               IF(B.EQ.AMU(KMAX)) KAKA=-1
               KMAX2=KMAX-K1+3+KAKA
            endif
         endif


      endif


c!      WRITE(60,100) A,B,I
      



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
C     ALFA=frazione di WDs che si trova in sistemi binari che originano
C     esplosioni di nova;
C     XLI1=frazione di massa espulsa sotto forma di 7Li;
C     XC13=frazione di massa espulsa sotto forma di 13C;
C     XC12=frazione di massa espulsa sotto forma di 12C;
C     XN15=frazione di massa espulsa sotto forma di 15N;
C     XN14=frazione di massa espulsa sotto forma di 14N;
C     XO16=frazione di massa espulsa sotto forma di 16O;
C     AVEJ=massa media espulsa da un sistema di nova durante un outburst
C     moltiplicata per il numero medio di ricorrenze dell'evento di nova.

      ALFA=0.0115
      BBETA=0.429
      GGAMM=0.243
      DDELT=0.188
      SETA=0.14
      IF (T(I).LE.0.045) THEN
C     Media su tutte le sequenze evolutive:
         XLI1=9.24E-7
C     XC13=2.38E-2
C     XC12=2.16E-2
C     XN15=5.34E-2
C     XN14=3.11E-2
C     XO16=7.28E-1
C     Minima produzione di 13C:
         XC13=1.5E-2
         XC12=2.1E-2
         XN15=1.2E-1
         XN14=4.6E-2
         XO16=2.2E-2
      ELSE
C     Media su tutte le sequenze evolutive:
         XLI1=2.85E-6
C     XC13=7.75E-2
C     XC12=3.45E-2
C     XN15=2.66E-2
C     XN14=8.69E-2
C     XO16=1.26
C     Minima produzione di 13C:
         XC13=2.83E-2
         XC12=1.61E-2
         XN15=3.67E-2
         XN14=8.03E-2
         XO16=9.76E-2
      END IF

C     IF (T(I).LE.0.045) THEN
C     Media su tutte le sequenze evolutive:
C     AVEJ=1.95E-1
C     ELSE
C     AVEJ=2.63E-1
C     END IF
C     Minima massa espulsa:
      AVEJ=1.95E-1
      AMLI1=XLI1*AVEJ
      AMC13=XC13*AVEJ
      AMC12=XC12*AVEJ
      AMN15=XN15*AVEJ
      AMN14=XN14*AVEJ
      AMO16=XO16*AVEJ
C     BBETA e` la frazione di sistemi di nova con delay time di 1 Gyr;
C     GGAMM e` la frazione di sistemi di nova con delay time di 2.5 Gyr;
C     DDELT e` la frazione di sistemi di nova con delay time di 3.5 Gyr;
C     SETA e` la frazione di sistemi di nova con delay time di 5 Gyr.
      TDEL1=2.0
      TDEL2=2.0
      TDEL3=2.0
      TDEL4=2.0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF(A.LE.8.)CALL SUPER

      DO J=1,NMAX
         
         IF(A.le.8.) then
            
            IF(B.le.8.) then

               IF(ITER.ne.1) then

                  L4=KB-1
                  SIM1=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
                  RAWD1=0.
                  RAWDM1=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  
                  DO JJ=1,L4
                     AMED=(PAMU(JJ)+PAMU(JJ+1))*0.5
                     
                     TME=TAU(AMED)
                     TMT=TMEDIO-TME
                     TMT=ABS(TMT)

                     IF (I.EQ.1) THEN
                        ISUP=1
                        GINF=GO(ISUP)
                        DO L=1,NMAX
                           XTIN(L)=X(ISUP,L)
                        ENDDO
                        TIS=T(ISUP)
                        RAINF=RATE(ISUP)*GINF**CK

                     ELSE
                        
                        DO  KO=1,I
                           IF(TMT.LT.T(KO)) exit
                        enddo
                       
                        
                        ISUP=KO
                        INF=KO-1
                        DTMT=TMT-T(INF)
                        TIS=T(ISUP)
                        TIN=T(INF)
                        DT=TIS-TIN
                        DG=GO(ISUP)-GO(INF)
                        GINF=GO(INF)+DG/DT*DTMT

                        DO L=1,NMAX
                           XTIN(L)=X(INF,L)+(X(ISUP,L)-X(INF,L))/DT*DTMT
                        enddo
                        
                        IF(TMT.LT.TAUIN) THEN
                           br=0.
                           bro=0.
                        else
                           br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
                           bro=br
                        end if
                        
                        SRTM=AR*T1*(1.-EXP(-TMT/T1))+BR*T2*(1.-EXP(-
     $                       (TMT-TAUIN)/T2))+GAS

                        SMSM=ARO*T1*(1.-EXP(-TMT/T1))+BRO*T2*(1.-EXP(-
     $                       (TMT-TAUIN)/T2))  +GAS
                        GRATM=(GRAT(INF)+GRAT(ISUP))*0.5
                        RAINF=GRATM*(SRTM/SMSM)**EXPO*
     $                       (SG/SRTM)**(CK-1.)*GINF**CK
                     endif
                     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
                     ZETAIN=0.
                     DO L=5,NMAX
                        ZETAIN=ZETAIN+XTIN(L)     
                     enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                     SIJ=0.0
                     
                     DO JS=1,NMAX
                        SIJ=SIJ+Q(AMED,J,JS,ZETAIN)*XTIN(JS)
                     enddo
                     
                     XL=TAU(PAMU(JJ))
                     XU=TAU(PAMU(JJ+1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
C     Rate di WDs in numero al tempo T:
                     RAWD1=RAWD1+RAINF*FII(XL,XU)
C     Rate di WDs in massa al tempo T:
                     RAWDM1=RAWDM1+RAINF*FI(XL,XU)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     
                     
                     SIM1=SIM1+SIJ*RAINF*FI(XL,XU)
                  enddo


                  PRZW1(J)=SIM1



               endif
               L5=KB3-1
               SIM2=0.0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
               RAWD2=0.
               RAWDM2=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               DO JJ=KB,L5
                  TDIZ=T(I)
                  ASTR=(PAMU(JJ)+PAMU(JJ+1))*0.5
                  TSTAR=TMEDIO-TAU(ASTR)
                  GSTAR=GO(I)+(GASP-GO(I))/DTEM*(TSTAR-TDIZ)
                  DO L=1,NMAX
C     HO SOSTITUITO XSTAR
                     XTIN(L)=X(I,L)+(XP(L)-X(I,L))/DTEM*(TSTAR-TDIZ)
                  ENDDO
                  
                  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
                  ZETAIN=0.
                  DO L=5,NMAX
                     ZETAIN=ZETAIN+XTIN(L)
                  ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  
                  PX=TAU(PAMU(JJ))
                  PY=TAU(PAMU(JJ+1))
                  SIJ1=0.0
                  DO JS=1,NMAX
                     SIJ1=SIJ1+Q(ASTR,J,JS,ZETAIN)*XTIN(JS)
                  ENDDO
                  IF(TSTAR.LT.TAUIN) THEN
                     br=0.
                     bro=0.
                  else
                     br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
                     bro=br
                  end if
                  SSTAR=AR*T1*(1.-EXP(-TSTAR/T1))+BR*T2*
     $                 (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                 GAS
                  SMTAR=ARO*T1*(1.-EXP(-TSTAR/T1))+BRO*T2*
     $                 (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                 GAS
                  RAINF=GRAT(I)*(SSTAR/SMTAR)**EXPO*(SG/SMTAR)**(CK-1.)*
     $                 GSTAR**CK
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
                  RAWD2=RAWD2+RAINF*FII(PX,PY)
                  RAWDM2=RAWDM2+RAINF*FI(PX,PY)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  SIM2=SIM2+RAINF*SIJ1*FI(PX,PY)
               ENDDO
               



               PRZW2(J)=SIM2
               L6=KMAX2-1
               SIM3=0.
               ASN2=0.
               DO LL=KB3,L6
                  TDIZ=T(I)
                  ASTR=(PAMU(LL)+PAMU(LL+1))*0.5
                  TSTAR=TMEDIO-TAU(ASTR)
                  GSTAR=GO(I)+(GASP-GO(I))/DTEM*(TSTAR-TDIZ)
                  DO L=1,NMAX
                     XTIN(L)=X(I,L)+(XP(L)-X(I,L))/DTEM*(TSTAR-TDIZ)
                  enddo
                  SSTAR=AR*T1*(1.-EXP(-TSTAR/T1))+BR*T2*
     $                 (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                 GAS
                  SMTAR=ARO*T1*(1.-EXP(-TSTAR/T1))+BRO*T2*
     $                 (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                 GAS
                  RAINF=GRAT(I)*(SSTAR/SMTAR)**EXPO*
     $                 (SG/SSTAR)**(CK-1.)*GSTAR**CK
                  PX=TAU(PAMU(LL))
                  PY=TAU(PAMU(LL+1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
                  ZETAIN=0.
                  DO L=5,NMAX
                     ZETAIN=ZETAIN+XTIN(L)
                  ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  SIJ2=0.
                  DO  JS=1,NMAX
                     SIJ2=SIJ2+Q(ASTR,J,JS,ZETAIN)*XTIN(JS)
                  ENDDO
                  
                  ASN2=ASN2+RAINF*FII(PX,PY)
                  SIM3=SIM3+RAINF*SIJ2*FI(PX,PY)
               ENDDO
               
               PRZW3(J)=SIM3
               ASNII=ASN2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
               WD=RAWD1+RAWD2
               WDM=RAWDM1+RAWDM2

               RM1A=0.
               RWD1=0.
               RM2A=0.
               RWD2=0.
               RM3A=0.
               RWD3=0.
               RM4A=0.
               RWD4=0.

               
               IF(TMEDIO.GT.TDEL1) then
                  TNOV1=TMEDIO-TDEL1
                  TNOV1=ABS(TNOV1)
                  DO KO1=19,I
                     IF(TNOV1.lt.T(KO1)) exit
                  ENDDO
                  ISO1=KO1
                  INO1=KO1-1
                  DTNOV1=TNOV1-T(INO1)
                  DTN1=T(ISO1)-T(INO1)
                  RM1A=RMAS(INO1)+(RMAS(ISO1)-RMAS(INO1))/DTN1*DTNOV1
                  RWD1=RNUM(INO1)+(RNUM(ISO1)-RNUM(INO1))/DTN1*DTNOV1
                  
             
                  IF(TMEDIO.GT.TDEL2) then
                     TNOV2=TMEDIO-TDEL2
                     TNOV2=ABS(TNOV2)
                     DO KO2=19,I
                        IF(TNOV2.lt.T(KO2)) exit
                     enddo
                     ISO2=KO2
                     INO2=KO2-1
                     DTNOV2=TNOV2-T(INO2)
                     DTN2=T(ISO2)-T(INO2)
                     RM2A=RMAS(INO2)+(RMAS(ISO2)-RMAS(INO2))/DTN2*DTNOV2
                     RWD2=RNUM(INO2)+(RNUM(ISO2)-RNUM(INO2))/DTN2*DTNOV2
                     
                     
                     IF(TMEDIO.gt.TDEL3) then
                        TNOV3=TMEDIO-TDEL3
                        TNOV3=ABS(TNOV3)
                        DO KO3=19,I
                           IF(TNOV3.lt.T(KO3)) exit
                        enddo
                        ISO3=KO3
                        INO3=KO3-1
                        DTNOV3=TNOV3-T(INO3)
                        DTN3=T(ISO3)-T(INO3)
                        RM3A=RMAS(INO3)+(RMAS(ISO3)-
     $                       RMAS(INO3))/DTN3*DTNOV3
                        RWD3=RNUM(INO3)+(RNUM(ISO3)-
     $                       RNUM(INO3))/DTN3*DTNOV3
                        
                        IF(TMEDIO.gt.TDEL4) then
                           TNOV4=TMEDIO-TDEL4
                           TNOV4=ABS(TNOV4)
                           DO KO4=19,I
                              IF(TNOV4.lt.T(KO4)) exit
                           enddo
                           ISO4=KO4
                           INO4=KO4-1
                           DTNOV4=TNOV4-T(INO4)
                           DTN4=T(ISO4)-T(INO4)
                           RM4A=RMAS(INO4)+(RMAS(ISO4)-
     $                          RMAS(INO4))/DTN4*DTNOV4
                           RWD4=RNUM(INO4)+(RNUM(ISO4)-
     $                          RNUM(INO4))/DTN4*DTNOV4
                           
                        end if
                        
                     endif
                     
                  endif
                  
               endif
               

C     Rate di novae in numero al tempo T:
               
               RANOV=ALFA*(RWD1*BBETA+RWD2*GGAMM+RWD3*DDELT+RWD4*SETA)
               RN1=BBETA*RM1A*AMLI1
               RN2=GGAMM*RM2A*AMLI1
               RN3=DDELT*RM3A*AMLI1
               RN4=SETA*RM4A*AMLI1
               RC131=BBETA*RM1A*AMC13
               RC132=GGAMM*RM2A*AMC13
               RC133=DDELT*RM3A*AMC13
               RC134=SETA*RM4A*AMC13
               RC121=BBETA*RM1A*AMC12
               RC122=GGAMM*RM2A*AMC12
               RC123=DDELT*RM3A*AMC12
               RC124=SETA*RM4A*AMC12
               RN151=BBETA*RM1A*AMN15
               RN152=GGAMM*RM2A*AMN15
               RN153=DDELT*RM3A*AMN15
               RN154=SETA*RM4A*AMN15
               RN141=BBETA*RM1A*AMN14
               RN142=GGAMM*RM2A*AMN14
               RN143=DDELT*RM3A*AMN14
               RN144=SETA*RM4A*AMN14
               RO161=BBETA*RM1A*AMO16
               RO162=GGAMM*RM2A*AMO16
               RO163=DDELT*RM3A*AMO16
               RO164=SETA*RM4A*AMO16
               DO KKK=1,NMAX
                  RNOV(KKK)=0.
               enddo
               RNOV(5)=ALFA*(RC121+RC122+RC123+RC124)
               RNOV(6)=ALFA*(RO161+RO162+RO163+RO164)
               RNOV(7)=ALFA*(RN141+RN142+RN143+RN144)
               RNOV(8)=ALFA*(RC131+RC132+RC133+RC134)
               RNOV(20)=ALFA*(RN1+RN2+RN3+RN4)
               RNOV(21)=ALFA*(RN151+RN152+RN153+RN154)
               
               CALL COSMIC
               WI1(J)=RR2(J)+PRZW1(J)+PRZW2(J)+PRZW3(J)+RNOV(J)
     $              +COSM(J)*GO(I)*SG*X(I,1)*.8

C     Se si escludono le novae WI1(J) va calcolato con la formula seguente:
C     WI1(J)=RR2(J)+PRZW1(J)+PRZW2(J)+PRZW3(J)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               WI2(J)=XINF(J)*FAP
               WI(J)=WI1(J)+WI2(J)



            else

               IF(ITER.ne.1) then
                  

                  L4=KB3-1
                  SOM1=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
                  RAWD3=0.
                  RAWDM3=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  
                  DO LL=1,L4
                     AMED=(PAMU(LL)+PAMU(LL+1))*0.5
                     TME=TAU(AMED)
                     TMT=ABS(TMEDIO-TME)
                     IF(I.EQ.1) then 
                        
                        ISUP=1
                        GINF=GO(ISUP)
                        DO LK=1,NMAX
                           XTIN(LK)=X(ISUP,LK)
                        enddo
                        TIS=T(ISUP)
                        RAINF=RATE(ISUP)*GINF**CK
                     else
                        DO  KO=1,I
                           IF(TMT.lt.T(KO)) exit
                        enddo
                        
                        
                        ISUP=KO
                        INF=KO-1
                        DTMT=TMT-T(INF)
                        TIS=T(ISUP)
                        TIN=T(INF)
                        DT=TIS-TIN
                        DG=GO(ISUP)-GO(INF)
                        GINF=GO(INF)+DG/DT*DTMT
                        DO L=1,NMAX
                           XTIN(L)=X(INF,L)+(X(ISUP,L)-X(INF,L))/DT*DTMT
                        enddo

                        IF(TMT.LT.TAUIN) THEN
                           br=0.
                           bro=0.
                        else
                           br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
                           bro=br
                        end if
                        SRTM=AR*T1*(1.-EXP(-TMT/T1))+BR*T2*(1.-EXP(-
     $                       (TMT-TAUIN)/T2)) +GAS
                        SMSM=ARO*T1*(1.-EXP(-TMT/T1))+BRO*T2*(1.-EXP(-
     $                       (TMT-TAUIN)/T2))
     $                       +GAS
                        GRATM=(GRAT(INF)+GRAT(ISUP))*0.5
                        RAINF=GRATM*(SRTM/SMSM)**EXPO*
     $                       (SG/SRTM)**(CK-1.)*GINF**CK
                     endif
                     
                     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
                     ZETAIN=0.
                     DO L=5,NMAX
                        ZETAIN=ZETAIN+XTIN(L)
                     enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     SIJ=0.
                     DO  JS=1,NMAX
                        SIJ=SIJ+Q(AMED,J,JS,ZETAIN)*XTIN(JS)
                     enddo
                     
                     XL=TAU(PAMU(LL))
                     XU=TAU(PAMU(LL+1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
                     RAWD3=RAWD3+RAINF*FII(XL,XU)
                     RAWDM3=RAWDM3+RAINF*FI(XL,XU)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     SOM1=SOM1+SIJ*RAINF*FI(XL,XU)
                  enddo
                  PRZW1(J)=SOM1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
C     Calcolo rate rest in pezzi.
                  WDM=RAWDM3
                  WD=RAWD3
                  RM13=0.
                  RWD13=0.
                  RM23=0.
                  RWD23=0.
                  RM33=0.
                  RWD33=0.
                  RM43=0.
                  RWD43=0.

                  
                  IF(TMEDIO.gt.TDEL1) then
                     
                     TNOV1=TMEDIO-TDEL1
                     TNOV1=ABS(TNOV1)
                     DO KO1=19,I
                        IF(TNOV1.lt.T(KO1)) exit
                     enddo
                     
                     ISO1=KO1
                     INO1=KO1-1
                     DTNOV1=TNOV1-T(INO1)
                     DTN1=T(ISO1)-T(INO1)
                     RM13=RMAS(INO1)+(RMAS(ISO1)-
     $                    RMAS(INO1))/DTN1*DTNOV1
                     RWD13=RNUM(INO1)+(RNUM(ISO1)-
     $                    RNUM(INO1))/DTN1*DTNOV1
                     
                     IF(TMEDIO.gt.TDEL2) then
                        TNOV2=TMEDIO-TDEL2
                        TNOV2=ABS(TNOV2)
                        DO KO2=19,I
                           IF(TNOV2.lt.T(KO2)) exit
                        enddo
                        
                        
                        ISO2=KO2
                        INO2=KO2-1
                        DTNOV2=TNOV2-T(INO2)
                        DTN2=T(ISO2)-T(INO2)
                        RM23=RMAS(INO2)+(RMAS(ISO2)-
     $                       RMAS(INO2))/DTN2*DTNOV2
                        RWD23=RNUM(INO2)+(RNUM(ISO2)-
     $                       RNUM(INO2))/DTN2*DTNOV2
                        
                        IF(TMEDIO.gt.TDEL3) then
                           TNOV3=TMEDIO-TDEL3
                           TNOV3=ABS(TNOV3)
                           
                           DO KO3=19,I
                              IF(TNOV3.lt.T(KO3)) exit
                           enddo
                           
                           
                           ISO3=KO3
                           INO3=KO3-1
                           DTNOV3=TNOV3-T(INO3)
                           DTN3=T(ISO3)-T(INO3)
                           RM33=RMAS(INO3)+(RMAS(ISO3)-
     $                          RMAS(INO3))/DTN3*DTNOV3
                           RWD33=RNUM(INO3)+(RNUM(ISO3)-
     $                          RNUM(INO3))/DTN3*DTNOV3
                           
                           IF(TMEDIO.gt.TDEL4) then
                              
                              TNOV4=TMEDIO-TDEL4
                              TNOV4=ABS(TNOV4)
                              
                              DO KO4=19,I
                                 IF(TNOV4.lt.T(KO4)) exit
                              enddo
                              ISO4=KO4
                              INO4=KO4-1
                              DTNOV4=TNOV4-T(INO4)
                              DTN4=T(ISO4)-T(INO4)
                              RM43=RMAS(INO4)+(RMAS(ISO4)-
     $                             RMAS(INO4))/DTN4*DTNOV4
                              RWD43=RNUM(INO4)+(RNUM(ISO4)-
     $                             RNUM(INO4))/DTN4*DTNOV4
                              
                           endif
                        endif
                     endif
                  endif
                  
                  
C     Rate di novae in numero al tempo T:
                  
                  RANOV=ALFA*(RWD13*BBETA+
     $                 RWD23*GGAMM+RWD33*DDELT+RWD43*SETA)
                  RN13=BBETA*RM13*AMLI1
                  RN23=GGAMM*RM23*AMLI1
                  RN33=DDELT*RM33*AMLI1
                  RN43=SETA*RM43*AMLI1
                  RC1313=BBETA*RM13*AMC13
                  RC1323=GGAMM*RM23*AMC13
                  RC1333=DDELT*RM33*AMC13
                  RC1343=SETA*RM43*AMC13
                  RC1213=BBETA*RM13*AMC12
                  RC1223=GGAMM*RM23*AMC12
                  RC1233=DDELT*RM33*AMC12
                  RC1243=SETA*RM43*AMC12
                  RN1513=BBETA*RM13*AMN15
                  RN1523=GGAMM*RM23*AMN15
                  RN1533=DDELT*RM33*AMN15
                  RN1543=SETA*RM43*AMN15
                  RN1413=BBETA*RM13*AMN14
                  RN1423=GGAMM*RM23*AMN14
                  RN1433=DDELT*RM33*AMN14
                  RN1443=SETA*RM43*AMN14
                  RO1613=BBETA*RM13*AMO16
                  RO1623=GGAMM*RM23*AMO16
                  RO1633=DDELT*RM33*AMO16
                  RO1643=SETA*RM43*AMO16
                  
                  DO  KJK=1,NMAX
                     RNOV(KJK)=0.
                  enddo
                  
                  RNOV(5)=ALFA*(RC1213+RC1223+RC1233+RC1243)
                  RNOV(6)=ALFA*(RO1613+RO1623+RO1633+RO1643)
                  RNOV(7)=ALFA*(RN1413+RN1423+RN1433+RN1443)
                  RNOV(8)=ALFA*(RC1313+RC1323+RC1333+RC1343)
                  RNOV(20)=ALFA*(RN13+RN23+RN33+RN43)
                  RNOV(21)=ALFA*(RN1513+RN1523+RN1533+RN1543)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  


                  KBS=KB3
                  L5=KB-1
                  ASN1=0.
                  SOM2=0.


                  DO LL=KBS,L5
                     
                     AMED=(PAMU(LL)+PAMU(LL+1))*0.5
                     TME=TAU(AMED)
                     TMT=ABS(TMEDIO-TME)
                     IF(I.EQ.1) then
                        ISUP=1
                        GINF=GO(ISUP)
                        DO LK=1,NMAX
                           XTIN(LK)=X(ISUP,LK)
                        enddo
                        TIS=T(ISUP)
                        RAINF=RATE(ISUP)*GINF**CK
                     else
                        
                        DO KO=1,I
                           IF(TMT.lt.T(KO)) exit
                        enddo
                        

                        
                        ISUP=KO
                        INF=KO-1
                        
                        DTMT=TMT-T(INF)
                        TIS=T(ISUP)
                        TIN=T(INF)
                        DT=TIS-TIN
                        DG=GO(ISUP)-GO(INF)
                        GINF=GO(INF)+DG/DT*DTMT
                        DO L=1,NMAX
                           XTIN(L)=X(INF,L)+(X(ISUP,L)-X(INF,L))/DT*DTMT
                        enddo

                        IF(TMT.LT.TAUIN) THEN
                           br=0.
                           bro=0.
                        else
                           br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
                           bro=br
                        end if
                        SRTM=AR*T1*(1.-EXP(-TMT/T1))+BR*T2*(1.-
     $                       EXP(-(TMT-TAUIN)/T2)) +GAS
                        SMSM=ARO*T1*(1.-EXP(-TMT/T1))+BRO*T2*(1.-
     $                       EXP(-(TMT-TAUIN)/T2)) +GAS
                        GRATM=(GRAT(INF)+GRAT(ISUP))*0.5
                        RAINF=GRATM*(SRTM/SMSM)**EXPO*
     $                       (SG/SRTM)**(CK-1.)*GINF**CK
                     endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
                     ZETAIN=0.
                     DO  L=5,NMAX
                        ZETAIN=ZETAIN+XTIN(L)
                     enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     SIJ1=0.
                     DO  JS=1,NMAX
                        SIJ1=SIJ1+Q(AMED,J,JS,ZETAIN)*XTIN(JS)
                     enddo
                     XL=TAU(PAMU(LL))
                     XU=TAU(PAMU(LL+1))
                     ASN1=ASN1+RAINF*FII(XL,XU)
                     SOM2=SOM2+SIJ1*RAINF*FI(XL,XU)
                  enddo



                  PRZW2(J)=SOM2


               endif

               ASN2=0.
               SOM3=0.
               L6=KMAX2-1


               IF(B.ne.AMU(KMAX)) then

                  
                  DO LL=KB,L6
                     TDIZ=T(I)
                     ASTR=(PAMU(LL)+PAMU(LL+1))*0.5
                     TSTAR=TMEDIO-TAU(ASTR)
                     GSTAR=GO(I)+(GASP-GO(I))/DTEM*(TSTAR-TDIZ)
                     
                     DO L=1,NMAX
                        XTIN(L)=X(I,L)+(XP(L)-X(I,L))/DTEM*(TSTAR-TDIZ)
                     enddo
                     
                     IF(TMT.LT.TAUIN) THEN
                        br=0.
                        bro=0.
                     else
                        br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
                        bro=br
                     end if
                     SSTAR=AR*T1*(1.-EXP(-TSTAR/T1))+BR*T2*
     $                    (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                    GAS
                     SMTAR=ARO*T1*(1.-EXP(-TSTAR/T1))+BRO*T2*
     $                    (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                    GAS
                     RAINF=GRAT(I)*(SSTAR/SMTAR)**EXPO*
     $                    (SG/SSTAR)**(CK-1.)*GSTAR**CK
                     PX=TAU(PAMU(LL))
                     PY=TAU(PAMU(LL+1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
                     ZETAIN=0.
                     DO  L=5,NMAX
                        ZETAIN=ZETAIN+XTIN(L)
                     enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     SIJ2=0.
                     DO JS=1,NMAX
                        SIJ2=SIJ2+Q(ASTR,J,JS,ZETAIN)*XTIN(JS)
                     enddo
                     ASN2=ASN2+RAINF*FII(PX,PY)

                     SOM3=SOM3+RAINF*SIJ2*FI(PX,PY)
                  enddo
                  PRZW3(J)=SOM3




               endif
               
               ASNII=ASN1+ASN2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.

               WD=RAWD3
               WDM=RAWDM3
               


               CALL COSMIC
               WI1(J)=RR1(J)+PRZW1(J)+PRZW2(J)+PRZW3(J)+RNOV(J)
     $              +COSM(J)*GO(I)*SG*X(I,1)*.8

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               WI2(J)=XINF(J)*FAP
               WI(J)=WI1(J)+WI2(J)

            endif

C     

         else


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     Parte se A è falso
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




            KBS=1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.
            RAWD3=0.
            RAWDM3=0.
            RANOV=0.
            DO KIJ=1,NMAX
               RNOV(KIJ)=0.
            enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            L5=KB-1
            ASN1=0.
            SOM2=0.

            DO LL=KBS,L5
               
               AMED=(PAMU(LL)+PAMU(LL+1))*0.5
               TME=TAU(AMED)
               TMT=ABS(TMEDIO-TME)
               IF(I.EQ.1) then
                  ISUP=1
                  GINF=GO(ISUP)
                  DO LK=1,NMAX
                     XTIN(LK)=X(ISUP,LK)
                  enddo
                  TIS=T(ISUP)
                  RAINF=RATE(ISUP)*GINF**CK
               else
                  
                  DO KO=1,I
                     IF(TMT.lt.T(KO)) exit
                  enddo
                  

                  
                  ISUP=KO
                  INF=KO-1
                  
                  DTMT=TMT-T(INF)
                  TIS=T(ISUP)
                  TIN=T(INF)
                  DT=TIS-TIN
                  DG=GO(ISUP)-GO(INF)
                  GINF=GO(INF)+DG/DT*DTMT
                  DO L=1,NMAX
                     XTIN(L)=X(INF,L)+(X(ISUP,L)-X(INF,L))/DT*DTMT
                  enddo

                  IF(TMT.LT.TAUIN) THEN
                     br=0.
                     bro=0.
                  else
                     br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
                     bro=br
                  end if
                  SRTM=AR*T1*(1.-EXP(-TMT/T1))+BR*T2*(1.-
     $                 EXP(-(TMT-TAUIN)/T2)) +GAS
                  SMSM=ARO*T1*(1.-EXP(-TMT/T1))+BRO*T2*(1.-
     $                 EXP(-(TMT-TAUIN)/T2)) +GAS
                  GRATM=(GRAT(INF)+GRAT(ISUP))*0.5
                  RAINF=GRATM*(SRTM/SMSM)**EXPO*
     $                 (SG/SRTM)**(CK-1.)*GINF**CK
               endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
               ZETAIN=0.
               DO  L=5,NMAX
                  ZETAIN=ZETAIN+XTIN(L)
               enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               SIJ1=0.
               DO  JS=1,NMAX
                  SIJ1=SIJ1+Q(AMED,J,JS,ZETAIN)*XTIN(JS)
               enddo
               XL=TAU(PAMU(LL))
               XU=TAU(PAMU(LL+1))
               ASN1=ASN1+RAINF*FII(XL,XU)
               SOM2=SOM2+SIJ1*RAINF*FI(XL,XU)
            enddo



            PRZW2(J)=SOM2


            

            ASN2=0.
            SOM3=0.
            L6=KMAX2-1


            IF(B.ne.AMU(KMAX)) then

               
               DO LL=KB,L6
                  TDIZ=T(I)
                  ASTR=(PAMU(LL)+PAMU(LL+1))*0.5
                  TSTAR=TMEDIO-TAU(ASTR)
                  GSTAR=GO(I)+(GASP-GO(I))/DTEM*(TSTAR-TDIZ)
                  
                  DO L=1,NMAX
                     XTIN(L)=X(I,L)+(XP(L)-X(I,L))/DTEM*(TSTAR-TDIZ)
                  enddo
                  
                  IF(TMT.LT.TAUIN) THEN
                     br=0.
                     bro=0.
                  else
                     br=(vmr-gas)/T2*(1.-exp(-(ETA-TAUIN)/T2))
                     bro=br
                  end if
                  SSTAR=AR*T1*(1.-EXP(-TSTAR/T1))+BR*T2*
     $                 (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                 GAS
                  SMTAR=ARO*T1*(1.-EXP(-TSTAR/T1))+BRO*T2*
     $                 (1.-EXP(-(TSTAR-TAUIN)/T2))+
     $                 GAS
                  RAINF=GRAT(I)*(SSTAR/SMTAR)**EXPO*
     $                 (SG/SSTAR)**(CK-1.)*GSTAR**CK
                  PX=TAU(PAMU(LL))
                  PY=TAU(PAMU(LL+1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Modifica ZETAIN
                  ZETAIN=0.
                  DO  L=5,NMAX
                     ZETAIN=ZETAIN+XTIN(L)
                  enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  SIJ2=0.
                  DO JS=1,NMAX
                     SIJ2=SIJ2+Q(ASTR,J,JS,ZETAIN)*XTIN(JS)
                  enddo
                  ASN2=ASN2+RAINF*FII(PX,PY)

                  SOM3=SOM3+RAINF*SIJ2*FI(PX,PY)
               enddo
               PRZW3(J)=SOM3




            endif
            
            ASNII=ASN1+ASN2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Parte per tener conto dei bursts di nova.

            WD=RAWD3
            WDM=RAWDM3
            


            CALL COSMIC
            WI1(J)=RR1(J)+PRZW1(J)+PRZW2(J)+PRZW3(J)+RNOV(J)
     $           +COSM(J)*GO(I)*SG*X(I,1)*.8

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            WI2(J)=XINF(J)*FAP
            WI(J)=WI1(J)+WI2(J)

         endif
         


         
      enddo

 4375 FORMAT(1X,1E12.5)
      
 100  FORMAT(30X,'A E B=',2E15.5,i5)

 240  FORMAT(1X,11E12.5,/,E12.5)


      RETURN
      END

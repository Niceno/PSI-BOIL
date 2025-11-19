C
C Simple example fortran program to write a
C binary datafile for tecplot.  This example
C does the following:
C
C   1.  Open a datafile called "t.plt"
C   2.  Assign values for X,Y, and P
C   3.  Write out a zone dimensioned 4x5
C   4.  Close the datafile.
C
C
      program test

      INCLUDE 'tecio.inc'

      character*1 NULLCHR
      Integer*4   Debug,III,NPts,NElm

      Dimension X(4,5), Y(4,5), P(4,5)
      Real*8    SolTime
      Integer*4 VIsDouble
      Integer*4 ZoneType,StrandID,ParentZn,IsBlock
      Integer*4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)

      NULLCHR = CHAR(0)
      NullPtr = 0
      Debug   = 1
      VIsDouble = 0
      IMax    = 4
      JMax    = 5
      KMax    = 1
      ZoneType = 0
      SolTime = 360.0
      StrandID = 0
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0
C
C... Open the file and write the tecplot datafile 
C... header information.
C
      I = TecIni110('SIMPLE DATASET'//NULLCHR,
     &              'X Y P'//NULLCHR,
     &              't.plt'//NULLCHR,
     &              '.'//NULLCHR,
     &              Debug,
     &              VIsDouble)

      Do 10 I = 1,4
      Do 10 J = 1,5
        X(I,J) = I
        Y(I,J) = J
        P(I,J) = I*J
   10 Continue
C
C... Write the zone header information.
C
      I = TecZne110('Simple Zone'//NULLCHR,
     &              ZoneType,
     &              IMax,
     &              JMax,
     &              KMax,
     &              ICellMax,
     &              JCellMax,
     &              KCellMax,
     &              SolTime,
     &              StrandID,
     &              ParentZn,
     &              IsBlock,
     &              NFConns,
     &              FNMode,
     &              Null,
     &              Null,
     &              Null,
     &              ShrConn)
C
C... Write out the field data.
C
      III = IMax*JMax
      I   = TecDat110(III,X,0)
      I   = TecDat110(III,Y,0)
      I   = TecDat110(III,P,0)

      I = TecEnd110()
      Stop
      End

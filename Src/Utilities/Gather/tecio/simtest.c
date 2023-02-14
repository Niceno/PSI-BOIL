/*
 * Simple example c program to write a
 * binary datafile for tecplot.  This example
 * does the following:
 *
 *   1.  Open a datafile called "t.plt"
 *   2.  Assign values for X,Y, and P
 *   3.  Write out a zone dimensioned 4x5
 *   4.  Close the datafile.
 */

#include "TECIO.h"

#ifndef NULL
#define NULL 0
#endif

main ()
{
  float X[5][4], Y[5][4], P[5][4];
  double SolTime;
  INTEGER4 Debug,I,J,III,DIsDouble,VIsDouble,IMax,JMax,KMax,ZoneType,StrandID,ParentZn,IsBlock;
  INTEGER4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn;

  Debug     = 1;
  VIsDouble = 0;
  DIsDouble = 0;
  IMax      = 4;
  JMax      = 5;
  KMax      = 1;
  ZoneType  = 0;      /* Ordered */
  SolTime   = 360.0;
  StrandID  = 0;     /* StaticZone */
  ParentZn  = 0;      /* No Parent */
  IsBlock   = 1;      /* Block */
  ICellMax  = 0;
  JCellMax  = 0;
  KCellMax  = 0;
  NFConns   = 0;
  FNMode    = 0;
  ShrConn   = 0;

/*
 * Open the file and write the tecplot datafile 
 * header information 
 */
  I = TECINI110("SIMPLE DATASET",
                "X Y P",
                "t.plt",
                ".",
                &Debug,
                &VIsDouble);

  for (J = 0; J < 5; J++)
  for (I = 0; I < 4; I++)
    {
      X[J][I] = (float)(I+1);
      Y[J][I] = (float)(J+1);
      P[J][I] = (float)((I+1)*(J+1));
    }
/*
 * Write the zone header information.
 */
  I = TECZNE110("Simple Zone",
                &ZoneType,
                &IMax,
                &JMax,
                &KMax,
                &ICellMax,
                &JCellMax,
                &KCellMax,
                &SolTime,
                &StrandID,
                &ParentZn,
                &IsBlock,
                &NFConns,
                &FNMode,
                NULL,           /* PassiveVarList */
                NULL,           /* ValueLocation = Nodal */
                NULL,           /* SharVarFromZone */
                &ShrConn);
/*
 * Write out the field data.
 */
  III = IMax*JMax;
  I   = TECDAT110(&III,&X[0][0],&DIsDouble);
  I   = TECDAT110(&III,&Y[0][0],&DIsDouble);
  I   = TECDAT110(&III,&P[0][0],&DIsDouble);

  I = TECEND110();
}

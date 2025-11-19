      INTERFACE

      INTEGER(4) FUNCTION tecini
     & (Title,
     &  Variables,
     &  FName,
     &  ScratchDir,
     &  Debug,
     &  VIsDouble)
        !MS$ATTRIBUTES STDCALL :: tecini
        !MS$ATTRIBUTES REFERENCE :: Title,Variables,FName
        !MS$ATTRIBUTES REFERENCE :: ScratchDir,Debug,VIsDouble
          character(len=*) Title
          character(len=*) Variables
          character(len=*) FName
          character(len=*) ScratchDir
          INTEGER(4) Debug
          INTEGER(4) VIsDouble
      END FUNCTION tecini

      INTEGER(4) FUNCTION teczne
     & (ZoneTitle,
     &  IMx,
     &  JMx,
     &  KMx,
     &  ZFormat,
     &  DupList)
        !MS$ATTRIBUTES STDCALL :: teczne
        !MS$ATTRIBUTES REFERENCE :: ZoneTitle,IMx,JMx,KMx
        !MS$ATTRIBUTES REFERENCE :: ZFormat,DupList
        character(len=*) ZoneTitle
        INTEGER(4) IMx
        INTEGER(4) JMx
        INTEGER(4) KMx
        character(len=*) ZFormat
        character(len=*) DupList
      END FUNCTION teczne

      INTEGER(4) FUNCTION tecdat
     & (N,
     &  FieldData,
     &  IsDouble)
        !MS$ATTRIBUTES STDCALL :: tecdat
        !MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
        INTEGER(4)  N
        REAL(4)     FieldData(*)
        INTEGER(4)  IsDouble
      END FUNCTION tecdat

      INTEGER(4) FUNCTION tecnod
     & (NData)
        !MS$ATTRIBUTES STDCALL :: tecnod
        !MS$ATTRIBUTES REFERENCE :: NData
        INTEGER(4)  NData(*)
      END FUNCTION tecnod

      INTEGER(4) FUNCTION teclab
     & (S)
        !MS$ATTRIBUTES STDCALL :: teclab
        !MS$ATTRIBUTES REFERENCE :: S
        character(len=*) S
      END FUNCTION teclab

      INTEGER(4) FUNCTION tecusr
     & (S)
        !MS$ATTRIBUTES STDCALL :: tecusr
        !MS$ATTRIBUTES REFERENCE :: S
        character(len=*) S
      END FUNCTION tecusr

      INTEGER(4) FUNCTION tecend()
        !MS$ATTRIBUTES STDCALL :: tecend
      END FUNCTION tecend

      INTEGER(4) FUNCTION tecgeo
     & (XPos,
     &  YPos,
     &  ZPos,
     &  PosCoordMode,
     &  AttachToZone,
     &  Zone,
     &  Color,
     &  FillColor,
     &  IsFilled,
     &  GeomType,
     &  LinePattern,
     &  PatternLength,
     &  LineThickness,
     &  NumEllipsePts,
     &  ArrowheadStyle,
     &  ArrowheadAttachment,
     &  ArrowheadSize,
     &  ArrowheadAngle,
     &  Scope,
     &  NumSegments,
     &  NumSegPts,
     &  XGeomData,
     &  YGeomData,
     &  ZGeomData,
     &  mfc)
        !MS$ATTRIBUTES STDCALL :: tecgeo
        !MS$ATTRIBUTES REFERENCE :: XPos,YPos,ZPos,PosCoordMode
        !MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Color,FillColor
        !MS$ATTRIBUTES REFERENCE :: IsFilled,GeomType,LinePattern
        !MS$ATTRIBUTES REFERENCE :: PatternLength,LineThickness
        !MS$ATTRIBUTES REFERENCE :: NumEllipsePts,ArrowheadStyle
        !MS$ATTRIBUTES REFERENCE :: ArrowheadAttachment,ArrowheadSize
        !MS$ATTRIBUTES REFERENCE :: ArrowheadAngle,Scope,NumSegments
        !MS$ATTRIBUTES REFERENCE :: NumSegPts,XGeomData,YGeomData
        !MS$ATTRIBUTES REFERENCE :: ZGeomData,mfc
        REAL(8)        XPos
        REAL(8)        YPos
        REAL(8)        ZPos
        INTEGER(4)     PosCoordMode
        INTEGER(4)     AttachToZone
        INTEGER(4)     Zone
        INTEGER(4)     Color
        INTEGER(4)     FillColor
        INTEGER(4)     IsFilled
        INTEGER(4)     GeomType
        INTEGER(4)     LinePattern
        REAL(8)        PatternLength
        REAL(8)        LineThickness
        INTEGER(4)     NumEllipsePts
        INTEGER(4)     ArrowheadStyle
        INTEGER(4)     ArrowheadAttachment
        REAL(8)        ArrowheadSize
        REAL(8)        ArrowheadAngle
        INTEGER(4)     Scope
        INTEGER(4)     NumSegments
        INTEGER(4)     NumSegPts(*)
        REAL(4)        XGeomData(*)
        REAL(4)        YGeomData(*)
        REAL(4)        ZGeomData(*)
        character(len=*) mfc
      END FUNCTION tecgeo

      INTEGER(4) FUNCTION tectxt
     & (XPos,
     &  YPos,
     &  PosCoordMode,
     &  AttachToZone,
     &  Zone,
     &  BFont,
     &  FontHeightUnits,
     &  FontHeight,
     &  BoxType,
     &  BoxMargin,
     &  BoxLineThickness,
     &  BoxColor,
     &  BoxFillColor,
     &  Angle,
     &  Anchor,
     &  LineSpacing,
     &  TextColor,
     &  Scope,
     &  Text,
     &  mfc)
        !MS$ATTRIBUTES STDCALL :: tectxt
        !MS$ATTRIBUTES REFERENCE :: XPos,YPos,PosCoordMode,AttachToZone
        !MS$ATTRIBUTES REFERENCE :: Zone,BFont,FontHeightUnits
        !MS$ATTRIBUTES REFERENCE :: FontHeight,BoxType,BoxMargin
        !MS$ATTRIBUTES REFERENCE :: BoxLineThickness,BoxColor
        !MS$ATTRIBUTES REFERENCE :: BoxFillColor,Angle,Anchor
        !MS$ATTRIBUTES REFERENCE :: LineSpacing,TextColor,Scope
        !MS$ATTRIBUTES REFERENCE :: Text,mfc
        REAL(8)     XPos
        REAL(8)     YPos
        INTEGER(4)  PosCoordMode
        INTEGER(4)  AttachToZone
        INTEGER(4)  Zone
        INTEGER(4)  BFont
        INTEGER(4)  FontHeightUnits
        REAL(8)     FontHeight
        INTEGER(4)  BoxType
        REAL(8)     BoxMargin
        REAL(8)     BoxLineThickness
        INTEGER(4)  BoxColor
        INTEGER(4)  BoxFillColor
        REAL(8)     Angle
        INTEGER(4)  Anchor
        REAL(8)     LineSpacing
        INTEGER(4)  TextColor
        INTEGER(4)  Scope
        character(len=*) Text
        character(len=*) mfc
      END FUNCTION tectxt

      INTEGER(4) FUNCTION tecfil
     & (F)
        !MS$ATTRIBUTES STDCALL :: tecfil
        !MS$ATTRIBUTES REFERENCE :: F
        INTEGER(4)  F
      END FUNCTION tecfil

      INTEGER(4) FUNCTION tecini100
     & (Title,
     &  Variables,
     &  FName,
     &  ScratchDir,
     &  Debug,
     &  VIsDouble)
        !MS$ATTRIBUTES STDCALL :: tecini100
        !MS$ATTRIBUTES REFERENCE :: Title,Variables,FName
        !MS$ATTRIBUTES REFERENCE :: ScratchDir,Debug,VIsDouble
          character(len=*) Title
          character(len=*) Variables
          character(len=*) FName
          character(len=*) ScratchDir
          INTEGER(4) Debug
          INTEGER(4) VIsDouble
      END FUNCTION tecini100

      INTEGER FUNCTION teczne100(ZoneTitle,
     &                           ZoneType,
     &                           IMxOrNumPts,
     &                           JMxOrNumElements,
     &                           KMx,
     &                           ICellMax,
     &                           JCellMax,
     &                           KCellMax,
     &                           IsBlock,
     &                           NumFaceConnections,
     &                           FaceNeighborMode,
     &                           ValueLocation,
     &                           ShareVarFromZone,
     &                           ShareConnectivityFromZone)
        !MS$ATTRIBUTES STDCALL :: teczne100
        !MS$ATTRIBUTES REFERENCE :: ZoneTitle,ZoneType,IMxOrNumPts
        !MS$ATTRIBUTES REFERENCE :: JMxOrNumElements,KMx
        !MS$ATTRIBUTES REFERENCE :: ICellMax,JCellMax,KCellMax,IsBlock
        !MS$ATTRIBUTES REFERENCE :: NumFaceConnections,FaceNeighborMode
        !MS$ATTRIBUTES REFERENCE :: ValueLocation,ShareVarFromZone
        !MS$ATTRIBUTES REFERENCE :: ShareConnectivityFromZone
      character(len=*) ZoneTitle
      INTEGER(4) ZoneType
      INTEGER(4) IMxOrNumPts
      INTEGER(4) JMxOrNumElements
      INTEGER(4) KMx
      INTEGER(4) ICellMax
      INTEGER(4) JCellMax
      INTEGER(4) KCellMax
      INTEGER(4) IsBlock
      INTEGER(4) NumFaceConnections
      INTEGER(4) FaceNeighborMode
      INTEGER(4) ValueLocation(*)
      INTEGER(4) ShareVarFromZone(*)
      INTEGER(4) ShareConnectivityFromZone
      END FUNCTION teczne100

      INTEGER(4) FUNCTION tecdat100
     & (N,
     &  FieldData,
     &  IsDouble)
        !MS$ATTRIBUTES STDCALL :: tecdat100
        !MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
        INTEGER(4)  N
        REAL(4)     FieldData(*)
        INTEGER(4)  IsDouble
      END FUNCTION tecdat100

      INTEGER(4) FUNCTION tecnod100
     & (NData)
        !MS$ATTRIBUTES STDCALL :: tecnod100
        !MS$ATTRIBUTES REFERENCE :: NData
        INTEGER(4)  NData(*)
      END FUNCTION tecnod100

      INTEGER(4) FUNCTION teclab100
     & (S)
        !MS$ATTRIBUTES STDCALL :: teclab100
        !MS$ATTRIBUTES REFERENCE :: S
        character(len=*) S
      END FUNCTION teclab100

      INTEGER(4) FUNCTION tecusr100
     & (S)
        !MS$ATTRIBUTES STDCALL :: tecusr100
        !MS$ATTRIBUTES REFERENCE :: S
        character(len=*) S
      END FUNCTION tecusr100

      INTEGER(4) FUNCTION tecend100()
        !MS$ATTRIBUTES STDCALL :: tecend100
      END FUNCTION tecend100

      INTEGER(4) FUNCTION tecgeo100
     & (XPos,
     &  YPos,
     &  ZPos,
     &  PosCoordMode,
     &  AttachToZone,
     &  Zone,
     &  Color,
     &  FillColor,
     &  IsFilled,
     &  GeomType,
     &  LinePattern,
     &  PatternLength,
     &  LineThickness,
     &  NumEllipsePts,
     &  ArrowheadStyle,
     &  ArrowheadAttachment,
     &  ArrowheadSize,
     &  ArrowheadAngle,
     &  Scope,
     &  Clipping,
     &  NumSegments,
     &  NumSegPts,
     &  XGeomData,
     &  YGeomData,
     &  ZGeomData,
     &  mfc)
        !MS$ATTRIBUTES STDCALL :: tecgeo100
        !MS$ATTRIBUTES REFERENCE :: XPos,YPos,ZPos,PosCoordMode
        !MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Color,FillColor
        !MS$ATTRIBUTES REFERENCE :: IsFilled,GeomType,LinePattern
        !MS$ATTRIBUTES REFERENCE :: PatternLength,LineThickness
        !MS$ATTRIBUTES REFERENCE :: NumEllipsePts,ArrowheadStyle
        !MS$ATTRIBUTES REFERENCE :: ArrowheadAttachment,ArrowheadSize
        !MS$ATTRIBUTES REFERENCE :: ArrowheadAngle,Scope,Clipping
        !MS$ATTRIBUTES REFERENCE :: NumSegments,NumSegPts
        !MS$ATTRIBUTES REFERENCE :: XGeomData,YGeomData
        !MS$ATTRIBUTES REFERENCE :: ZGeomData,mfc
        REAL(8)        XPos
        REAL(8)        YPos
        REAL(8)        ZPos
        INTEGER(4)     PosCoordMode
        INTEGER(4)     AttachToZone
        INTEGER(4)     Zone
        INTEGER(4)     Color
        INTEGER(4)     FillColor
        INTEGER(4)     IsFilled
        INTEGER(4)     GeomType
        INTEGER(4)     LinePattern
        REAL(8)        PatternLength
        REAL(8)        LineThickness
        INTEGER(4)     NumEllipsePts
        INTEGER(4)     ArrowheadStyle
        INTEGER(4)     ArrowheadAttachment
        REAL(8)        ArrowheadSize
        REAL(8)        ArrowheadAngle
        INTEGER(4)     Scope
        INTEGER(4)     Clipping
        INTEGER(4)     NumSegments
        INTEGER(4)     NumSegPts(*)
        REAL(4)        XGeomData(*)
        REAL(4)        YGeomData(*)
        REAL(4)        ZGeomData(*)
        character(len=*) mfc
      END FUNCTION tecgeo100

      INTEGER(4) FUNCTION tectxt100
     & (XOrThetaPos,
     &  YOrRPos,
     &  ZOrUnusedPos,
     &  PosCoordMode,
     &  AttachToZone,
     &  Zone,
     &  Font,
     &  FontHeightUnits,
     &  FontHeight,
     &  BoxType,
     &  BoxMargin,
     &  BoxLineThickness,
     &  BoxColor,
     &  BoxFillColor,
     &  Angle,
     &  Anchor,
     &  LineSpacing,
     &  TextColor,
     &  Scope,
     &  Clipping,
     &  Text,
     &  mfc)
        !MS$ATTRIBUTES STDCALL :: tectxt100
        !MS$ATTRIBUTES REFERENCE :: XOrThetaPos,YOrRPos
        !MS$ATTRIBUTES REFERENCE :: ZOrUnusedPos,PosCoordMode
        !MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Font
        !MS$ATTRIBUTES REFERENCE :: FontHeightUnits
        !MS$ATTRIBUTES REFERENCE :: FontHeight,BoxType,BoxMargin
        !MS$ATTRIBUTES REFERENCE :: BoxLineThickness,BoxColor
        !MS$ATTRIBUTES REFERENCE :: BoxFillColor,Angle,Anchor
        !MS$ATTRIBUTES REFERENCE :: LineSpacing,TextColor,Scope,Clipping
        !MS$ATTRIBUTES REFERENCE :: Text,mfc
        REAL(8)     XOrThetaPos
        REAL(8)     YOrRPos
        REAL(8)     ZOrUnusedPos
        INTEGER(4)  PosCoordMode
        INTEGER(4)  AttachToZone
        INTEGER(4)  Zone
        INTEGER(4)  Font
        INTEGER(4)  FontHeightUnits
        REAL(8)     FontHeight
        INTEGER(4)  BoxType
        REAL(8)     BoxMargin
        REAL(8)     BoxLineThickness
        INTEGER(4)  BoxColor
        INTEGER(4)  BoxFillColor
        REAL(8)     Angle
        INTEGER(4)  Anchor
        REAL(8)     LineSpacing
        INTEGER(4)  TextColor
        INTEGER(4)  Scope
        INTEGER(4)  Clipping
        character(len=*) Text
        character(len=*) mfc
      END FUNCTION tectxt100

      INTEGER(4) FUNCTION tecfil100
     & (F)
        !MS$ATTRIBUTES STDCALL :: tecfil100
        !MS$ATTRIBUTES REFERENCE :: F
        INTEGER(4)  F
      END FUNCTION tecfil100

      INTEGER FUNCTION tecauxstr100(Name,
     &                            Value)
        !MS$ATTRIBUTES STDCALL :: tecauxstr100
        !MS$ATTRIBUTES REFERENCE :: Name, Value
        character(len=*) Name 
        character(len=*) Value
      END FUNCTION tecauxstr100

      INTEGER FUNCTION teczauxstr100(Name,
     &                             Value)
        !MS$ATTRIBUTES STDCALL :: teczauxstr100
        !MS$ATTRIBUTES REFERENCE :: Name, Value
        character(len=*) Name 
        character(len=*) Value
      END FUNCTION teczauxstr100

      INTEGER(4) FUNCTION tecface100(FaceConnections)
        !MS$ATTRIBUTES STDCALL :: tecface100
        !MS$ATTRIBUTES REFERENCE :: FaceConnections
        INTEGER(4) FACECONNECTIONS
      END FUNCTION tecface100


      INTEGER(4) FUNCTION tecini110
     & (Title,
     &  Variables,
     &  FName,
     &  ScratchDir,
     &  Debug,
     &  VIsDouble)
        !MS$ATTRIBUTES STDCALL :: tecini110
        !MS$ATTRIBUTES REFERENCE :: Title,Variables,FName
        !MS$ATTRIBUTES REFERENCE :: ScratchDir,Debug,VIsDouble
          character(len=*) Title
          character(len=*) Variables
          character(len=*) FName
          character(len=*) ScratchDir
          INTEGER(4) Debug
          INTEGER(4) VIsDouble
      END FUNCTION tecini110


      INTEGER(4) FUNCTION TECZNE110(ZoneTitle,
     &                             ZoneType,
     &                             IMxOrNumPts,
     &                             JMxOrNumElements,
     &                             KMx,
     &                             ICellMax,
     &                             JCellMax,
     &                             KCellMax,
     &                             SolutionTime,
     &                             StrandID,
     &                             ParentZone,
     &                             IsBlock,
     &                             NumFaceConnections,
     &                             FaceNeighborMode,
     &                             PassiveVarList,
     &                             ValueLocation,
     &                             ShareVarFromZone,
     &                             ShareConnectivityFromZone)
        !MS$ATTRIBUTES STDCALL :: teczne110
        !MS$ATTRIBUTES REFERENCE :: ZoneTitle,ZoneType,IMxOrNumPts
        !MS$ATTRIBUTES REFERENCE :: JMxOrNumElements,KMx
        !MS$ATTRIBUTES REFERENCE :: ICellMax,JCellMax,KCellMax
        !MS$ATTRIBUTES REFERENCE :: SolutionTime,StrandID,ParentZone
        !MS$ATTRIBUTES REFERENCE :: IsBlock,PassiveVarList
        !MS$ATTRIBUTES REFERENCE :: NumFaceConnections,FaceNeighborMode
        !MS$ATTRIBUTES REFERENCE :: ValueLocation,ShareVarFromZone
        !MS$ATTRIBUTES REFERENCE :: ShareConnectivityFromZone
      character(len=*) ZoneTitle
      INTEGER(4) ZoneType
      INTEGER(4) IMxOrNumPts
      INTEGER(4) JMxOrNumElements
      INTEGER(4) KMx
      INTEGER(4) ICellMax
      INTEGER(4) JCellMax
      INTEGER(4) KCellMax
      DOUBLE PRECISION SolutionTime
      INTEGER(4) StrandID
      INTEGER(4) ParentZone
      INTEGER(4) IsBlock
      INTEGER(4) NumFaceConnections
      INTEGER(4) FaceNeighborMode
      INTEGER(4) PassiveVarList(*)
      INTEGER(4) ValueLocation(*)
      INTEGER(4) ShareVarFromZone(*)
      INTEGER(4) ShareConnectivityFromZone
      END FUNCTION teczne110

      INTEGER(4) FUNCTION tecdat110
     & (N,
     &  FieldData,
     &  IsDouble)
        !MS$ATTRIBUTES STDCALL :: tecdat110
        !MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
        INTEGER(4)  N
        REAL(4)     FieldData(*)
        INTEGER(4)  IsDouble
      END FUNCTION tecdat110

      INTEGER(4) FUNCTION tecnod110
     & (NData)
        !MS$ATTRIBUTES STDCALL :: tecnod110
        !MS$ATTRIBUTES REFERENCE :: NData
        INTEGER(4)  NData(*)
      END FUNCTION tecnod110

      INTEGER(4) FUNCTION teclab110
     & (S)
        !MS$ATTRIBUTES STDCALL :: teclab110
        !MS$ATTRIBUTES REFERENCE :: S
        character(len=*) S
      END FUNCTION teclab110

      INTEGER(4) FUNCTION tecusr110
     & (S)
        !MS$ATTRIBUTES STDCALL :: tecusr110
        !MS$ATTRIBUTES REFERENCE :: S
        character(len=*) S
      END FUNCTION tecusr110

      INTEGER(4) FUNCTION tecend110()
        !MS$ATTRIBUTES STDCALL :: tecend110
      END FUNCTION tecend110

      INTEGER(4) FUNCTION tecgeo110
     & (XPos,
     &  YPos,
     &  ZPos,
     &  PosCoordMode,
     &  AttachToZone,
     &  Zone,
     &  Color,
     &  FillColor,
     &  IsFilled,
     &  GeomType,
     &  LinePattern,
     &  PatternLength,
     &  LineThickness,
     &  NumEllipsePts,
     &  ArrowheadStyle,
     &  ArrowheadAttachment,
     &  ArrowheadSize,
     &  ArrowheadAngle,
     &  Scope,
     &  Clipping,
     &  NumSegments,
     &  NumSegPts,
     &  XGeomData,
     &  YGeomData,
     &  ZGeomData,
     &  mfc)
        !MS$ATTRIBUTES STDCALL :: tecgeo110
        !MS$ATTRIBUTES REFERENCE :: XPos,YPos,ZPos,PosCoordMode
        !MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Color,FillColor
        !MS$ATTRIBUTES REFERENCE :: IsFilled,GeomType,LinePattern
        !MS$ATTRIBUTES REFERENCE :: PatternLength,LineThickness
        !MS$ATTRIBUTES REFERENCE :: NumEllipsePts,ArrowheadStyle
        !MS$ATTRIBUTES REFERENCE :: ArrowheadAttachment,ArrowheadSize
        !MS$ATTRIBUTES REFERENCE :: ArrowheadAngle,Scope,Clipping
        !MS$ATTRIBUTES REFERENCE :: NumSegments,NumSegPts
        !MS$ATTRIBUTES REFERENCE :: XGeomData,YGeomData
        !MS$ATTRIBUTES REFERENCE :: ZGeomData,mfc
        REAL(8)        XPos
        REAL(8)        YPos
        REAL(8)        ZPos
        INTEGER(4)     PosCoordMode
        INTEGER(4)     AttachToZone
        INTEGER(4)     Zone
        INTEGER(4)     Color
        INTEGER(4)     FillColor
        INTEGER(4)     IsFilled
        INTEGER(4)     GeomType
        INTEGER(4)     LinePattern
        REAL(8)        PatternLength
        REAL(8)        LineThickness
        INTEGER(4)     NumEllipsePts
        INTEGER(4)     ArrowheadStyle
        INTEGER(4)     ArrowheadAttachment
        REAL(8)        ArrowheadSize
        REAL(8)        ArrowheadAngle
        INTEGER(4)     Scope
        INTEGER(4)     Clipping
        INTEGER(4)     NumSegments
        INTEGER(4)     NumSegPts(*)
        REAL(4)        XGeomData(*)
        REAL(4)        YGeomData(*)
        REAL(4)        ZGeomData(*)
        character(len=*) mfc
      END FUNCTION tecgeo110

      INTEGER(4) FUNCTION tectxt110
     & (XOrThetaPos,
     &  YOrRPos,
     &  ZOrUnusedPos,
     &  PosCoordMode,
     &  AttachToZone,
     &  Zone,
     &  Font,
     &  FontHeightUnits,
     &  FontHeight,
     &  BoxType,
     &  BoxMargin,
     &  BoxLineThickness,
     &  BoxColor,
     &  BoxFillColor,
     &  Angle,
     &  Anchor,
     &  LineSpacing,
     &  TextColor,
     &  Scope,
     &  Clipping,
     &  Text,
     &  mfc)
        !MS$ATTRIBUTES STDCALL :: tectxt110
        !MS$ATTRIBUTES REFERENCE :: XOrThetaPos,YOrRPos
        !MS$ATTRIBUTES REFERENCE :: ZOrUnusedPos,PosCoordMode
        !MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Font
        !MS$ATTRIBUTES REFERENCE :: FontHeightUnits
        !MS$ATTRIBUTES REFERENCE :: FontHeight,BoxType,BoxMargin
        !MS$ATTRIBUTES REFERENCE :: BoxLineThickness,BoxColor
        !MS$ATTRIBUTES REFERENCE :: BoxFillColor,Angle,Anchor
        !MS$ATTRIBUTES REFERENCE :: LineSpacing,TextColor,Scope,Clipping
        !MS$ATTRIBUTES REFERENCE :: Text,mfc
        REAL(8)     XOrThetaPos
        REAL(8)     YOrRPos
        REAL(8)     ZOrUnusedPos
        INTEGER(4)  PosCoordMode
        INTEGER(4)  AttachToZone
        INTEGER(4)  Zone
        INTEGER(4)  Font
        INTEGER(4)  FontHeightUnits
        REAL(8)     FontHeight
        INTEGER(4)  BoxType
        REAL(8)     BoxMargin
        REAL(8)     BoxLineThickness
        INTEGER(4)  BoxColor
        INTEGER(4)  BoxFillColor
        REAL(8)     Angle
        INTEGER(4)  Anchor
        REAL(8)     LineSpacing
        INTEGER(4)  TextColor
        INTEGER(4)  Scope
        INTEGER(4)  Clipping
        character(len=*) Text
        character(len=*) mfc
      END FUNCTION tectxt110

      INTEGER(4) FUNCTION tecfil110
     & (F)
        !MS$ATTRIBUTES STDCALL :: tecfil110
        !MS$ATTRIBUTES REFERENCE :: F
        INTEGER(4)  F
      END FUNCTION tecfil110

      INTEGER FUNCTION tecauxstr110(Name,
     &                            Value)
        !MS$ATTRIBUTES STDCALL :: tecauxstr110
        !MS$ATTRIBUTES REFERENCE :: Name, Value
        character(len=*) Name 
        character(len=*) Value
      END FUNCTION tecauxstr110

      INTEGER FUNCTION teczauxstr110(Name,
     &                             Value)
        !MS$ATTRIBUTES STDCALL :: teczauxstr110
        !MS$ATTRIBUTES REFERENCE :: Name, Value
        character(len=*) Name 
        character(len=*) Value
      END FUNCTION teczauxstr110

      INTEGER(4) FUNCTION tecface110(FaceConnections)
        !MS$ATTRIBUTES STDCALL :: tecface110
        !MS$ATTRIBUTES REFERENCE :: FaceConnections
        INTEGER(4) FACECONNECTIONS(*)
      END FUNCTION tecface110



      END INTERFACE

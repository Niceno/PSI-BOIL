 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <map>
#include <string>
#include <sstream>
#include <stdarg.h>
#include <boost/assign.hpp>
#include "ThirdPartyHeadersEnd.h"
namespace tecplot { namespace ___3931 { typedef std::map<uint16_t, std::string> ___3940; typedef std::map<uint16_t, uint64_t> ___3943; typedef std::map<std::string, std::string> TagDescriptionToDescriptionMap; std::string getFileNameSansFolder(std::string const& ___1392); }}
 #define SZPLT_FILE_TYPE_NAME          "SZPLT files"
 #define SZPLT_FILE_FILTER             "*.szplt"
 #define SZPLT_FILE_TYPE               "SZPLT"
 #define SZPLT_MIN_READ_VERSION        104
 #define SZPLT_MAX_READ_VERSION        105
 #define SZPLT_CUR_WRITE_VERSION       105
 #define SZPLT_FILE_MARKER             0x01020304
 #define SZPLT_ZONE_MARKER             0x02030405
 #define SZPLT_CSZ_CONNECTIVITY_MARKER 0x03040506
 #define SZPLT_NSZ_CONNECTIVITY_MARKER 0x04050607
 #define SZPLT_CSZ_FIELD_DATA_MARKER   0x05060708
 #define SZPLT_NSZ_FIELD_DATA_MARKER   0x06070809
 #define SZPLT_CC_FIELD_DATA_MARKER    0x07080910
 #define SZPLT_NODAL_FIELD_DATA_MARKER 0x08091011
 #define SZPLT_NODE_MAP_MARKER         0x09101112
 #define SZPLT_USER_FACE_NEIGHBOR_MARKER 0x10111213
 #define SZPLT_AUX_DATA_MARKER         0x11121314
 #define SZPLT_GEOMETRY_MARKER         0x12131415
 #define SZPLT_TEXT_MARKER             0x13141516
 #define SZPLT_CUSTOM_LABELS_MARKER    0x14151617
 #define SZPLT_PARTITION_MARKER        0x15161718
 #define SZPLT_ZONE_MARKER_DESCRIPTION               "zoneMarker*"
 #define SZPLT_ZONE_NUM_DESCRIPTION                  "zoneNum*"
 #define SZPLT_VARIABLE_NUM_DESCRIPTION              "variableNum*"
 #define SZPLT_CC_FIELD_DATA_MARKER_DESCRIPTION      "ccFieldDataMarker*"
 #define SZPLT_NODAL_FIELD_DATA_MARKER_DESCRIPTION   "nodalFieldDataMarker*"
 #define SZPLT_CSZ_FIELD_DATA_MARKER_DESCRIPTION     "cszFieldDataMarker*"
 #define SZPLT_NSZ_FIELD_DATA_MARKER_DESCRIPTION     "nszFieldDataMarker*"
 #define SZPLT_SUBZONE_NUM_DESCRIPTION               "subzoneNum*"
 #define SZL_ISTART_DESCRIPTION                      "iStart*"
 #define SZL_IEND_DESCRIPTION                        "iEnd*"
 #define SZL_JSTART_DESCRIPTION                      "jStart*"
 #define SZL_JEND_DESCRIPTION                        "jEnd*"
 #define SZL_KSTART_DESCRIPTION                      "kStart*"
 #define SZL_KEND_DESCRIPTION                        "kEnd*"
 #define SZPLT_CSZ_FIELD_DATA_DESCRIPTION            "cszFieldData"
 #define SZPLT_NSZ_FIELD_DATA_DESCRIPTION            "nszFieldData"
 #define SZPLT_NODE_MAP_MARKER_DESCRIPTION           "nodeMapMarker*"
 #define SZPLT_NSZ_CONNECTIVITY_MARKER_DESCRIPTION   "nszConnectivityMarker*"
 #define SZPLT_CSZ_CONNECTIVITY_MARKER_DESCRIPTION   "cszConnectivityMarker*"
 #define SZPLT_USER_FACE_NEIGHBORS_COUNT_DESCRIPTION "userFaceNeighborCount"
 #define SZPLT_USER_FACE_NEIGHBOR_MARKER_DESCRIPTION "userFaceNeighborMarker*"
 #define SZPLT_USER_FACE_NEIGHBORS_DESCRIPTION       "UserFaceNeighbors"
 #define SZPLT_FIELD_DATA_DESCRIPTION                "FieldData"
 #define SZPLT_NODEMAP_DESCRIPTION                   "NodeMap"
 #define SZPLT_AUX_DATA_MARKER_DESCRIPTION "auxData*"
 #define SZPLT_GEOMETRY_MARKER_DESCRIPTION "geometry*"
 #define SZPLT_TEXT_MARKER_DESCRIPTION "text*"
 #define SZPLT_CUSTOM_LABELS_MARKER_DESCRIPTION "customLabels*"
 #define SZPLT_PARTITION_MARKER_DESCRIPTION     "partitionMarker*"
 #define SZPLT_PARTITION_NUM_DESCRIPTION        "partitionNum*"
 #define SZPLT_NUM_INTERMEDIATE_HEADERS "numIntermediateHeaders"
 #define SZPLT_FILE_HEADER_MIN_SIZE 16
 #define ___206 40
 #define VALID_DESCRIPTION_OR_NULL(___970) (((___970)==NULL) || (VALID_REF((___970)) && strlen((___970))<___206))
 #define VALID_DESCRIPTION(___970) (VALID_REF((___970)) && strlen((___970)) < ___206)
 #define ___194      1  
 #define ASCII_INT16_LEN     6  
 #define ___208    5  
 #define ASCII_INT32_LEN     11 
 #define ___209    10 
 #define ASCII_INT64_LEN     20 
 #define ___210    19 
 #define ___211     3  
 #define ___196     15 
 #define ___195    23 
 #define ASCII_HEXUINT8_LEN  2  
 #define ASCII_HEXUINT16_LEN 4  
 #define ASCII_HEXUINT32_LEN 8  
 #define ___197 16 
 #define ASCII_SPACING_LEN   2  
 #define ASCII_VAR_LEN       16 
 #define ASCII_UINT8_MINMAX_LEN  (2*___211 + 1)
 #define ASCII_INT16_MINMAX_LEN  (2*ASCII_INT16_LEN + 1)
 #define ASCII_INT32_MINMAX_LEN  (2*ASCII_INT32_LEN + 1)
 #define ASCII_FLOAT_MINMAX_LEN  (2*___196 + 1)
 #define ASCII_DOUBLE_MINMAX_LEN (2*___195 + 1)
 #define ___1383 0x0102030405060708
 #define ___1384 "fileByteOrderIndicator"
 #define ___893 "dataSetHeaderFileLoc"
static uint16_t const  ___2847 = 0x0001; static uint16_t const  ___2843 = 0x0002; static uint16_t const  SUBZONE_MAX_I_SIZE_TAG = 0x0003; static uint16_t const  SUBZONE_MAX_J_SIZE_TAG = 0x0004; static uint16_t const  SUBZONE_MAX_K_SIZE_TAG = 0x0005; static uint16_t const  SUBZONE_MAX_FE_SIZE_TAG = 0x0006; static uint16_t const  ___241 = 0x0007; static uint16_t const  ___1579 = 0x0008; static uint16_t const  ___793 = 0x0009; static uint16_t const  ___4081 = 0x000A; static uint16_t const  NUM_PARTITION_FILES_TAG = 0x000B; static uint16_t const  DATA_SET_TITLE_TAG = 0X000C; static uint16_t const  ___4649 = 0x0100; static uint16_t const  ZONE_TYPE_FILE_LOC_TAG = 0x0101; static uint16_t const  ___4628 = 0x0102; static uint16_t const  ___4608 = 0x0103; static uint16_t const  ___4674 = 0x0104; static uint16_t const  ___4684 = 0x0105; static uint16_t const  ___4668 = 0x0106; static uint16_t const  ZONE_IMAX_OR_NUM_NODES_FILE_LOC_TAG = 0x0107; static uint16_t const  ZONE_JMAX_OR_NUM_CELLS_FILE_LOC_TAG = 0x0108; static uint16_t const  ZONE_KMAX_OR_NUM_CORNERS_FILE_LOC_TAG = 0x0109; static uint16_t const  ___4672 = 0x010C; static uint16_t const  ___4363 = 0x0200; static uint16_t const  VAR_ZONE_MIN_MAX_FILE_LOC_TAG = 0x0300; static uint16_t const  ___4391 = 0x0301; static uint16_t const  ___4382 = 0x0302; static uint16_t const  ___4385 = 0x0303; static uint16_t const  ___4388 = 0x0304;
 #define ___2837 "numTags"
 #define ___3941 "tags"
 #define DATA_SET_TITLE_LENGTH "dataSetTitleLength"
 #define DATA_SET_TITLE_DESCRIPTION "dataSetTitle"
 #define NUM_ZONES_DESCRIPTION "numberOfZones"
 #define NUM_VARS_DESCRIPTION "numberOfVars"
 #define SUBZONE_MAX_I_SIZE_DESCRIPTION "subzoneMaxISize"
 #define SUBZONE_MAX_J_SIZE_DESCRIPTION "subzoneMaxJSize"
 #define SUBZONE_MAX_K_SIZE_DESCRIPTION "subzoneMaxKSize"
 #define SUBZONE_MAX_FE_SIZE_DESCRIPTION "subzoneMaxFESize"
 #define AUX_DATA_FILE_LOC_DESCRIPTION "auxDataFileLoc"
 #define GEOMETRY_FILE_LOC_DESCRIPTION "geometry"
 #define CUSTOM_LABELS_FILE_LOC_DESCRIPTION "customLabels"
 #define TEXT_FILE_LOC_DESCRIPTION "text"
 #define NUM_PARTITION_FILES_DESCRIPTION "numPartitionFiles"
 #define ZONE_NAME_FILE_LOC_DESCRIPTION "zoneNames"
 #define ZONE_TYPE_FILE_LOC_DESCRIPTION "zoneTypes"
 #define ZONE_HEADER_FILE_LOC_DESCRIPTION "zoneHeaderFileLocs"
 #define ZONE_FACE_NEIGHBOR_MODE_FILE_LOC_DESCRIPTION "zoneFaceNbrModes"
 #define ZONE_SOLUTION_TIME_FILE_LOC_DESCRIPTION "zoneSolutionTimes"
 #define ZONE_STRAND_ID_FILE_LOC_DESCRIPTION "zoneStrandIDs"
 #define ZONE_PARENT_ZONE_FILE_LOC_DESCRIPTION "zoneParentZones"
 #define ZONE_IMAX_OR_NUM_NODES_FILE_LOC_DESCRIPTION "zoneIMaxOrNumNodes"
 #define ZONE_JMAX_OR_NUM_CELLS_FILE_LOC_DESCRIPTION "zoneJMaxOrNumCells"
 #define ZONE_KMAX_OR_NUM_CORNERS_FILE_LOC_DESCRIPTION "zoneKMaxOrNumFaces"
 #define ZONE_SHARE_CONNECTIVITY_FILE_LOC_DESCRIPTION "zoneShareConnectivity"
 #define VAR_NAME_FILE_LOC_DESCRIPTION "varNames"
 #define VAR_ZONE_MIN_MAX_FILE_LOC_DESCRIPTION "zoneMinMaxes"
 #define VAR_ZONE_VALUE_LOCATION_FILE_LOC_DESCRIPTION "valueLocations"
 #define VAR_ZONE_FIELD_DATA_TYPE_FILE_LOC_DESCRIPTION "fieldDataTypes"
 #define VAR_ZONE_IS_PASSIVE_FILE_LOC_DESCRIPTION "passiveVar"
 #define VAR_ZONE_SHARE_FILE_LOC_DESCRIPTION "shareVar"
 #define ___233 "numAuxData"
 #define NUM_GEOMS_DESCRIPTION "numGeometries"
 #define ___791 "numCustomLabelSets"
 #define ___2818 "numLabels"
 #define ___792 "customLabelSet"
 #define NUM_TEXTS_DESCRIPTION "numTexts"
 #define ___4648 "zoneNamesArray"
 #define ___4692 "zoneTypesArray"
 #define ___4626 "zoneHeaderFileLocsArray"
 #define ___4607 "zoneFaceNbrModesArray"
 #define ___4673 "zoneSolutionTimesArray"
 #define ___4683 "zoneStrandIDsArray"
 #define ___4667 "zoneParentZonesArray"
 #define ZONE_IMAX_OR_NUM_NODES_ARRAY_DESCRIPTION "zoneIMaxOrNumNodesArray"
 #define ZONE_JMAX_OR_NUM_CELLS_ARRAY_DESCRIPTION "zoneJMaxOrNumCellsArray"
 #define ZONE_KMAX_OR_NUM_CORNERS_ARRAY_DESCRIPTION "zoneKMaxOrNumFacesArray"
 #define ___4671 "zoneShareConnectivityArray"
 #define ___4362 "varNamesArray"
 #define VAR_ZONE_MIN_MAX_ARRAY_DESCRIPTION "zoneMinMaxesForVar"
 #define ___4390 "valueLocationsForVar"
 #define ___4381 "fieldDataTypesForVar"
 #define ___4384 "passiveForVar"
 #define ___4387 "shareForVar"
 #define ___252 "auxDataLocation"
 #define ___239 "auxDataEntity"
 #define ___262 "auxDataNameLength"
 #define ___261 "auxDataName"
 #define ___268 "auxDataType"
 #define ___272 "auxDataValueLength"
 #define ___271 "auxDataValue"
 #define ___1651 "xOrThetaPos"
 #define ___1653 "yOrRPos"
 #define ___1656 "zPos"
 #define ___1625 "posCoordMode"
 #define ___1609 "isAttached"
 #define ___1655 "zone"
 #define ___1569 "color"
 #define ___1580 "fillColor"
 #define ___1610 "isFilled"
 #define GEOM_GEOM_TYPE "geomType"
 #define ___1611 "linePattern"
 #define ___1623 "patternLength"
 #define ___1612 "lineThickness"
 #define ___1619 "numEllipsePts"
 #define ___1567 "arrowheadStyle"
 #define ___1561 "arrowheadAttachment"
 #define ___1566 "arrowheadSize"
 #define ___1560 "arrowheadAngle"
 #define ___1631 "scope"
 #define GEOM_CLIPPING "clipping"
 #define ___1614 "macroFunctionCommandLength"
 #define ___1613 "macroFunctionCommand"
 #define ___1628 "rectangleWidth"
 #define ___1627 "rectangleHeight"
 #define ___1647 "squareSize"
 #define ___1606 "circleRadius"
 #define ___1576 "ellipseHorizontalAxis"
 #define ___1577 "ellipseVerticalAxis"
 #define ___1621 "numSegments"
 #define ___1620 "numPointsInSegment"
 #define ___1652 "xValuesForSegment"
 #define ___1654 "yValuesForSegment"
 #define ___1657 "zValuesForSegment"
 #define ___4134 "xOrThetaPos"
 #define ___4135 "yOrRPos"
 #define ___4138 "zPos"
 #define ___4080 "textColor"
 #define ___4058 "textBoxColor"
 #define ___4060 "textBoxFillColor"
 #define ___4070 "textBoxLineThickness"
 #define ___4072 "textBoxMargin"
 #define TEXT_BOX_TYPE "textBoxType"
 #define ___4042 "textAnchor"
 #define ___4055 "textAngle"
 #define TEXT_CLIPPING "textClipping"
 #define TEXT_HEIGHT "textHeight"
 #define ___4106 "textLineSpacing"
 #define ___4109 "textMacroFunctionCmdLength"
 #define ___4108 "textMacroFunctionCmd"
 #define ___4114 "textPositionCoordSys"
 #define ___4118 "textScope"
 #define ___4123 "textSizeUnits"
 #define ___4126 "textStringLength"
 #define ___4125 "textString"
 #define ___4129 "textTypefaceFamilyLength"
 #define ___4128 "textTypefaceFamily"
 #define ___4131 "textTypefaceIsBold"
 #define ___4133 "textTypefaceIsItalic"
 #define ___4137 "textZoneOrMap"
 #define ___4104 "textIsAttached"
 #define DATASET_HEADER_TAG_MAP_INITIALIZER boost::assign::map_list_of<uint16_t, uint64_t> \
 (DATA_SET_TITLE_TAG, 0) \
 (___2847, 0) \
 (___2843, 0) \
 (SUBZONE_MAX_I_SIZE_TAG, 0) \
 (SUBZONE_MAX_J_SIZE_TAG, 0) \
 (SUBZONE_MAX_K_SIZE_TAG, 0) \
 (SUBZONE_MAX_FE_SIZE_TAG, 0) \
 (___241, 0) \
 (___1579, 0) \
 (___793, 0) \
 (___4081, 0) \
 (___4649, 0) \
 (ZONE_TYPE_FILE_LOC_TAG, 0) \
 (___4628, 0) \
 (___4608, 0) \
 (___4674, 0) \
 (___4684, 0) \
 (___4668, 0) \
 (ZONE_IMAX_OR_NUM_NODES_FILE_LOC_TAG, 0) \
 (ZONE_JMAX_OR_NUM_CELLS_FILE_LOC_TAG, 0) \
 (ZONE_KMAX_OR_NUM_CORNERS_FILE_LOC_TAG, 0) \
 (___4672, 0) \
 (___4363, 0) \
 (VAR_ZONE_MIN_MAX_FILE_LOC_TAG, 0) \
 (___4391, 0) \
 (___4382, 0) \
 (___4385, 0) \
 (___4388, 0)
 #define DATSET_HEADER_TAG_DESCRIPTIONS_MAP_INITIALIZER boost::assign::map_list_of<uint16_t, std::string> \
 (DATA_SET_TITLE_TAG, DATA_SET_TITLE_DESCRIPTION) \
 (___2847, NUM_ZONES_DESCRIPTION) \
 (___2843, NUM_VARS_DESCRIPTION) \
 (SUBZONE_MAX_I_SIZE_TAG, SUBZONE_MAX_I_SIZE_DESCRIPTION) \
 (SUBZONE_MAX_J_SIZE_TAG, SUBZONE_MAX_J_SIZE_DESCRIPTION) \
 (SUBZONE_MAX_K_SIZE_TAG, SUBZONE_MAX_K_SIZE_DESCRIPTION) \
 (SUBZONE_MAX_FE_SIZE_TAG, SUBZONE_MAX_FE_SIZE_DESCRIPTION) \
 (___241, AUX_DATA_FILE_LOC_DESCRIPTION) \
 (___1579, GEOMETRY_FILE_LOC_DESCRIPTION) \
 (___793, CUSTOM_LABELS_FILE_LOC_DESCRIPTION) \
 (___4081, TEXT_FILE_LOC_DESCRIPTION) \
 (NUM_PARTITION_FILES_TAG, NUM_PARTITION_FILES_DESCRIPTION) \
 (___4649, ZONE_NAME_FILE_LOC_DESCRIPTION) \
 (ZONE_TYPE_FILE_LOC_TAG, ZONE_TYPE_FILE_LOC_DESCRIPTION) \
 (___4628, ZONE_HEADER_FILE_LOC_DESCRIPTION) \
 (___4608, ZONE_FACE_NEIGHBOR_MODE_FILE_LOC_DESCRIPTION) \
 (___4674, ZONE_SOLUTION_TIME_FILE_LOC_DESCRIPTION) \
 (___4684, ZONE_STRAND_ID_FILE_LOC_DESCRIPTION) \
 (___4668, ZONE_PARENT_ZONE_FILE_LOC_DESCRIPTION) \
 (ZONE_IMAX_OR_NUM_NODES_FILE_LOC_TAG, ZONE_IMAX_OR_NUM_NODES_FILE_LOC_DESCRIPTION) \
 (ZONE_JMAX_OR_NUM_CELLS_FILE_LOC_TAG, ZONE_JMAX_OR_NUM_CELLS_FILE_LOC_DESCRIPTION) \
 (ZONE_KMAX_OR_NUM_CORNERS_FILE_LOC_TAG, ZONE_KMAX_OR_NUM_CORNERS_FILE_LOC_DESCRIPTION) \
 (___4672, ZONE_SHARE_CONNECTIVITY_FILE_LOC_DESCRIPTION) \
 (___4363, VAR_NAME_FILE_LOC_DESCRIPTION) \
 (VAR_ZONE_MIN_MAX_FILE_LOC_TAG, VAR_ZONE_MIN_MAX_FILE_LOC_DESCRIPTION) \
 (___4391, VAR_ZONE_VALUE_LOCATION_FILE_LOC_DESCRIPTION) \
 (___4382, VAR_ZONE_FIELD_DATA_TYPE_FILE_LOC_DESCRIPTION) \
 (___4385, VAR_ZONE_IS_PASSIVE_FILE_LOC_DESCRIPTION) \
 (___4388, VAR_ZONE_SHARE_FILE_LOC_DESCRIPTION)
static uint16_t const ___2782 = 0x0400; static uint16_t const ___2822 = 0x0401;
 #define NUM_CELL_SUBZONES_DESCRIPTION "numCellSubzones"
 #define NUM_NODE_SUBZONES_DESCRIPTION "numNodeSubzones"
static uint16_t const NSZ_CONNECT_FILE_LOC_TAG = 0x0500; static uint16_t const CSZ_CONNECT_FILE_LOC_TAG = 0x0501; static uint16_t const UNUSED1_ZONE_TAG = 0x0502; static uint16_t const NUM_REF_PARTITIONS_TAG = 0x0503; static uint16_t const REF_PARTITIONS_TAG = 0x0504; static uint16_t const NUM_REF_NODE_SUBZONES_TAG = 0x0505; static uint16_t const NUM_REF_CELL_SUBZONES_TAG = 0x0506; static uint16_t const CELL_SUBZONE_INCLUDES_PARTITIONS_TAG = 0x0507; static uint16_t const NODE_SUBZONE_INCLUDES_PARTITIONS_TAG = 0x0508;
 #define NSZ_CONNECT_FILE_LOC_DESCRIPTION "nszConnFileLocs"
 #define CSZ_CONNECT_FILE_LOC_DESCRIPTION "cszConnFileLocs"
 #define ZONE_UNUSED1_TAG_DESCRIPTION "zoneUnusedTag1"
 #define NUM_REF_PARTITIONS_DESCRIPTION "numRefPartitions"
 #define REF_PARTITIONS_DESCRIPTION "refPartitions"
 #define NUM_REF_NODE_SUBZONES_DESCRIPTION "cszNumRefNszs"
 #define NUM_REF_CELL_SUBZONES_DESCRIPTION "nszNumRefCszs"
 #define NUM_REF_NODE_SUBZONES_ARRAY_DESCRIPTION "cszNumRefNszsForZone"
 #define NUM_REF_CELL_SUBZONES_ARRAY_DESCRIPTION "nszNumRefCszsForZone"
 #define CELL_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION "cszInclPtnBitArray"
 #define NODE_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION "nszInclPtnBitArray"
static uint16_t const ___4340 = 0x0600; static uint16_t const ___2723 = 0x0601; static uint16_t const USER_FACE_NBR_FILE_LOC_TAG = 0x0602;
 #define VAR_FILE_LOC_DESCRIPTION "varFileLocs"
 #define NODE_MAP_FILE_LOC_DESCRIPTION "nodeMapFileLoc"
 #define USER_FACE_NBR_FILE_LOC_DESCRIPTION "userFaceNbrFileLoc"
static uint16_t const NUM_PARTITIONS_TAG = 0x0700; static uint16_t const PARTITION_FILE_NUM_TAG = 0x0701; static uint16_t const PARTITION_FILE_LOC_TAG = 0x0702; static uint16_t const PARTITION_NUM_CELLS_TAG = 0x0703; static uint16_t const PARTITION_NUM_NODES_TAG = 0x0704; static uint16_t const VAR_PARTITION_MIN_MAX_TAG = 0x0705; static uint16_t const PARTITION_MIN_NODE_NUMBERS_TAG = 0x0706; static uint16_t const PARTITION_MAX_NODE_NUMBERS_TAG = 0x0707;
 #define NUM_PARTITIONS_DESCRIPTION "numPartitions"
 #define PARTITION_FILE_NUM_DESCRIPTION "partitionFileNums"
 #define PARTITION_FILE_LOC_DESCRIPTION "partitionFileLocs"
 #define PARTITION_NUM_CELLS_DESCRIPTION "partitionNumCells"
 #define PARTITION_NUM_NODES_DESCRIPTION "partitionNumNodes"
 #define VAR_PARTITION_MIN_MAX_DESCRIPTION "partitionMinMaxes"
 #define PARTITION_MIN_NODE_NUMBERS_DESCRIPTION "partitionMinNodeNumbers"
 #define PARTITION_MAX_NODE_NUMBERS_DESCRIPTION "partitionMaxNodeNumbers"
 #define ___4687 boost::assign::map_list_of<uint16_t, std::string> \
 ((uint16_t)___2782, NUM_CELL_SUBZONES_DESCRIPTION) \
 ((uint16_t)___2822, NUM_NODE_SUBZONES_DESCRIPTION) \
 ((uint16_t)CSZ_CONNECT_FILE_LOC_TAG, CSZ_CONNECT_FILE_LOC_DESCRIPTION) \
 ((uint16_t)NSZ_CONNECT_FILE_LOC_TAG, NSZ_CONNECT_FILE_LOC_DESCRIPTION) \
 ((uint16_t)UNUSED1_ZONE_TAG, ZONE_UNUSED1_TAG_DESCRIPTION) \
 ((uint16_t)NUM_REF_PARTITIONS_TAG, NUM_REF_PARTITIONS_DESCRIPTION) \
 ((uint16_t)REF_PARTITIONS_TAG, REF_PARTITIONS_DESCRIPTION) \
 ((uint16_t)NUM_REF_NODE_SUBZONES_TAG, NUM_REF_NODE_SUBZONES_DESCRIPTION) \
 ((uint16_t)NUM_REF_CELL_SUBZONES_TAG, NUM_REF_CELL_SUBZONES_DESCRIPTION) \
 ((uint16_t)CELL_SUBZONE_INCLUDES_PARTITIONS_TAG, CELL_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION) \
 ((uint16_t)NODE_SUBZONE_INCLUDES_PARTITIONS_TAG, NODE_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION) \
 ((uint16_t)___4340, VAR_FILE_LOC_DESCRIPTION) \
 ((uint16_t)___2723, NODE_MAP_FILE_LOC_DESCRIPTION) \
 ((uint16_t)USER_FACE_NBR_FILE_LOC_TAG, USER_FACE_NBR_FILE_LOC_DESCRIPTION) \
 ((uint16_t)NUM_PARTITIONS_TAG, NUM_PARTITIONS_DESCRIPTION) \
 ((uint16_t)PARTITION_FILE_NUM_TAG, PARTITION_FILE_NUM_DESCRIPTION) \
 ((uint16_t)PARTITION_FILE_LOC_TAG, PARTITION_FILE_LOC_DESCRIPTION) \
 ((uint16_t)PARTITION_NUM_CELLS_TAG, PARTITION_NUM_CELLS_DESCRIPTION) \
 ((uint16_t)PARTITION_NUM_NODES_TAG, PARTITION_NUM_NODES_DESCRIPTION) \
 ((uint16_t)VAR_PARTITION_MIN_MAX_TAG, VAR_PARTITION_MIN_MAX_DESCRIPTION) \
 ((uint16_t)PARTITION_MIN_NODE_NUMBERS_TAG, PARTITION_MIN_NODE_NUMBERS_DESCRIPTION) \
 ((uint16_t)PARTITION_MAX_NODE_NUMBERS_TAG, PARTITION_MAX_NODE_NUMBERS_DESCRIPTION)
 #define PARTITION_TAG_DESCRIPTIONS_MAP_INITIALIZER boost::assign::map_list_of<uint16_t, std::string> \
 ((uint16_t)CSZ_CONNECT_FILE_LOC_TAG, CSZ_CONNECT_FILE_LOC_DESCRIPTION) \
 ((uint16_t)NSZ_CONNECT_FILE_LOC_TAG, NSZ_CONNECT_FILE_LOC_DESCRIPTION) \
 ((uint16_t)NUM_REF_PARTITIONS_TAG, NUM_REF_PARTITIONS_DESCRIPTION) \
 ((uint16_t)REF_PARTITIONS_TAG, REF_PARTITIONS_DESCRIPTION) \
 ((uint16_t)NUM_REF_NODE_SUBZONES_TAG, NUM_REF_NODE_SUBZONES_DESCRIPTION) \
 ((uint16_t)NUM_REF_CELL_SUBZONES_TAG, NUM_REF_CELL_SUBZONES_DESCRIPTION) \
 ((uint16_t)CELL_SUBZONE_INCLUDES_PARTITIONS_TAG, CELL_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION) \
 ((uint16_t)NODE_SUBZONE_INCLUDES_PARTITIONS_TAG, NODE_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION) \
 ((uint16_t)___4340, VAR_FILE_LOC_DESCRIPTION) \
 ((uint16_t)___2782, NUM_CELL_SUBZONES_DESCRIPTION) \
 ((uint16_t)___2822, NUM_NODE_SUBZONES_DESCRIPTION)
 #define NSZ_CONNECT_FILE_LOC_ARRAY_DESCRIPTION "nszConnFileLocsForZone"
 #define CSZ_CONNECT_FILE_LOC_ARRAY_DESCRIPTION "cszConnFileLocsForZone"
 #define ___4339 "varFileLocsForZone"
 #define CELL_SUBZONE_DATA_MIN_MAX_ARRAY_DESCRIPTION_INPUT "cszDataMinMaxes"
 #define NODE_SUBZONE_DATA_MIN_MAX_ARRAY_DESCRIPTION_INPUT "nszDataMinMaxes"
 #define CELL_SUBZONE_DATA_MIN_MAX_ARRAY_DESCRIPTION "cszDataMinMaxesForVar"
 #define NODE_SUBZONE_DATA_MIN_MAX_ARRAY_DESCRIPTION "nszDataMinMaxesForVar"
 #define SZ_CONNECT_REF_PTN_OFFSET_2BIT_ARRAY_DESCRIPTION "refPtn2BitOffsets"
 #define SZ_CONNECT_REF_PTN_OFFSET_4BIT_ARRAY_DESCRIPTION "refPtn4BitOffsets"
 #define SZ_CONNECT_REF_PTN_OFFSET_8BIT_ARRAY_DESCRIPTION "refPtn8BitOffsets"
 #define SZ_CONNECT_REF_PTN_OFFSET_16BIT_ARRAY_DESCRIPTION "refPtn16BitOffsets"
 #define NSZ_CONNECT_NUM_REF_CSZ "numRefCellSubzones"
 #define NSZ_CONNECT_REF_CSZS "refCellSubzones"
 #define CSZ_CONNECT_NUM_REF_NSZ "numRefNodeSubzones"
 #define CSZ_CONNECT_REF_NSZ_ARRAY_DESCRIPTION "refNodeSubzones"
 #define CSZ_CONNECT_REF_NSZ_OFFSET_2BIT_ARRAY_DESCRIPTION "refNsz2BitOffsets"
 #define CSZ_CONNECT_REF_NSZ_OFFSET_4BIT_ARRAY_DESCRIPTION "refNsz4BitOffsets"
 #define CSZ_CONNECT_REF_NSZ_OFFSET_8BIT_ARRAY_DESCRIPTION "refNsz8BitOffsets"
 #define CSZ_CONNECT_REF_NSZ_OFFSET_16BIT_ARRAY_DESCRIPTION "refNsz16BitOffsets"
 #define CSZ_CONNECT_NSZ_OFFSET_ARRAY_DESCRIPTION "nszOffsets"
 #define ___4703 'o'
 #define ZONETYPE_FETRIANGLE_CHAR 'r'
 #define ZONETYPE_FEQUAD_CHAR 'q'
 #define ZONETYPE_FETETRA_CHAR 't'
 #define ZONETYPE_FEBRICK_CHAR 'b'
 #define ZONETYPE_FELINESEG_CHAR 'l'
 #define ZONETYPE_FEPOLYGON_CHAR 'g'
 #define ZONETYPE_FEPOLYHEDRON_CHAR 'p'
 #define ZONETYPE_FEMIXEDVOLUME_CHAR 'm'
 #define ALT_TAG_DESCRIPTIONS_MAP_INITIALIZER boost::assign::map_list_of<std::string, std::string> \
 (NSZ_CONNECT_FILE_LOC_DESCRIPTION, "nszConnectivityFileLocs")   \
 (CSZ_CONNECT_FILE_LOC_DESCRIPTION, "cszConnectivityFileLocs")   \
 (NSZ_CONNECT_FILE_LOC_ARRAY_DESCRIPTION, "nszConnectivityFileLocsForZone")   \
 (CSZ_CONNECT_FILE_LOC_ARRAY_DESCRIPTION, "cszConnectivityFileLocsForZone")   \
 (VAR_ZONE_IS_PASSIVE_FILE_LOC_DESCRIPTION, "passive")   \
 (VAR_ZONE_SHARE_FILE_LOC_DESCRIPTION, "share")

 #pragma once
#include "MASTER.h"
#include "GLOBAL.h"
#include "CodeContract.h"
#include "ItemAddress.h"
namespace tecplot { namespace ___3931 { typedef int32_t ___4634; ___4634 const MAX_NUM_ZONES = ___2389; ___4634 const BAD_ZONE_INDEX = ___4634(0x7FFFFFFF); typedef int32_t ___4350; ___4350 const MAX_NUM_VARS = ___2389; ___4350 const BAD_VAR_INDEX = ___4350(0x7FFFFFFF); typedef int64_t ___81; ___81 const MAX_ANY_INDEX = ___81(0x7FFFFFFFFFFFFFFEll); ___81 const BAD_ANY_INDEX = ___81(0x7FFFFFFFFFFFFFFFll); typedef int64_t ___463; ___463 const MAX_NUM_CELLS = ___463(0x7FFFFFFFFFFFFFFEll); ___463 const BAD_CELL_INDEX = ___463(0x7FFFFFFFFFFFFFFFll); typedef int64_t ___2716; ___2716 const MAX_NUM_NODES = ___2716(0x7FFFFFFFFFFFFFFEll); ___2716 const BAD_NODE_INDEX = ___2716(0x7FFFFFFFFFFFFFFFll); typedef int32_t ___680; ___680 const MAX_NUM_CELL_CORNERS = 8; ___680 const BAD_CORNER_INDEX = ___680(-1); ___680 const NUM_IJK_CELL_CORNERS = 8; ___680 const NUM_TETRA_CELL_CORNERS = 4; ___680 const NUM_BRICK_CELL_CORNERS = 8; typedef int32_t FaceIndex_t; FaceIndex_t const MAX_NUM_CELL_FACES = 6; FaceIndex_t const BAD_FACE_INDEX = FaceIndex_t(-1); FaceIndex_t const NUM_IJK_CELL_FACES = 6; FaceIndex_t const NUM_TETRA_CELL_FACES = 4; FaceIndex_t const NUM_BRICK_CELL_FACES = 6; typedef uint64_t ___1391; ___1391 const ___330 = ___1391(-1);
 #define VALID_FILE_LOC(fileLoc) ( (fileLoc) != ___330 )
typedef uint16_t RefSubzoneOffset_t; RefSubzoneOffset_t const BAD_REFSZ_INDEX = RefSubzoneOffset_t(-1); }}

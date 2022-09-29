#include "ClassicFEZoneFaceNeighborGenerator.h"
#include "AltTecUtil.h"
#include "CodeContract.h"
namespace tecplot { namespace ___3931 { namespace { int32_t cellCountForZone(___37& ___36, ___4634 zone) { REQUIRE(zone >= 0); ___1842 ___1841; ___36.___4613(zone + 1, ___1841); return static_cast<int32_t>(___1841.___2103()); } int32_t faceCountForZone(___37& ___36, ___4634 zone) { REQUIRE(zone >= 0); ___1842 ___1841; ___36.___4613(zone + 1, ___1841); return static_cast<int32_t>(___1841.___2132()); } } ClassicFEZoneFaceNeighborGenerator::ClassicFEZoneFaceNeighborGenerator(class ___37& ___36) : FaceNeighborGeneratorAbstract(___36) {} ___372 ClassicFEZoneFaceNeighborGenerator::gatherUserFaceNeighbors(std::vector<int32_t>& userFaceNeighbors, ___4634 zone) const { REQUIRE(userFaceNeighbors.empty()); ___3499 ___1151 = 0; ___372 ___2037 = ___2335.___4612(&___1151); if (___2037) { ___1290 ___1272 = ___2335.___837(zone + 1); int32_t cellCount = cellCountForZone(___2335, zone); int32_t faceCount = faceCountForZone(___2335, zone); for(int32_t ___447 = 1; ___447 <= cellCount; ++___447) { for(int32_t face = 1; face <= faceCount; ++face) appendUserFaceNeighborsForCellFace(userFaceNeighbors, ___1272, ___1151, zone, ___447, face); } } ___2335.___3482(&___1151); return ___2037; } }}
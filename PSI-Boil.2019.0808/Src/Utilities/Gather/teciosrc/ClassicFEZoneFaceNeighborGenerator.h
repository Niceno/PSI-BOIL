 #pragma once
#include "FaceNeighborGeneratorAbstract.h"
namespace tecplot { namespace ___3931 { class ClassicFEZoneFaceNeighborGenerator : public FaceNeighborGeneratorAbstract { public: explicit ClassicFEZoneFaceNeighborGenerator(class ___37& ___36); virtual ___372 gatherUserFaceNeighbors(std::vector<int32_t>& userFaceNeighbors, ___4634 zone) const; }; }}

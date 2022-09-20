 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3931 { class FaceNeighborGeneratorAbstract; class ClassicZoneFaceNeighborWriter { public: ClassicZoneFaceNeighborWriter( FaceNeighborGeneratorAbstract& faceNeighborGenerator, ___4634 zone, ___4634 ___341); ___372 write(class FileWriterInterface& file); uint64_t sizeInFile(bool ___2000); private: FaceNeighborGeneratorAbstract const& m_faceNeighborGenerator; ___4634 const ___2675; ___4634 const m_baseZone; std::string const m_zoneNumberLabel; std::vector<int32_t> m_userFaceNeighbors; }; }}

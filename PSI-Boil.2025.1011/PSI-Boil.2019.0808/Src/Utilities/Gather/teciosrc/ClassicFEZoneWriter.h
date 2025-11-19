 #pragma once
#include "ClassicFEZoneConnectivityWriter.h"
#include "ClassicFEZoneFaceNeighborGenerator.h"
#include "ClassicZoneFaceNeighborWriter.h"
#include "ClassicZoneWriterAbstract.h"
namespace tecplot { namespace ___3931 { class ClassicFEZoneWriter : public ClassicZoneWriterAbstract { UNCOPYABLE_CLASS(ClassicFEZoneWriter); public: ClassicFEZoneWriter( ItemSetIterator&              varIter, ___4634                   zone, ___4634                   ___341, std::vector<___372> const& ___4562, ___372                     ___4497, ___37&                   ___36); virtual ~ClassicFEZoneWriter(); private: virtual uint64_t zoneConnectivityFileSize(bool ___2000); virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile); ClassicFEZoneConnectivityWriter m_connectivityWriter; ClassicFEZoneFaceNeighborGenerator m_faceNeighborGenerator; ClassicZoneFaceNeighborWriter m_faceNeighborWriter; std::string m_zoneNumberLabel; }; }}

 #pragma once
#include "ClassicZoneFaceNeighborWriter.h"
#include "ClassicOrderedZoneFaceNeighborGenerator.h"
#include "ClassMacros.h"
#include "ClassicZoneWriterAbstract.h"
namespace tecplot { namespace ___3931 { class ItemSetIterator; class ClassicOrderedZoneWriter : public ClassicZoneWriterAbstract { UNCOPYABLE_CLASS(ClassicOrderedZoneWriter); public: ClassicOrderedZoneWriter( ItemSetIterator&              varIter, ___4634                   zone, ___4634                   ___341, std::vector<___372> const& ___4562, ___372                     ___4497, ___37&                   ___36); virtual ~ClassicOrderedZoneWriter(); private: virtual uint64_t zoneConnectivityFileSize(bool ___2000); virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile); ClassicOrderedZoneFaceNeighborGenerator m_faceNeighborGenerator; ClassicZoneFaceNeighborWriter m_faceNeighborWriter; ___372 ___4512( FileWriterInterface& file, ValueLocation_e      ___4324, ___4350           ___4334); }; }}

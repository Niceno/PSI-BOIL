 #pragma once
#include "ClassicZoneHeaderWriter.h"
#include "ClassicZoneVariableWriter.h"
#include "ZoneWriterAbstract.h"
#include "ClassicZoneFileLocations.h"
namespace tecplot { namespace ___3931 { class ClassicZoneWriterAbstract : public ___4707 { public: ClassicZoneWriterAbstract( ItemSetIterator&              varIter, ___4634                   zone, ___4634                   ___341, std::vector<___372> const& ___4562, ___372                     ___4497, ___37&                   ___36); virtual ~ClassicZoneWriterAbstract(); protected: ClassicZoneVariableWriter m_variableWriter; ClassicZoneHeaderWriter m_headerWriter; ClassicZoneFileLocations m_zoneFileLocations; private: virtual uint64_t zoneConnectivityFileSize(bool ___2000) = 0; virtual uint64_t zoneDataFileSize(bool ___2000); virtual uint64_t zoneHeaderFileSize(bool ___2000); virtual ___372 writeZoneData(FileWriterInterface& szpltFile); virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile) = 0; virtual ___372 writeZoneHeader(FileWriterInterface& szpltFile); }; }}

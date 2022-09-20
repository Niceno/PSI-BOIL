#include "ClassicOrderedZoneWriter.h"
#include "ThirdPartyHeadersBegin.h"
#include <boost/assign.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "AltTecUtil.h"
#include "checkPercentDone.h"
#include "FieldData.h"
#include "fileStuff.h"
#include "ItemSetIterator.h"
#include "writeValueArray.h"
namespace tecplot { namespace ___3931 { ClassicOrderedZoneWriter::ClassicOrderedZoneWriter( ItemSetIterator&              varIter, ___4634                   zone, ___4634                   ___341, std::vector<___372> const& ___4562, ___372                     ___4497, ___37&                   ___36) : ClassicZoneWriterAbstract(varIter, zone, ___341, ___4562, ___4497, ___36) , m_faceNeighborGenerator(___36) , m_faceNeighborWriter(m_faceNeighborGenerator, zone, ___341) {} ClassicOrderedZoneWriter::~ClassicOrderedZoneWriter() {} ___372 ClassicOrderedZoneWriter::writeZoneConnectivity(FileWriterInterface& szpltFile) { ___372 ___2037 = ___4224; m_zoneFileLocations.___2496 = ___330; if (m_writeConnectivity) { m_zoneFileLocations.___2661 = szpltFile.fileLoc(); ___2037 = m_faceNeighborWriter.write(szpltFile); } else { m_zoneFileLocations.___2661 = ___330; } return ___2037; } uint64_t ClassicOrderedZoneWriter::zoneConnectivityFileSize(bool ___2000) { if (m_writeConnectivity) return m_faceNeighborWriter.sizeInFile(___2000); else return 0; } }}

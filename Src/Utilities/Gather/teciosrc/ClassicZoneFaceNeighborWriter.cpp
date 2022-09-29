#include "ClassicZoneFaceNeighborWriter.h"
#include "FaceNeighborGeneratorAbstract.h"
#include "FileWriterInterface.h"
#include "fileStuff.h"
#include "SzlFileLoader.h"
#include "writeValueArray.h"
namespace tecplot { namespace ___3931 { namespace { char const* USER_FACE_NEIGHBOR_COUNT_LABEL = "userFaceNeighborCount"; char const* USER_FACE_NEIGHBOR_MARKER_LABEL = "userFaceNeighborMarker*"; char const* USER_FACE_NEIGHBORS_LABEL = "UserFaceNeighbors"; } ClassicZoneFaceNeighborWriter::ClassicZoneFaceNeighborWriter( FaceNeighborGeneratorAbstract& faceNeighborGenerator, ___4634 zone, ___4634 ___341) : m_faceNeighborGenerator(faceNeighborGenerator) , ___2675(zone) , m_baseZone(___341) , m_zoneNumberLabel(SZPLT_ZONE_NUM_DESCRIPTION) {} ___372 ClassicZoneFaceNeighborWriter::write( FileWriterInterface& szpltFile) { REQUIRE(szpltFile.___2039()); ___372 ___2037 = ___4224; if (m_userFaceNeighbors.empty()) m_faceNeighborGenerator.gatherUserFaceNeighbors(m_userFaceNeighbors, ___2675); if (szpltFile.___2000()) { ___2037 = writeValue<uint32_t, false, 0>(szpltFile, USER_FACE_NEIGHBOR_MARKER_LABEL, SZPLT_USER_FACE_NEIGHBOR_MARKER) && writeValue<uint32_t, false, 0>(szpltFile, m_zoneNumberLabel.c_str(), (___2675 - m_baseZone + 1)); } ___1963 faceNeighborVector; ___2037 = ___2037 && writeValue<uint64_t, false, 0>(szpltFile, USER_FACE_NEIGHBOR_COUNT_LABEL, m_userFaceNeighbors.size()); if (___2037 && m_userFaceNeighbors.size() > 0) { ___2037 = ___2037 && faceNeighborVector.alloc(m_userFaceNeighbors.size()); if (___2037) { for(size_t i = 0; i < m_userFaceNeighbors.size(); ++i) faceNeighborVector[i] = m_userFaceNeighbors[i]; ___2037 = ___2037 && ___4561<int32_t, false, 0>(szpltFile, USER_FACE_NEIGHBORS_LABEL, ___2743, m_userFaceNeighbors.size(), &faceNeighborVector[0]); } } return ___2037; } uint64_t ClassicZoneFaceNeighborWriter::sizeInFile(bool ___2000) { if (m_userFaceNeighbors.empty()) m_faceNeighborGenerator.gatherUserFaceNeighbors(m_userFaceNeighbors, ___2675); uint64_t ___3356 = 0; if (___2000) ___3356 += 2 * valueSizeInFile<uint32_t, false>(___2000); ___3356 += valueSizeInFile<uint64_t, false>(___2000); if (m_userFaceNeighbors.size() > 0) ___3356 += arraySizeInFile<int32_t, false>(m_userFaceNeighbors.size(), ___2000); return ___3356; } }}
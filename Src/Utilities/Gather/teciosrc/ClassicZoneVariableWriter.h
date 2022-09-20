 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <string>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3931 { class ___37; class FileWriterInterface; class ___1350; class ItemSetIterator; class ClassicZoneVariableWriter { public: ClassicZoneVariableWriter( ItemSetIterator& varIter, ___4634      zone, ___4634      ___341, ___37&      ___36); static uint64_t varHeaderSizeInFile(bool ___2000); uint64_t varSizeInFile(___4350 ___4334, bool ___2000) const; ___372 writeVarHeader( FileWriterInterface& file, ValueLocation_e      ___4324, ___4350           ___4334); ___372 writeVar( FileWriterInterface& szpltFile, ___1350 const&     ___1349); private: ItemSetIterator&  m_varIter; ___4634 const ___2675; ___4634 const m_baseZone; ___37&       ___2335; std::string const m_zoneNumberLabel; }; }}

#include "ZoneHeaderWriterAbstract.h"
#include "ThirdPartyHeadersBegin.h"
#include <boost/assign.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "AltTecUtil.h"
#include "FileWriterInterface.h"
#include "ItemSetIterator.h"
#include "SzlFileLoader.h"
#include "writeValueArray.h"
namespace tecplot { namespace ___3931 { ZoneHeaderWriterAbstract::ZoneHeaderWriterAbstract( ItemSetIterator&   varIter, ___4634        zone, ___4634        ___341, ___37&        ___36, std::string const& zoneMarkerLabel  , uint32_t           zoneMarker  , std::string const& zoneNumberLabel  , std::string const& zoneNumberSuffix  ) : m_varIter(varIter) , ___2675(zone) , m_baseZone(___341) , ___2335(___36) , m_zoneMarkerLabel(zoneMarkerLabel) , m_zoneMarker(zoneMarker) , m_zoneNumberLabel(zoneNumberLabel) , m_zoneNumberSuffix(zoneNumberSuffix) {} ZoneHeaderWriterAbstract::~ZoneHeaderWriterAbstract() {} void ZoneHeaderWriterAbstract::setMarkerAndLabels( std::string const& zoneMarkerLabel, uint32_t zoneMarker, std::string const& zoneNumberLabel, std::string const& zoneNumberSuffix) { m_zoneMarkerLabel = zoneMarkerLabel; m_zoneMarker = zoneMarker; m_zoneNumberLabel = zoneNumberLabel; m_zoneNumberSuffix = zoneNumberSuffix; } ___372 ZoneHeaderWriterAbstract::___4563(FileWriterInterface& file, ___3943 const& ___3942) const { REQUIRE(file.___2039()); ___372 ___2037 = ___4224; if (file.___2000()) { ___2037 = writeValue<uint32_t, false, 0>(file, m_zoneMarkerLabel.c_str(), m_zoneMarker) && writeValue<uint32_t, false, 0>(file, m_zoneNumberLabel.c_str(), (___2675 - m_baseZone) + 1); } UInt16Array tagArray; ___2037 = ___2037 && tagArray.alloc(___3942.size()); if (___2037) { size_t numTags = 0; for (___3943::const_iterator ___4312 = ___3942.begin(); ___4312 != ___3942.end(); ++___4312) tagArray[numTags++] = ___4312->first; ___2037 = ___2037 && writeValue<uint16_t, false, 0>(file, ___2837, static_cast<uint16_t>(___3942.size())) && ___4561<uint16_t, false, 0>(file, ___3941, ___2743, numTags, &tagArray[0]); } if (___2037) { ___3940 tagDescriptionMap = ___4687; for (___3943::const_iterator ___4312 = ___3942.begin(); ___4312 != ___3942.end(); ++___4312) ___2037 = ___2037 && writeValue<uint64_t, true, 0>(file, tagDescriptionMap[___4312->first].c_str(), ___4312->second); } ENSURE(VALID_BOOLEAN(___2037)); return ___2037; } uint64_t ZoneHeaderWriterAbstract::zoneHeaderTagsSizeInFile(uint16_t numTags, bool ___2000) const { uint64_t ___3356 = 0; if (___2000) ___3356 += 2 * valueSizeInFile<uint32_t, false>(___2000); ___3356 += valueSizeInFile<uint16_t, false>(___2000); ___3356 += arraySizeInFile<uint16_t, false>(static_cast<size_t>(numTags), ___2000); ___3356 += numTags * valueSizeInFile<uint64_t, true  >(___2000); return ___3356; } std::string ZoneHeaderWriterAbstract::appendZoneSuffix(std::string const& str) const { return str + m_zoneNumberSuffix; } }}
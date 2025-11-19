 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_ptr.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
#include "ZoneInfoCache.h"
#include "ZoneVarMetadata.h"
namespace tecplot { namespace ___3931 { class ___37; class FileWriterInterface; class ItemSetIterator; class DataSetWriter { public: DataSetWriter( ___37*               ___36, ___3499                    varsToWrite, ___3499                    zonesToWrite, ___1842 const&                maxIJKSubzoneSize, ___2088::ItemOffset_t maxFESubzoneSize, bool                      flushToDisk = false); virtual ~DataSetWriter(); virtual ___372 writeDataSet( FileWriterInterface& szpltFile, ___1390&        szpltZoneHeaderFileLocs); void replaceDataSource( ___37* ___36, ___3499      varsToWrite, ___3499      zonesToWrite); void clear( int32_t        numZonesToRetain, int32_t const* zonesToRetain); ___4705 const& ___4704() { return *m_zoneVarMetadata; } protected: void getZoneSharing( std::vector<___372>& ___4562, ___372&              ___4497, ___4634             zone, ___4634             ___341, DataFileType_e          ___842) const; ___37*                        ___2335; boost::scoped_ptr<ItemSetIterator> m_varIter; boost::scoped_ptr<ItemSetIterator> m_zoneIter; ZoneInfoCache                      ___2678; boost::scoped_ptr<___4705> m_zoneVarMetadata; bool const                         m_flushingToDisk; }; }}

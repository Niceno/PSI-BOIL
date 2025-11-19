 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <vector>
#include <boost/shared_ptr.hpp>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3931 { class ___37; class ZoneInfoCache; class ___4707; class ItemSetIterator; class ___4709 { public: ___4709( ZoneInfoCache& zoneInfoCache, ___37& ___36); boost::shared_ptr<___4707> ___4706( ItemSetIterator&              varIter, ___4634                   zone, ___4634                   ___341, std::vector<___372> const& ___4562, ___372                     ___4497); protected: ZoneInfoCache& ___2678; ___37& ___2335; }; }}

#include "checkPercentDone.h"
#include "stringformat.h"
#include "SzlFileLoader.h"
#include "zoneUtil.h"
namespace tecplot { namespace ___3931 { ___372 ___483( SzPltWriteOperation_e szPltWriteOperation, ___37& ___36, ___4350 ___4334, ___4634 zone) { REQUIRE(VALID_ENUM(szPltWriteOperation, SzPltWriteOperation_e)); ___4350  const ___2841 = ___36.___888(); ___4634 const ___2844 = ___36.___889(); REQUIRE(___4334<___2841); REQUIRE(zone<___2844); ___372 const isOrderedZone = ___4639(___36, zone); ___372 const isSzlZone = isOrderedZone ? ___4640(___36, zone) : ___4642(___36, zone); int const determineSubzonesWork = (isSzlZone && !isOrderedZone ? 1 : 0); int const writeConnectivityWork = (isOrderedZone ? 0 : 1); int const calcSzMinMaxWork = (isSzlZone ? 1 : 0); int const writeFieldDataWork = 1; int completedWork = 0; int currentWork = 0; std::string status; if ( szPltWriteOperation ==  SzPltWriteOperation_DetermineSubzones) { ___476(!isOrderedZone); ___476(isSzlZone); ___476(___4334==0); status = "Determining subzones for zone " + ___4185(uint64_t(zone)+1) + "..."; completedWork = 0; currentWork = determineSubzonesWork; } else if ( szPltWriteOperation == SzPltWriteOperation_WriteConnectivity ) { ___476(!isOrderedZone); ___476(___4334==0); status = "Writing connectivity for zone " + ___4185(uint64_t(zone)+1) + "..."; if ( isSzlZone ) completedWork = determineSubzonesWork; else completedWork = determineSubzonesWork + writeConnectivityWork + calcSzMinMaxWork; currentWork = writeConnectivityWork; } else if ( szPltWriteOperation == SzPltWriteOperation_CalcSubzoneMinMax ) { ___476(isSzlZone); status = "Calculating subzone min-maxes for zone " + ___4185(uint64_t(zone)+1) + " variable " + ___4185(uint64_t(___4334)+1) + "..."; completedWork = determineSubzonesWork + writeConnectivityWork; currentWork = calcSzMinMaxWork; } else { ___476(szPltWriteOperation == SzPltWriteOperation_WriteFieldData); status = "Writing field data for zone " + ___4185(uint64_t(zone)+1) + " variable " + ___4185(uint64_t(___4334)+1) + "..."; if ( isSzlZone ) completedWork = determineSubzonesWork + writeConnectivityWork + calcSzMinMaxWork; else completedWork = 0; currentWork = writeFieldDataWork; } int const totalPerVarWork = determineSubzonesWork + writeConnectivityWork + calcSzMinMaxWork + writeFieldDataWork; ___476(completedWork+currentWork <= totalPerVarWork); int const percentDone = ( (100*completedWork*___2841 + 100*currentWork*___4334)/totalPerVarWork + 100*zone*___2841)/(___2841*___2844); ___476(0<=percentDone && percentDone<=100); ___36.___3776(status.c_str()); ___372 isInterrupted = ___36.___3767(percentDone); REQUIRE(VALID_BOOLEAN(isInterrupted)); return isInterrupted; } }}

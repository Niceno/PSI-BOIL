 #pragma once
#include "basicTypes.h"
#include "ItemAddress.h"
namespace tecplot { namespace ___3931 { class FESubzonePartitionerInterface { public: virtual ~FESubzonePartitionerInterface() {} virtual ___463                  numCellsInZone() const = 0; virtual ___2088::SubzoneOffset_t ___2781() const = 0; virtual ___2088::ItemOffset_t    ___2780(___2088::SubzoneOffset_t ___467) const = 0; virtual ___463                  ___4606(___2088 ___449) const = 0; virtual ___2088                  szCoordinateAtZoneCell(___463 zoneCell) const = 0; virtual ___2716                  numNodesInZone() const = 0; virtual ___2088::SubzoneOffset_t ___2821() const = 0; virtual ___2088::ItemOffset_t    ___2820(___2088::SubzoneOffset_t ___2732) const = 0; virtual ___2716                  ___4655(___2088 nodeAddress) const = 0; virtual ___2088                  ___3922(___2716 ___4654) const = 0; virtual void                         setNodeSubzoneCoordinate(___2716 ___4654, ___2088 ___2757) = 0; }; }}

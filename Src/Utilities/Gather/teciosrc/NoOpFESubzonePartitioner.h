 #pragma once
#include "FESubzonePartitionerInterface.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3931 { class ___37; class NoOpFESubzonePartitioner : public FESubzonePartitionerInterface { public: NoOpFESubzonePartitioner( ___37& ___36, ___4634 zone) { ___1842 ___1841; ___36.___4613(zone + 1, ___1841); m_numCells = static_cast<___463>(___1841.___2103()); m_numNodes = static_cast<___2716>(___1841.i()); } virtual ~NoOpFESubzonePartitioner() {} virtual ___463                  numCellsInZone() const { return m_numCells; } virtual ___2088::SubzoneOffset_t ___2781() const { return 0; } virtual ___2088::ItemOffset_t    ___2780(___2088::SubzoneOffset_t  ) const { ___476(___1303); return 0; } virtual ___463                  ___4606(___2088  ) const { ___476(___1303); return 0; } virtual ___2088                  szCoordinateAtZoneCell(___463  ) const { ___476(___1303); return ___2088(0,0,0); } virtual ___2716                  numNodesInZone() const { return m_numNodes; } virtual ___2088::SubzoneOffset_t ___2821() const { return 0; } virtual ___2088::ItemOffset_t    ___2820(___2088::SubzoneOffset_t  ) const { ___476(___1303); return 0; } virtual ___2716                  ___4655(___2088  ) const { ___476(___1303); return 0; } virtual ___2088                  ___3922(___2716  ) const { ___476(___1303); return ___2088(0,0,0); } virtual void                         setNodeSubzoneCoordinate(___2716  , ___2088 /*___2757*/) { ___476(___1303); } private: ___463 m_numCells; ___2716 m_numNodes; }; }}

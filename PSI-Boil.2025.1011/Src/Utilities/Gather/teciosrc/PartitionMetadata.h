 #pragma once
#include "ClassMacros.h"
#include "SzlFileLoader.h"
namespace tecplot { namespace ___3931 { class PartitionMetadata { UNCOPYABLE_CLASS(PartitionMetadata); public: ___1390 m_cszConnectivityFileLocs; ___1390 m_nszConnectivityFileLocs; ___2088::___2978 m_numRefPartitions; PartitionArray           m_refPartitions; RefSubzoneOffsetArray m_cszNumRefNszs; RefSubzoneOffsetArray m_nszNumRefCszs; UInt8Array m_cszIncludesPtn; UInt8Array m_nszIncludesPtn; ___1390 m_cszMinMaxFileLocs; ___1390 m_nszMinMaxFileLocs; ___1390 m_szDataStartFileLocs; FileLoc2DArray m_varSzFileLoc; public: PartitionMetadata() {} }; }}

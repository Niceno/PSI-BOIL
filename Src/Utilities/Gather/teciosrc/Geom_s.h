 #pragma once
#include "SzlFileLoader.h"
#include "fileio.h"
namespace tecplot { namespace tecioszl { struct ___1554 : public tecplot::___3931::___1554 { public: static ___1554 invalidGeom() { return ___1554( 0.0, 0.0, 0.0, ___661, ___1303, 1, ___364, ___4453, ___1303, GeomType_Invalid, ___2269, 2.0, 0.1, 72, ___192, ___181, 5.0, 12.0, ___3444, ___496, std::vector<std::vector<tecplot::___3931::___4578> >(), "" ); } ___1554( double ___4572, double ___4589, double ___4713, CoordSys_e ___3157, ___372 ___2002, ___1170 zone, ___514 color, ___514 ___1409, ___372 ___2020, GeomType_e ___1649, LinePattern_e ___2261, double ___2984, double ___2287, uint16_t ___2791, ArrowheadStyle_e arrowheadStyle, ArrowheadAttachment_e arrowheadAttachment, double arrowheadSize, double arrowheadAngle, Scope_e ___3440, Clipping_e ___493, std::vector<std::vector<tecplot::___3931::___4578> > const& ___1570, std::string const& ___2328) : tecplot::___3931::___1554(___4572, ___4589, ___4713, ___3157, ___2002, zone, color, ___1409, ___2020, ___1649, ___2261, ___2984, ___2287, ___2791, arrowheadStyle, arrowheadAttachment, arrowheadSize, arrowheadAngle, ___3440, ___493, ___1570, ___2328) {} bool ___2065() { return VALID_ENUM(___2614, CoordSys_e) && VALID_ENUM(___2464, GeomType_e); } void writeToFile(___3931::FileWriterInterface& outputFile, bool ___4478) const { writeScalar(outputFile, ___2615[0], ___4478); writeScalar(outputFile, ___2615[1], ___4478); writeScalar(outputFile, ___2615[2], ___4478); writeScalar(outputFile, (uint32_t)___2614, ___4478); writeScalar(outputFile, ___2482, ___4478); writeScalar(outputFile, ___2675, ___4478); writeScalar(outputFile, ___2393, ___4478); writeScalar(outputFile, ___2460, ___4478); writeScalar(outputFile, ___2484, ___4478); writeScalar(outputFile, (uint32_t)___2464, ___4478); writeScalar(outputFile, (uint32_t)___2487, ___4478); writeScalar(outputFile, ___2613, ___4478); writeScalar(outputFile, ___2488, ___4478); writeScalar(outputFile, ___2500, ___4478); writeScalar(outputFile, (uint32_t)___2341, ___4478); writeScalar(outputFile, (uint32_t)___2339, ___4478); writeScalar(outputFile, ___2340, ___4478); writeScalar(outputFile, ___2338, ___4478); writeScalar(outputFile, (uint32_t)___2617, ___4478); writeScalar(outputFile, (uint32_t)___2392, ___4478); writeScalar(outputFile, (uint64_t)___2462.size(), ___4478); for (std::vector<std::vector<tecplot::___3931::___4578> >::const_iterator vec = ___2462.begin(); vec != ___2462.end(); ++vec) { writeScalar(outputFile, (uint64_t)vec->size(), ___4478); for (std::vector<tecplot::___3931::___4578>::const_iterator xyz = vec->begin(); xyz != vec->end(); ++xyz) { writeScalar(outputFile, xyz->x(), ___4478); writeScalar(outputFile, xyz->___4581(), ___4478);
writeScalar(outputFile, xyz->z(), ___4478); } } ___4542(outputFile, ___2489, ___4478); } uint64_t sizeInFile(bool ___4478) const { uint64_t sizeInFile = 0; sizeInFile += scalarSizeInFile(___2615[0], ___4478); sizeInFile += scalarSizeInFile(___2615[1], ___4478); sizeInFile += scalarSizeInFile(___2615[2], ___4478); sizeInFile += scalarSizeInFile((uint32_t)___2614, ___4478); sizeInFile += scalarSizeInFile(___2482, ___4478); sizeInFile += scalarSizeInFile(___2675, ___4478); sizeInFile += scalarSizeInFile(___2393, ___4478); sizeInFile += scalarSizeInFile(___2460, ___4478); sizeInFile += scalarSizeInFile(___2484, ___4478); sizeInFile += scalarSizeInFile((uint32_t)___2464, ___4478); sizeInFile += scalarSizeInFile((uint32_t)___2487, ___4478); sizeInFile += scalarSizeInFile(___2613, ___4478); sizeInFile += scalarSizeInFile(___2488, ___4478); sizeInFile += scalarSizeInFile(___2500, ___4478); sizeInFile += scalarSizeInFile((uint32_t)___2341, ___4478); sizeInFile += scalarSizeInFile((uint32_t)___2339, ___4478); sizeInFile += scalarSizeInFile(___2340, ___4478); sizeInFile += scalarSizeInFile(___2338, ___4478); sizeInFile += scalarSizeInFile((uint32_t)___2617, ___4478); sizeInFile += scalarSizeInFile((uint32_t)___2392, ___4478); sizeInFile += scalarSizeInFile((uint64_t)___2462.size(), ___4478); for (std::vector<std::vector<tecplot::___3931::___4578> >::const_iterator vec = ___2462.begin(); vec != ___2462.end(); ++vec) { sizeInFile += scalarSizeInFile((uint64_t)vec->size(), ___4478); for (std::vector<tecplot::___3931::___4578>::const_iterator xyz = vec->begin(); xyz != vec->end(); ++xyz) { sizeInFile += scalarSizeInFile(xyz->x(), ___4478); sizeInFile += scalarSizeInFile(xyz->___4581(), ___4478); sizeInFile += scalarSizeInFile(xyz->z(), ___4478); } } sizeInFile += stringSizeInFile(___2489, ___4478); return sizeInFile; } ___1554(___3931::___1397& inputFile, bool readASCII) { readScalar(inputFile, ___2615[0], readASCII); readScalar(inputFile, ___2615[1], readASCII); readScalar(inputFile, ___2615[2], readASCII); READ_ENUM(___2614, CoordSys_e, inputFile, readASCII); readScalar(inputFile, ___2482, readASCII); readScalar(inputFile, ___2675, readASCII); readScalar(inputFile, ___2393, readASCII); readScalar(inputFile, ___2460, readASCII); readScalar(inputFile, ___2484, readASCII); READ_ENUM(___2464, GeomType_e, inputFile, readASCII); READ_ENUM(___2487, LinePattern_e, inputFile, readASCII); readScalar(inputFile, ___2613, readASCII); readScalar(inputFile, ___2488, readASCII); readScalar(inputFile, ___2500, readASCII); READ_ENUM(___2341, ArrowheadStyle_e, inputFile, readASCII); READ_ENUM(___2339, ArrowheadAttachment_e, inputFile, readASCII); readScalar(inputFile, ___2340, readASCII);
readScalar(inputFile, ___2338, readASCII); READ_ENUM(___2617, Scope_e, inputFile, readASCII); READ_ENUM(___2392, Clipping_e, inputFile, readASCII); uint64_t length; readScalar(inputFile, length, readASCII); ___2462.resize((size_t)length); for(size_t i = 0; i < ___2462.size(); ++i) { readScalar(inputFile, length, readASCII); ___2462[i].reserve((size_t)length); double x, ___4581, z; for(uint64_t ___2103 = 0; ___2103 < length; ++___2103) { readScalar(inputFile, x, readASCII); readScalar(inputFile, ___4581, readASCII); readScalar(inputFile, z, readASCII); ___2462[i].push_back(___3931::___4578(x, ___4581, z)); } } readString(inputFile, ___2489, readASCII); } }; }}

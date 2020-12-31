#include "recorder.h"
//#include "sompackexporter.h"
#include <sstream>

namespace dredviz {

void Recorder::record(const DataMatrix& data)
{
   std::ostringstream filename;
   filename << filename_stem << counter << ".dat";

   //SOMPackExporter exporter(filename.str());
   //exporter.exportData(data);
   counter++;
}
}

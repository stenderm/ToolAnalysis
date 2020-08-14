#ifndef TestLAPPDSimWaveforms_H
#define TestLAPPDSimWaveforms_H

#include "Tool.h"

class TestLAPPDSimWaveforms: public Tool{

public:

  TestLAPPDSimWaveforms();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

private:


};


#endif

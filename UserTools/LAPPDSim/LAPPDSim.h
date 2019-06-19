#ifndef LAPPDSim_H
#define LAPPDSim_H

#include <string>
#include <iostream>

#include "Geometry.h"
#include "Detector.h"
#include "Tool.h"
#include "TFile.h"
#include "TTree.h"
#include "wcsimT.h"
#include "LAPPDresponse.hh"
#include "TBox.h"
#include "TApplication.h"
#include "LAPPDDisplay.h"
// #include "Hit.h"
// #include "LAPPDHit.h"

class LAPPDSim: public Tool {


 public:

  LAPPDSim();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();
  Waveform<double> SimpleGenPulse(vector<double> pulsetimes);

 private:

   //ROOT random number generator
   TRandom3* myTR;
   TString SimInput;
   TFile* _tf;
   std::vector <TH2D*> _LAPPD_histograms;
   int iter = 0;
   int _event_counter = 0;
   int _display_config;
   LAPPDDisplay* _display;
   Geometry* _geom = nullptr;
   std::map<unsigned long, Waveform<double> >* LAPPDWaveforms = nullptr;

};


#endif

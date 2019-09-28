#ifndef RingCounting_H
#define RingCounting_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "BeamStatus.h"
#include "TriggerClass.h"
#include "Detector.h"
#include "Geometry.h"
#include "Hit.h"
#include "LAPPDHit.h"
#include "CherenkovParameters.h"

class RingCounting: public Tool {


 public:

  RingCounting();
  //virtual ~RingCounting(); //ToDo: Implementation
  double returnMassInGeV(int PDGCode);
  double particleVelocity(int particleNumber);
  bool radiatesCherenkov(int particleNumber);
  double cherenkovAngle(int particleNumber);
  void readAndPrintConfigParameters();
  void readandPrintGeometryParameters();
  void readMCdata();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:
   TFile* _file = nullptr;
   TH1D* _ring_events = nullptr;
   TApplication* ringCountingApp = nullptr;
   TCanvas* ringCountingCanv=nullptr;
//Declare ConfigFile parameters
   int verbose;
   int _csv_config;
//Declare Geometry variables
  Position detector_center;
  double tank_height;
  double tank_radius;
  int n_tank_pmts;
  int n_mrd_pmts;
  int n_veto_pmts;
  int n_lappds;
  int channelkey_pmt[200];
  double charge_pmt[200];
  Position position_pmt[200];
  double x_pmt[200];
  double y_pmt[200];
  double z_pmt[200];
  double charge[200];
  double time[200];
  double charge_lappd[200];
  std::vector<std::vector<Position>> hits_LAPPDs;
  std::vector<std::vector<double>> time_lappd;
  Geometry* geom = nullptr;
//Declare MC parameters
  int evnum;
  int runnumber;
  int subrunnumber;
  int terminate_execution;  //for continuous execution of multiple events, prompt user input
  std::map<ChannelKey,vector<Hit>>* TDCData=nullptr;
  TimeClass* EventTime=nullptr;
  BeamStatusClass* BeamStatus=nullptr;
  std::vector<TriggerClass>* TriggerData=nullptr;
  std::map<unsigned long, std::vector<MCHit>> mchits;
  std::map<unsigned long, std::vector<MCLAPPDHit>>* mclappdhits=nullptr;
  std::vector<MCParticle>* MCParticles=nullptr;
  int total_hits_pmts;
  int total_hits_lappds;
  int num_lappds_hit;
  int lappds_selected;
  std::string lappds_file;
  std::vector<int> active_lappds;
  int act_lappds[200]={0};
  double maximum_lappds;
  const double refractive_index_water = 1.34;
  std::ofstream _csvfile;
  int _hit_number;
  };


#endif

#ifndef MonitorDataQuality_H
#define MonitorDataQuality_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TFile.h"
#include "TTree.h"

/**
 * \class MonitorDataQuality
 *
 * This tool is meant to extract metrics from real ANNIE data per event and per run in order to compare the detector's performance over long time periods.
 * There will be a script (maybe a tool) to do the actual plots. This tool should only write properties to file that are really needed.
 * Many other variables are already produced by the PhaseIITreeMaker.
*
* $Author: M.Stender $
* $Date: 2024/08/02 19:37:00 $
* Contact: malte.stender@desy.de
*/
class MonitorDataQuality: public Tool {


 public:
  // Simple constructor
  MonitorDataQuality();
  // Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Initialise(std::string configfile,DataModel &data);
  // Function used to perform Tool purpose.
  bool Execute();
  // Function used to clean up resources.
  bool Finalise();

 private:
  // Controls the output frequency of the tool
  int m_verbosity{0};
  // Name of the .root file to which the metrics are written
  std::string m_output_file_name{"Default"};
  // File for the output
  TFile* m_output_file{nullptr};
  // Trees for holding the metrics
  TTree* m_data_quality_tree_facc{nullptr};
  TTree* m_data_quality_tree_tank_pmt{nullptr};
  TTree* m_data_quality_tree_tank_lappd{nullptr};
  TTree* m_data_quality_tree_mrd{nullptr};
  TTree* m_data_quality_tree_trigger{nullptr};

  ///========================Variables for tank PMTs========================
  // sum of charge in tank per event over all cluster, presumably in coulumb
  double m_sum_cluster_charge_tank_pmt{-999.9};
  // sum of charge in tank per event over all cluster, in photoelectrons
  double m_sum_cluster_charge_tank_pmt_pe{-999.9};
  // sum of charge in tank per event regardless of clustering presumably in coulumb
  double m_sum_charge_tank_pmt{-999.9};
  // sum of charge in tank per event regardless of clustering in photoelectrons
  double m_sum_charge_tank_pmt_pe{-999.9};
  // ratio of cluster charge to regular charge
  double ratio_charge_cluster_vs_all{-999.9};
  // ratio of cluster charge to regular charge in photoelectrons
  double ratio_charge_cluster_vs_all_pe{-999.9};

  // Extract the necessary information from the maps to store them into the corresponding variable
  void ExtractMetricsTankHits(std::map<unsigned long, std::vector<Hit>>* t_regularHitMap, std::map<double,std::vector<Hit>>* t_clusterHitMap);
};


#endif

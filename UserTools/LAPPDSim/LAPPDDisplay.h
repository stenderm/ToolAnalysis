/*
 * LAPPDDisplay.h
 *
 *  Created on: Aril 25, 2019
 *      Author: stenderm
 */

#ifndef SRC_LAPPDDisplay_H_
#define SRC_LAPPDDisplay_H_

class LAPPDDisplay{
public:
  LAPPDDisplay(string filePath, int confignumber);
  ~LAPPDDisplay();
  void InitialiseHistoAllLAPPDs(int eventNumber);
  void MCTruthDrawing(int eventnumber, unsigned long actualTubeNo, vector <MCLAPPDHit> mchits);
  void FinaliseHistoAllLAPPDs();
  void RecoDrawing(int eventCounter, unsigned long tubeNumber, std::vector<Waveform<double>> waveformVector);
private:
  TApplication* _LAPPD_sim_app = nullptr;
  TCanvas* _LAPPD_sim_MC_all = nullptr;
  TCanvas* _LAPPD_sim_MC = nullptr;
  TCanvas* _LAPPD_sim_MC_time = nullptr;
  TCanvas* _LAPPD_sim_MC_waveform = nullptr;
  TH2D* _all_hits = nullptr;
  TFile* _output_file = nullptr;
  TH2D* _left = nullptr;
  TH2D* _right = nullptr;
  int _config_number = 0;
};


#endif /* SRC_LAPPDDisplay_H_ */

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
  void OpenNewFile(int filenumber);
  void MCTruthDrawing(int eventnumber, unsigned long actualTubeNo, vector <MCLAPPDHit> mchits);
  void FinaliseHistoAllLAPPDs();
  void RecoDrawing(int eventCounter, unsigned long tubeNumber, std::vector<Waveform<double>> waveformVector);
private:
  TApplication* _LAPPD_sim_app;
  TCanvas* _LAPPD_MC_all_canvas;
  TCanvas* _LAPPD_MC_canvas;
  TCanvas* _LAPPD_MC_time_canvas;
  TCanvas* _LAPPD_all_waveforms_canvas;
  TCanvas* _LAPPD_waveform_canvas;
  TH2D* _all_hits;
  TFile* _output_file;
  int _config_number;
  string _output_file_name;
};


#endif /* SRC_LAPPDDisplay_H_ */

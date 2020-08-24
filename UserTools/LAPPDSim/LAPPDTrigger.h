/*
 * LAPPDTrigger.h
 *
 *  Created on: May 20, 2020
 *      Author: stenderm
 */
#ifndef SRC_LAPPDTrigger_H_
#define SRC_LAPPDTrigger_H_

#include <vector>
#include "TF1.h"
#include "TRandom.h"
#include "Waveform.h"

class LAPPDTrigger{
public:
  LAPPDTrigger(int numberOfLAPPDs, int numberOfChannels, double threshold, int numberOfAdjacentTriggers);
  ~LAPPDTrigger();
  std::vector< std::vector< std::vector<Waveform<double>> > > TriggerWaveforms(std::vector<Waveform<double> >, int LAPPDNumber, double startTracingTime);
  void CreateSyncFunc(double startTime, double endTime, int numberOfSamples, int SampleSize, bool IsReference);
  inline Waveform<double> GetSyncFunc(){return _sync_waveform;}
  inline Waveform<double> GetSyncFuncReference(){return _sync_waveform_reference;}

private:
  Waveform<double> ResampleWaveforms(Waveform<double> untriggeredWaveforms, int LAPPDNumber, int channelNumber, int firstSample);
  std::vector<std::vector<std::vector<double> > > * _LAPPD_sample_times;
  std::vector<std::vector<double > > * _LAPPD_channel_lengths;
  Waveform<double> _sync_waveform;
  Waveform<double> _sync_waveform_reference;
  double _trigger_threshold;
  int _number_adjacent_samples;
  int _starting_sample;
  double UniformDeviate(int seed);
};

#endif /* SRC_LAPPDTrigger_H_ */

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

private:
  std::vector<std::vector<std::vector<double> > > * _LAPPD_sample_times;
  double _trigger_threshold;
  int _number_adjacent_samples;
  int _starting_sample;
  double UniformDeviate(int seed);
};

#endif /* SRC_LAPPDTrigger_H_ */

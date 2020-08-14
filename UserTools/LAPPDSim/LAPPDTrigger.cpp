/*
* LAPPDTrigger.cpp
* This class models the properties of the LAPPD trigger. It consists of the constructor and a method. In the following a short description of the concept of this class is given:
* In order to simulate a realistic LAPPD waveform three effects need to be considered:
* 1. The sample size of a waveform is not constant, but consistent. This means for each channel we have 256 samples with different sizes, but these sizes stay constant for each channel.
* 2. There is a so called wraparound delay, when reading out the LAPPD, which occurs getting from sample 255 to sample 0. This delay is aroung 560 ps long or about 7 samples and results in a shift. Sample 0 contains sample 7, sample 1 contains sample 8 and so on.
* 3. At the so called trigger location a current injection happens, which leads to high voltages in the channels. After a time of about 20 ns the trigger signal arrives and the LAPPD is read out. This leads to high noise voltages 20 ns after the edge, which exceeded the threshold occured.
*
*
* These three effects are simulated in the following way:
* First of all, in the constructor the sample sizes are pulled of a Gaussian distribution for each sample of each channel of each LAPPD and pushed into a vector, which can be later used to have the waveform sampled in a realistic way.
* Furthermore, a number is pulled from a flat distribution representing the number of the sample at the time 0. This means, all LAPPDs have at the start time the same sample number for example 100 and then the cycle begins from 100 to 255 to 0 to 100 and so on.
*
* The TriggerWaveforms method starts than with an implementation to fill the waveform coming from the LAPPDresponse class, which uses the same constant sample size for all samples in the waveforms with different sample sizes. With this procedure the first effect is handled.
* This is followed by an implementation, which looks for voltages on the channels above a threshold, which can be set in the config file. If a number of adjacent strips exceed the threshold at the same timem, a trigger is identified. (Number of adjacent triggers is also adjustable in the config file.)
* This triggering mechanism is done for both LAPPD sides resulting in two vectors of triggers for each LAPPD. The then following code is also done for each LAPPD side individually.
* The next step is to loop over all triggers and to look at the trigger location. This is done by pulling a trigger time from a Gaussian distribution. The sample, which is the one corresponding to the trigger time is here called trigger sample.
* At this sample and the following few (this number is also pulled from a Gaussian distribution) random high voltages are injected.
* At last then the wraparound delay is implemented. The delay is also taken from a Gaussian distribution and all samples are shifted according to the delay. This means the sample 0 gets the content of sample 7 and so on.
*
*  Created on: May 20, 2020
*      Author: stenderm
*/

#include "LAPPDTrigger.h"

/**
* Constructor LAPPDTrigger: Initialises the trigger Class. Is used to create individual sample times for each sample on each channel and save them, so that ever channel has 256 samples with individual sizes,
*                           that stay constant for that channel.
* @param numberOfLAPPDs            (int) The number of LAPPDs. (For realistic data this would be 5.)
* @param numberOfChannels          (int) Number of Channels for each LAPPD. (Currently this is 60, but most likely will change to 56 for including 4 synchronisation channels per LAPPD.)
* @param threshold                 (double) Threshold in mV, which must be exceeded to trigger a LAPPD channel.
* @param numberOfAdjacentTriggers  (int) Number of simultaneous trigger on adjacent strips to trigger the LAPPD.
*/
LAPPDTrigger::LAPPDTrigger(int numberOfLAPPDs, int numberOfChannels, double threshold, int numberOfAdjacentTriggers):_LAPPD_sample_times(nullptr), _LAPPD_channel_lengths(nullptr), _trigger_threshold(threshold), _number_adjacent_samples(numberOfAdjacentTriggers), _starting_sample(0)
{
  //Create a distribution for the sample times. Sigma is 1.27 ps, mean is 97.15 ps.
  //ToDo: Check the sigma.
  TF1 *sampleTimesGaussian = new TF1("sampleTimesGaussian", "TMath::Gaus(x, 97.15,1.27)", 0, 200);

  //Vector of LAPPDs of vector of Channels of vector of 256 sample sizes
  _LAPPD_sample_times = new std::vector<std::vector<std::vector<double> > >;

  //Vector of LAPPDs of vector of Channels of the length of the waveform consisting of the 256 samples
  _LAPPD_channel_lengths = new std::vector<std::vector<double> >;

  //Loop over all LAPPDs
  for(int i = 0; i < numberOfLAPPDs; i++){
    //Vector of channels of vector of 256 sample sizes
    std::vector<std::vector<double>> channelWithSamples;
    //Loop over all Channels
    for(int j = 0; j < numberOfChannels; j++){
      //Vector of sample sizes
      std::vector<double> sampleTimes;
      //Loop over all samples (256)
      for(int k = 0; k < 256; k++){
        //Get a random sample size of the distribution defined above
        double aSampleTime = sampleTimesGaussian->GetRandom();
        sampleTimes.push_back(aSampleTime);
      }//End sample size loop
      channelWithSamples.push_back(sampleTimes);
      sampleTimes.clear();
    }//End channel loop
    _LAPPD_sample_times->push_back(channelWithSamples);
    channelWithSamples.clear();
  }//End LAPPD loop


  //This is the startSample for all Channels at the global event time 0
  //It is obtained by pulling a random sample number of a uniformly distributed distribution
  _starting_sample = UniformDeviate(rand()) * 256;

  // DEBUG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for(size_t i_LAPPD = 0; i_LAPPD < _LAPPD_sample_times->size(); i_LAPPD++){
    std::vector<double> lengths;
    for(size_t i_Channel = 0; i_Channel < _LAPPD_sample_times->at(i_LAPPD).size(); i_Channel++){
      double waveformLength = 0.0;
      for(size_t i_Sample = 0; i_Sample < _LAPPD_sample_times->at(i_LAPPD).at(i_Channel).size(); i_Sample++){
        waveformLength += _LAPPD_sample_times->at(i_LAPPD).at(i_Channel).at(i_Sample);
      }
        // std::cout << "LAPPDNumber: " << i_LAPPD << " ChannelNumber: " << i_Channel << " SampleNumber: " << waveformLength << std::endl;
        lengths.push_back(waveformLength);
    }
    _LAPPD_channel_lengths->push_back(lengths);
    lengths.clear();
  }



}


/**
* Method TriggerWaveforms:         This method is used to get triggered Waveforms. In order to do this the following steps are done.
*                                  1. Take the untriggered Waveforms and fill them in Waveforms with the right sample sizes.
*                                  2. Trigger the waveforms based on a threshold and the number of adjacent triggers.
*                                  3. Simulate the current injection at the trigger.
*                                  4.
* @param untriggeredWaveforms      (std::vector<Waveform<double> >) The untriggered Waveforms of the LAPPD with number LAPPDNumber.
* @param LAPPDNumber               (int) The number of the current LAPPD to trigger.
*/
std::vector< std::vector< std::vector<Waveform<double>> > > LAPPDTrigger::TriggerWaveforms(std::vector<Waveform<double> > untriggeredWaveforms, int LAPPDNumber, double startTracingTime){

  double originalSampleSize = 50;

  //Number of the sample at which the trigger occured
  int triggerSample = 0;
  //Number of adjacent strips, which fired simultaneously
  int adjacentCounter = 0;
  //Vector of pairs of <channel number, sample number>
  std::vector<std::pair<int, int> > vectorOfTriggerSamplesFirstSide;
  std::vector<std::pair<int, int> > vectorOfTriggerSamplesSecondSide;

  //Get the number of samples of the first of the waveform , since all waveforms have the same size.
  int numberOfSamples = untriggeredWaveforms.at(0).GetSamples()->size();

  //Loop over samples
  for (size_t i_sample = 0; i_sample < numberOfSamples; i_sample++) {
    adjacentCounter = 0;

    //Loop over channels (First side)
    for (size_t i_channel = 0; i_channel < (size_t)(untriggeredWaveforms.size()/2); i_channel++) {

      //if trigger threshold is exceeded add one to the adjacent Counter
      if(abs(untriggeredWaveforms.at(i_channel).GetSample(i_sample)) > _trigger_threshold){
        adjacentCounter++;

      //if trigger threshold is not exceeded, check, whether sufficient adjacent strips were triggered

        if(adjacentCounter >= _number_adjacent_samples){
          //Save at which sample the trigger fired at which channel
          triggerSample = i_sample;
          vectorOfTriggerSamplesFirstSide.push_back(std::pair<int,int>(i_channel,i_sample));
          //Since 256 samples will be read out, no trigger is possible for the next 256 samples, so skip these
          //Since we are currently looking at samples of the size of 50 ps, we need to skip about 500
          //ToDo: Finetune this number?
          i_sample += 500;
          //reset the adjacentCounter
          adjacentCounter = 0;
          //Leave the channel loop for this sample and the next 256 samples, since a trigger is found
          break;
        }
      }
    }//end of channel loop
  }//end of sample loop

  //Loop over samples
  for (size_t i_sample = 0; i_sample < numberOfSamples; i_sample++) {
    adjacentCounter = 0;
    //Loop over channels (Second side)
    for (size_t i_channel = (size_t)(untriggeredWaveforms.size()/2); i_channel < untriggeredWaveforms.size(); i_channel++) {
      //if trigger threshold is exceeded add one to the adjacent Counter
      if(abs(untriggeredWaveforms.at(i_channel).GetSample(i_sample)) > _trigger_threshold){
        adjacentCounter++;

      //if trigger threshold is not exceeded, check, whether sufficient adjacent strips were triggered
      if(adjacentCounter >= _number_adjacent_samples){
          //Save at which sample the trigger fired at which channel
          triggerSample = i_sample;
          vectorOfTriggerSamplesSecondSide.push_back(std::pair<int,int>(i_channel,i_sample));
          //Since 256 samples will be read out, no trigger is possible for the next 256 samples, so skip these
          //Since we are currently looking at samples of the size of 50 ps, we need to skip about 500
          //ToDo: Finetune this number?
          i_sample += 500;
          //reset the adjacentCounter
          adjacentCounter = 0;
          //Leave the channel loop for this sample and the next 256 samples, since a trigger is found
          break;
        }
      }
    }//end of channel loop
  }//end of sample loop

  //Vector of Channels of vector of channels of vector of samples
  std::vector<std::vector<std::pair<int, vector<double> > > > triggeredWaveformsLeft;
  std::vector<std::vector<std::pair<int, vector<double> > > > triggeredWaveformsRight;
  std::vector<std::pair<int, vector<double> > > oneTriggeredWaveform;

  //ToDo: Check Sigmas
  //Distribution for the time the LAPPD needs to start the readout after a trigger occured
  //Mean: 20 ns; Sigma: 2.5 ns
  TF1 *triggerTimeGaussian = new TF1("triggerTimeGaussian", "TMath::Gaus(x, 20.0, 2.5)", 0, 100);

  //Distribution of the number of samples, which get the curent injection
  //Originally, I used Mean: 8.5; Sigma: 1 for the already sampled waveform
  //Now, I use Mean: 16.5; Sigma: 1 for the constant sized waveforms
  TF1 *triggerSampleSize = new TF1("triggerSampleSize", "TMath::Gaus(x, 8.5, 1)", 0, 100);

  //Distribution of the voltage of the current injection
  //Mean: 80 mV; Sigma: 20 mV
  TF1 *triggerVoltage = new TF1("triggerSampleSize", "TMath::Gaus(x, 80, 20)", 0, 100);

  //Distribution of the wraparound delay
  //Mean: 560 ps; Sigma: 20 ps
  TF1 *wraparoundDelay = new TF1("wraparoundDelay", "TMath::Gaus(x, 560, 20)", 0, 1000);

  //ToDo: Ich muss nun noch für JEDEN Channel die Current Injection machen und mir Gedanken um den Wraparound machen.
  //Insbesondere wird es schwierig, die zweite Seite abzuspreichern, da ich ja nur für eine Seite triggere. Ggf. rechne ich
  //mithilfe der Trigger-Zeit einfach aus, wann ich die andere Seite abzuspeichern habe. Ich werde einfach beide Seiten
  //triggern und dann, wenn das alles fertig ist, nach Evans Ansagen ändern.

  std::vector< std::vector<Waveform<double>> > resultWaveformsAllTriggerLeftSide;
  std::vector<Waveform<double>> resultWaveformOneTriggerLeftSide;

  //Loop over the trigger on one side
  for (size_t i_trigger = 0; i_trigger < vectorOfTriggerSamplesFirstSide.size(); i_trigger++) {
    double sampleSize = 0;

    // for(size_t i_sample = 0; i_sample <= vectorOfTriggerSamplesFirstSide.at(i_trigger).second; i_sample++){
    //   sampleSize += sampledWaveforms.at(i_trigger).at(i_sample).second.at(0);
    // }

    for(size_t i_channel = 0; i_channel < untriggeredWaveforms.size()/2; i_channel++){

      //Get values
      int triggerSample = vectorOfTriggerSamplesFirstSide.at(i_trigger).second;
      int triggerChannel = vectorOfTriggerSamplesFirstSide.at(i_trigger).first;


      //Get a time for the readout
      double triggerTime = triggerTimeGaussian->GetRandom();
      //Time that is already gone since the trigger occured
      double timeGone = 0;
      //sample that is currently looked at
      int sampleLoop = triggerSample;

      //while-loop over the samples as long as the added up sample times to not exceed the trigger readout time

      while(timeGone < triggerTime){
        timeGone += originalSampleSize/1000;
        sampleLoop++;
      }

      //Get the number of samples that get a current injection
      int numberOfNoiseSamples = (int) triggerSampleSize->GetRandom();

      //loop over noise samples
      for(int i = 0; i <= numberOfNoiseSamples; i++){
        untriggeredWaveforms.at(i_channel).SetSample(sampleLoop-i, -triggerVoltage->GetRandom());
      }

      //Implementation of the wraparound effect
      double wraparoundTime = wraparoundDelay->GetRandom();
      double sumSampleTimes = 0;
      int numberOfShiftedSamples = 0;

      //Calculate start time
      double waveformStartTime = 0.0;
      Waveform<double> waveformVoltages;
      double lengthOfWaveform = _LAPPD_channel_lengths->at(LAPPDNumber).at(i_channel);
      double numberOfSamplesDouble = lengthOfWaveform/originalSampleSize;
      int numberOfSamplesInt = (int)numberOfSamplesDouble + 1;

      //ToDo: Find a way to get rid of the negative numbers here
      for(size_t i_sample = 0 ; i_sample < sampleLoop - numberOfSamplesInt; i_sample++){
        // std::cout << "i_sample " << i_sample << std::endl;
        waveformStartTime += originalSampleSize;
      }
      waveformStartTime = waveformStartTime + startTracingTime;

      //Loop over the samples to get a waveform with 256 samples
      int firstSample = (_starting_sample + sampleLoop - numberOfSamplesInt)%256;

      for(size_t i_sample = sampleLoop - numberOfSamplesInt; i_sample <= sampleLoop; i_sample++){
        //Look for the 0th sample, since at that point the wraparound delay needs to be applied.

        int currentSample = (_starting_sample + i_sample)%256;
        if(currentSample == 0){
         //Calculate the number of samples, which are overwritten due to the delay

          while(wraparoundTime > sumSampleTimes){
            sumSampleTimes += originalSampleSize;
            numberOfShiftedSamples++;
          }

        }
        //Create a waveform and fill the samples.
        waveformVoltages.PushSample(untriggeredWaveforms.at(i_channel).GetSample(i_sample + numberOfShiftedSamples));
      }
      //Set the start time for the waveform and push it to a vector for the one trigger currently looked at
      waveformVoltages.SetStartTime(waveformStartTime);
      waveformVoltages = ResampleWaveforms(waveformVoltages, LAPPDNumber, i_channel, firstSample);
      resultWaveformOneTriggerLeftSide.push_back(waveformVoltages);

    }
    //Save the waveforms for all triggers for the left side
    resultWaveformsAllTriggerLeftSide.push_back(resultWaveformOneTriggerLeftSide);
    resultWaveformOneTriggerLeftSide.clear();
  }


  std::vector< std::vector<Waveform<double>> > resultWaveformsAllTriggerRightSide;
  std::vector<Waveform<double>> resultWaveformOneTriggerRightSide;


  //Loop over the trigger other side
  //This is basically copied code from above.
  for (size_t i_trigger = 0; i_trigger < vectorOfTriggerSamplesSecondSide.size(); i_trigger++) {
    double sampleSize = 0;

    // for(size_t i_sample = 0; i_sample <= vectorOfTriggerSamplesFirstSide.at(i_trigger).second; i_sample++){
    //   sampleSize += sampledWaveforms.at(i_trigger).at(i_sample).second.at(0);
    // }

    for(size_t i_channel = untriggeredWaveforms.size()/2; i_channel < untriggeredWaveforms.size(); i_channel++){

      //Get values
      int triggerSample = vectorOfTriggerSamplesSecondSide.at(i_trigger).second;
      int triggerChannel = vectorOfTriggerSamplesSecondSide.at(i_trigger).first;


      //Get a time for the readout
      double triggerTime = triggerTimeGaussian->GetRandom();
      //Time that is already gone since the trigger occured
      double timeGone = 0;
      //sample that is currently looked at
      int sampleLoop = triggerSample;

      //while-loop over the samples as long as the added up sample times do not exceed the trigger readout time
      while(timeGone < triggerTime){
        timeGone += originalSampleSize/1000;
        sampleLoop++;
      }

      //Get the number of samples that get a current injection
      int numberOfNoiseSamples = (int) triggerSampleSize->GetRandom();

      //loop over noise samples
      for(int i = 0; i <= numberOfNoiseSamples; i++){
        untriggeredWaveforms.at(i_channel).SetSample(sampleLoop-i, -triggerVoltage->GetRandom());
      }

      //Implementation of the wraparound effect
      double wraparoundTime = wraparoundDelay->GetRandom();
      double sumSampleTimes = 0;
      int numberOfShiftedSamples = 0;

      //Calculate start time
      double waveformStartTime = 0.0;
      Waveform<double> waveformVoltages;
      double lengthOfWaveform = _LAPPD_channel_lengths->at(LAPPDNumber).at(i_channel);
      double numberOfSamplesDouble = lengthOfWaveform/originalSampleSize;
      int numberOfSamplesInt = (int)numberOfSamplesDouble + 1;

      for(size_t i_sample = 0 ; i_sample < sampleLoop - numberOfSamplesInt; i_sample++){
        waveformStartTime += originalSampleSize;
      }

      waveformStartTime = waveformStartTime + startTracingTime;

      int firstSample = (_starting_sample + sampleLoop - numberOfSamplesInt)%256;

      for(size_t i_sample = sampleLoop - numberOfSamplesInt; i_sample <= sampleLoop; i_sample++){
        int currentSample = (_starting_sample + i_sample)%256;
        if(currentSample == 0){
         //Calculate the number of samples, which are overwritten due to the delay
          while(wraparoundTime > sumSampleTimes){
            sumSampleTimes += originalSampleSize;
            numberOfShiftedSamples++;
          }
        }
        //Create a waveform and fill the samples.
        waveformVoltages.PushSample(untriggeredWaveforms.at(i_channel).GetSample(i_sample + numberOfShiftedSamples));
      }//end of sample loop
      //Set the start time for the waveform and push it to a vector for the one trigger currently looked at
      waveformVoltages.SetStartTime(waveformStartTime);
      waveformVoltages = ResampleWaveforms(waveformVoltages, LAPPDNumber, i_channel, firstSample);
      resultWaveformOneTriggerRightSide.push_back(waveformVoltages);
    }//end of channel loop
    //Save the waveforms for all triggers for the left side

    resultWaveformsAllTriggerRightSide.push_back(resultWaveformOneTriggerRightSide);
    resultWaveformOneTriggerRightSide.clear();
  }//end of trigger loop



  // std::cout << "Size of Waveforms left " << triggeredWaveformsLeft.size() << std::endl;
  // std::cout << "Size of Waveforms right " << triggeredWaveformsRight.size() << std::endl;
  //
  // std::cout << "Size of Waveforms left trigger " << resultWaveformsAllTriggerLeftSide.size() << std::endl;
  // std::cout << "Size of Waveforms right trigger " << resultWaveformsAllTriggerRightSide.size() << std::endl;

std::vector< std::vector< std::vector<Waveform<double>> > > resultWaveformsAllTriggerBothSides;
resultWaveformsAllTriggerBothSides.push_back(resultWaveformsAllTriggerLeftSide);
resultWaveformsAllTriggerBothSides.push_back(resultWaveformsAllTriggerRightSide);
return resultWaveformsAllTriggerBothSides;
}


Waveform<double> LAPPDTrigger::ResampleWaveforms(Waveform<double> untriggeredWaveforms, int LAPPDNumber, int channelNumber, int firstSample){

    Waveform<double> resampledWaveform;
    resampledWaveform.SetStartTime(untriggeredWaveforms.GetStartTime());

    //Vector of channels of vector of samples of pairs corresponding to the sample number and a vector of two doubles corresponding to the sample size and the sample voltage
    std::vector<std::vector<std::pair<int, vector<double> > > > sampledWaveforms;

      //Time of the start of the sample currently looked at in the sampled waveforms
      double startSampleTime = 0.0;
      //Voltage to fill in the sample currently looked at
      double voltageToFill = 0;

      //Number of the sample in the sampled waveforms
      int sampleCounter = firstSample;

      //Get the random but constant sample size, which was created in the constructor for this particular LAPPD at that particular channel for the starting sample
      double sampleSize = _LAPPD_sample_times->at(LAPPDNumber).at(channelNumber).at(firstSample%256);

      //Time of the end of the sample currently looked at in the samples waveforms
      double endSampleTime = startSampleTime + sampleSize;

      //Vector of all samples of one waveform, with <sample size, <sample size, sample voltage> >
      std::vector<std::pair<int,vector<double> > > sampledVoltages;

      //Loop over samples
      for(size_t i_Sample = 0; i_Sample < untriggeredWaveforms.GetSamples()->size(); i_Sample++){
        // std::cout << "ChannelNumber " << i_Channel << "SampleNumber " << i_Sample << " Voltage " << untriggeredWaveforms.at(i_Channel).GetSample(i_Sample) << std::endl;

        //Original sample size of the untriggered waveforms used for the GetTrace()-Method of the LAPPDresponse class
        //ToDo: Implement this as a parameter for this function
        double originalSize = 50;

        //Get the voltage of the untriggered waveform at the sample i_Sample
        double voltage = untriggeredWaveforms.GetSample(i_Sample);

        //Calculate the start time of the sample of the untriggered waveform
        double startTemplateTime = i_Sample * originalSize;

        //Calculate the end time of the sample of the untriggered waveform
        double endTemplateTime = (i_Sample + 1) * originalSize;

        //If the end time of the sample of the sampled waveform is bigger than the end time of the sample of the untriggered waveforms,
        //then the percentage of voltage from the untriggered waveforms can be filled in the sample of the sampled waveforms, since
        //the constant sized sample is at an end
        // |....|....|....|....| <-constant samples
        // |vvvv|
        // |.......|........|......| <- samples with different size
        //
        if(endSampleTime > endTemplateTime){
          //Calculate what percentage the constant size sample makes up in the sample with different size and add this to the voltage of the non-constant sample
          voltageToFill += ((endTemplateTime - startSampleTime)/sampleSize)*voltage;

          //When in last iteration the values need to be pushed into the vector
          if(i_Sample == untriggeredWaveforms.GetSamples()->size()-1){
            //vector of sample size and sample voltage
            std::vector<double> sizeVoltage;
            //Fill the size of the sample in the vector
            sizeVoltage.push_back(sampleSize);
            //Fill the calculated voltage in the vector
            sizeVoltage.push_back(voltageToFill);
            //Make a pair of the sample number an the vector above and push it in the sample vector
            sampledVoltages.push_back( std::pair<int,std::vector<double> >(sampleCounter%256, sizeVoltage));
            resampledWaveform.PushSample(voltageToFill);
          }
        }
        //Now the end time of the constant sample is bigger than the non-constant one, so that the calculation of
        //the voltage of the non-constant sample is finished
        // |....|....|....|....| <-constant samples
        //      |vv|
        // |.......|........|......| <- samples with different size
        else{
          //Caclulate what percentage the constant size sample makes up in the sample with different size and add this to the voltage of the non-constant sample
          voltageToFill += ((endSampleTime - startTemplateTime)/sampleSize)*voltage;
          //vector of sample size and sample voltage
          std::vector<double> sizeVoltage;
          //Fill the size of the sample in the vector
          sizeVoltage.push_back(sampleSize);
          //Fill the calculated voltage in the vector
          sizeVoltage.push_back(voltageToFill);
          //Make a pair of the sample number an the vector above and push it in the sample vector
          sampledVoltages.push_back( std::pair<int,std::vector<double> >(sampleCounter%256, sizeVoltage));
          resampledWaveform.PushSample(voltageToFill);
          //Reset values
          sizeVoltage.clear();
          voltageToFill = 0;
          //Get the current sample size
          sampleSize = _LAPPD_sample_times->at(LAPPDNumber).at(channelNumber).at(sampleCounter%256);
          //Calculate the new sample start time
          startSampleTime += sampleSize;
          //Go to the next sample in the non-constant sample waveform
          sampleCounter++;
          //Get the sample size of the next sample
          sampleSize = _LAPPD_sample_times->at(LAPPDNumber).at(channelNumber).at(sampleCounter%256);
          //Calculate the new sample end time
          endSampleTime = startSampleTime + sampleSize;
        }
      }//end of sample loop

      sampledWaveforms.push_back(sampledVoltages);
      sampledVoltages.clear();

      return resampledWaveform;


}





  /**
  * Method UniformDeviate: This method is used to create an uniformly distributed distribution based on a seed.
  * @param seed            (int) Seed used for the creation of the distribution.
  * @return                (double) A random number out of the distribution
  */
  double LAPPDTrigger::UniformDeviate(int seed){
    return seed * (1.0 / (RAND_MAX + 1.0));
  }

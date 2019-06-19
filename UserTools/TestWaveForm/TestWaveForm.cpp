#include "TestWaveForm.h"

TestWaveForm::TestWaveForm():Tool(){}


bool TestWaveForm::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////




  return true;
}


bool TestWaveForm::Execute(){
  bool testgeom = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", _geom);
  if (not testgeom)
  {
    std::cerr << "LAPPDSim Tool: Could not find Geometry in the ANNIEEvent!" << std::endl;
    return false;
  }

  bool testwav = m_data->Stores["ANNIEEvent"]->Get("LAPPDWaveforms", _waveforms);
  if (not testwav)
  {
    std::cerr << "TestWaveForm Tool: Could not find LAPPDWaveforms in the ANNIEEvent!" << std::endl;
    return false;
  }
   std::map<unsigned long, Waveform<double> >::iterator itr;
   for(itr = _waveforms->begin(); itr != _waveforms->end(); itr++){
     std::cout << "Channelkey " << itr->first << std::endl;
     Waveform<double> Penis = itr->second;
     std::vector<double> * samples = Penis.GetSamples();
     double time = Penis.GetStartTime();
     std::cout << "Start time " << time << std::endl;
     for(int i = 0; i < samples->size(); i++){
       std::cout << "Sample " << i << " Voltage " << samples->at(i) << std::endl;
     }
   }


  return true;
}


bool TestWaveForm::Finalise(){

  return true;
}

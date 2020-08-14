#include "TestLAPPDSimWaveforms.h"
#include "Waveform.h"
#include "Channel.h"
#include "Detector.h"

TestLAPPDSimWaveforms::TestLAPPDSimWaveforms()
{
  //Do nothing
}

bool TestLAPPDSimWaveforms::Initialise(std::string configfile, DataModel &data)
{
  std::cout << "Initialise TestLAPPDSimWaveforms" << std::endl;
  /////////////////// Usefull header ///////////////////////
  if (configfile != "")
  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();
  m_data = &data; //assigning transient data pointer
  return true;
}

bool TestLAPPDSimWaveforms::Execute()
{
  std::cout << "Execute TestLAPPDSimWaveforms" << std::endl;
  Geometry* AnnieGeometry;
  //Get the Geometry information
  bool testgeom = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", AnnieGeometry);
  if (not testgeom)
  {
    std::cerr << "TestLAPPDSimWaveforms Tool: Could not find Geometry in the ANNIEEvent!" << std::endl;
    return false;
  }

  std::map<unsigned long, Waveform<double> >* TriggeredLAPPDWaveforms;
  bool testwave = m_data->Stores["ANNIEEvent"]->Get("TriggeredLAPPDWaveforms", TriggeredLAPPDWaveforms);
  if(not testwave)
  {
    std::cerr << "TestLAPPDSimWaveforms Tool: No triggered Waveforms found for this event!" << std::endl;
    return false;
  }

  std::map<unsigned long, Waveform<double> >::iterator it;
  for(it = TriggeredLAPPDWaveforms->begin(); it != TriggeredLAPPDWaveforms->end(); ++it){
    unsigned long channelnumber = it->first;
    Detector* theLAPPD = AnnieGeometry->ChannelToDetector(channelnumber);
    unsigned long tubeNumber = theLAPPD->GetDetectorID();
    Channel* theChannel = AnnieGeometry->GetChannel(channelnumber);

    Waveform<double> aWaveform = it->second;
    // std::cout << "There was a hit at LAPPD number " << tubeNumber << std::endl;
    // std::cout << "Channelkey " << channelnumber << " represents channelID " << theChannel->GetChannelID() << " with strip number " << theChannel->GetStripNum() << " on side " << theChannel->GetStripSide() << std::endl;
  }


  return true;
}

bool TestLAPPDSimWaveforms::Finalise(){
  return true;
}

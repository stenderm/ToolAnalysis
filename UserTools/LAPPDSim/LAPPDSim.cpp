#include "LAPPDSim.h"
#include <unistd.h>

LAPPDSim::LAPPDSim() :
		Tool()
{
}

bool LAPPDSim::Initialise(std::string configfile, DataModel &data)
{

	/////////////////// Usefull header ///////////////////////
	if (configfile != "")
		m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	m_data = &data; //assigning transient data pointer

	m_variables.Get("EventDisplay", _display_config);
	std::cout << "DisplayNumber " << _display_config << std::endl;
	std::string outputFile;
	m_variables.Get("OutputFile", outputFile);
	std::cout << "OutputFile " << outputFile << std::endl;

	//Get the Geometry information
	bool testgeom = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", _geom);
	if (not testgeom)
	{
		std::cerr << "LAPPDSim Tool: Could not find Geometry in the ANNIEEvent!" << std::endl;
		return false;
	}
	/////////////////////////////////////////////////////////////////

	// here I would open the input file

	// here I would also store relevant geometry information

	//m_data->Stores["ANNIEEvent"] = new BoostStore(false, 2);

// This quantity should be set to false if we are working with real data later
	//bool isSim = true;
	//m_data->Stores["ANNIEEvent"]->Header->Set("isSim",isSim);

	// initialize the ROOT random number generator
	myTR = new TRandom3();

	m_variables.Get("SimInput", SimInput);

	iter = 0;
	_tf = new TFile("/nashome/m/mstender/WorkingToolAnalysis/UserTools/LAPPDSim/pulsecharacteristics.root", "READ");
	if (_display_config > 0)
	{
		_display = new LAPPDDisplay(outputFile, _display_config);
	}
	return true;
}

bool LAPPDSim::Execute()
{
	std::cout << "Executing LAPPDSim; event counter " << _event_counter << std::endl;

	//std::map<unsigned long, vector<Waveform<double>>> RawLAPPDData;
	LAPPDWaveforms = new std::map<unsigned long, Waveform<double> >;
	LAPPDWaveforms->clear();

	if (iter % 100 == 0)
		cout << "iteration: " << iter << endl;

	// get the MC Hits
	std::map<unsigned long, std::vector<MCLAPPDHit> >* lappdmchits;
	bool testval = m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHits", lappdmchits);

	if (not testval)
	{
		std::cerr << "LAPPDSim Tool: Could not find MCLAPPDHits in the ANNIEEvent!" << std::endl;
		return false;
	}

	std::map<unsigned long, std::vector<MCLAPPDHit> >::iterator itr;
	//Draw the histograms for the MC truth information
	if (_display_config > 0)
	{
		_display->InitialiseHistoAllLAPPDs(_event_counter);
	}

	// loop over the number of lappds
	for (itr = lappdmchits->begin(); itr != lappdmchits->end(); ++itr)
	{

		unsigned long tubeno = itr->first;
		Detector* thelappd = _geom->ChannelToDetector(tubeno);
		unsigned long actualTubeNo = thelappd->GetDetectorID();
		vector<MCLAPPDHit> mchits = itr->second;
		if (_display_config > 0)
		{
			_display->MCTruthDrawing(_event_counter, actualTubeNo, mchits);
		}

		std::vector<double> pulsetimes;
		//LAPPDresponse* response = new LAPPDresponse();  //SD
		LAPPDresponse response;
		response.Initialise(_tf);

		//loop over the pulses on each lappd
		for (int j = 0; j < mchits.size(); j++)
		{
			// Here we would input these pulses into our lappd model
			// and extract the signals on each of 60 channels...
			// For now we just extract the 2 Tpsec times
			// and input them in 5 channels

			LAPPDHit ahit = mchits.at(j);
			double atime = ahit.GetTime();      //*1000.;
			pulsetimes.push_back(atime);
			vector<double> localpos = ahit.GetLocalPosition();  //SD
			double trans = localpos.at(1) * 1000;         //SD
			double para = localpos.at(0) * 1000;               //SD
			response.AddSinglePhotonTrace(trans, para, atime);       //SD
		}
		vector<Waveform<double>> Vwavs;

		for (int i = -30; i < 31; i++)
		{
			if (i == 0)
			{
				continue;
			}
			Waveform<double> awav = response.GetTrace(i, 0.0, 100, 256, 1.0);
			Vwavs.push_back(awav);
			//RawLAPPDData.insert(pair<unsigned long, vector<Waveform<double>>>(actualTubeNo, Vwavs));
		}

		std::map<unsigned long, Channel>* lappdchannel = thelappd->GetChannels();
		int numberOfLAPPDChannels = lappdchannel->size();
		std::map<unsigned long, Channel>::iterator chitr;
		//std::cout << "TubeNumber " << tubeno << std::endl;
		for (chitr = lappdchannel->begin(); chitr != lappdchannel->end(); ++chitr)
		{
			//	std::cout << "TubeNo " << actualTubeNo << std::endl;
			//std::cout << "Channelkey " << chitr->first << std::endl;
			Channel achannel = chitr->second;
//			    std::cout << "Channel ID " << achannel.GetChannelID() << std::endl;
//			    std::cout << "Stripnumber " << achannel.GetStripNum() << std::endl;
//			    std::cout << "Stripside " << achannel.GetStripSide() << std::endl;
			if (achannel.GetStripSide() == 0)
			{
				LAPPDWaveforms->insert(pair<unsigned long, Waveform<double>>(achannel.GetChannelID(), Vwavs[achannel.GetStripNum()]));
			} else
			{
				LAPPDWaveforms->insert(pair<unsigned long, Waveform<double>>(achannel.GetChannelID(), Vwavs[numberOfLAPPDChannels - achannel.GetStripNum() - 1]));
			}

		}

		if (_display_config > 0)
		{
			_display->RecoDrawing(_event_counter, actualTubeNo, Vwavs);

			if (_display_config == 2)
			{
				do
				{
					std::cout << "Press a key to continue..." << std::endl;
				} while (cin.get() != '\n');

				std::cout << "Continuing" << std::endl;
			}
		}
	}				//end loop over LAPPDs

	if (_display_config > 0)
	{
		_display->FinaliseHistoAllLAPPDs();
	}
	// std::map<unsigned long, Waveform<double> >::iterator itera;
	// for(itera = LAPPDWaveforms->begin(); itera != LAPPDWaveforms->end(); itera++){
	// 	std::cout << "Channelkey " << itera->first << std::endl;
	// 	Waveform<double> Penis = itera->second;
	// 	std::vector<double> * samples = Penis.GetSamples();
	// 	for(int i = 0; i < samples->size(); i++){
	// 		std::cout << "Sample " << i << " Voltage " << samples->at(i) << std::endl;
	// 	}
	// }


	std::cout << "Saving waveforms to store" << std::endl;

	m_data->Stores.at("ANNIEEvent")->Set("LAPPDWaveforms", LAPPDWaveforms, true);
	//m_data->Stores["ANNIEEvent"]->Set("RawLAPPDData", RawLAPPDData);

	iter++;
	_event_counter++;
	return true;
}

bool LAPPDSim::Finalise()
{
	_tf->Close();
	_display->~LAPPDDisplay();
	return true;
}

Waveform<double> LAPPDSim::SimpleGenPulse(vector<double> pulsetimes)
{

	int npulses = pulsetimes.size();

	// generate gaussian TF1s for each pulse
	TF1** aGauss = new TF1*[npulses];
	for (int i = 0; i < npulses; i++)
	{
		TString gname;
		gname += "gaus";
		gname += i;
		aGauss[i] = new TF1(gname, "gaus", 0, 256);

		// random pulse amplitude chosen with from a gaussian distribution
		double theamp = fabs(myTR->Gaus(30., 30.)); // 30 mV mean, 30 mV sigma

		aGauss[i]->SetParameter(0, theamp); // amplitude of the pulse
		aGauss[i]->SetParameter(1, pulsetimes.at(i)); // peak location (in samples)
		aGauss[i]->SetParameter(2, 8.); // width (sigma) of the pulse
	}

	// loop over 256 samples
	// generate the trace, populated with the fake pulses
	Waveform<double> thewav;
	for (int i = 0; i < 256; i++)
	{

		double noise = myTR->Gaus(0., 2.0); //add in random baseline noise (2 mV sig)
		double signal = 0;

		//now add all the pulses to the signal
		for (int j = 0; j < npulses; j++)
		{
			signal += aGauss[j]->Eval(i, 0, 0, 0);
		}

		double thevoltage = signal + noise;
		thewav.PushSample(-thevoltage);
	}

	return thewav;
}

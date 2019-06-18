/*
 * LAPPDDisplay.cpp
 *
 *  Created on: April 25, 2019
 *      Author: stenderm
 */

#include "LAPPDDisplay.h"
#include <thread>
/**
 * Constructor LAPPDDisplay: Initialises the class with a TApplication, a canvas and the output file.
 * @param filePath      Path and name of the output root file.
 */

LAPPDDisplay::LAPPDDisplay(std::string filePath, int confignumber)
{
	_config_number = confignumber;
	int myargc = 0;
	char *myargv[] =
	{ (const char*) "somestring" };
	_LAPPD_sim_app = new TApplication("LAPPDSimApp", &myargc, myargv);
	if (_config_number == 2)
	{
		std::cout << "_config_number " << _config_number << std::endl;
		Double_t canvwidth = 1280;
		Double_t canvheight = 720;
		_LAPPD_sim_MC_all = new TCanvas("LAPPDSimMCAllCanvas", "LAPPDSimMCAllCanvas", canvwidth, canvheight);
		_LAPPD_sim_MC = new TCanvas("LAPPDSimMCCanvas", "LAPPDSimMCCanvas", canvwidth, canvheight);
		_LAPPD_sim_MC->SetRightMargin(0.15);
		_LAPPD_sim_MC_time = new TCanvas("LAPPDSimMCTimeCanvas", "LAPPDSimMCTimeCanvas", canvwidth, canvheight);
		_LAPPD_sim_MC_time->Divide(2, 1, 0.01, 0.01, 0);
		_LAPPD_sim_MC_waveform = new TCanvas("LAPPDSimWaveformCanvas", "LAPPDSimWaveformCanvas", canvwidth, canvheight);
		_LAPPD_sim_MC_waveform->Divide(2, 1, 0.015, 0.01, 0);
		_LAPPD_sim_MC_waveform->GetPad(1)->SetRightMargin(0.15);
		_LAPPD_sim_MC_waveform->GetPad(2)->SetRightMargin(0.15);
	}
	const char *filename = filePath.c_str();
	std::cout << "filename " << filename << std::endl;
	_output_file = new TFile(filename, "Recreate");
}

/**
 * Destructor LAPPDDisplay: Closes the output file, deletes the canvas and the TApplication.
 */

LAPPDDisplay::~LAPPDDisplay()
{
	if (_config_number == 2)
	{
		if (gROOT->FindObject("LAPPDSimCanvas") != nullptr)
		{
			delete _LAPPD_sim_MC_time;
		}
		if (gROOT->FindObject("LAPPDSimMCAllCanvas") != nullptr)
		{
			delete _LAPPD_sim_MC_all;
		}
		if (gROOT->FindObject("LAPPDSimMCCanvas") != nullptr)
		{
			delete _LAPPD_sim_MC;
		}
		if (gROOT->FindObject("LAPPDSimWaveformCanvas") != nullptr)
		{
			delete _LAPPD_sim_MC_waveform;
		}
	}
	delete _LAPPD_sim_app;
	_output_file->Close();

}


void LAPPDDisplay::InitialiseHistoAllLAPPDs(int eventNumber){
	//Create the name for the histogram for all LAPPDs
	std::string eventnumber = boost::lexical_cast < std::string > (eventNumber);
	string allHitsName = "event" + eventnumber + "AllLAPPDs";
	const char *allHitsNamec = allHitsName.c_str();
	//Initialisation of the histogram for all LAPPDs
	_all_hits = new TH2D(allHitsNamec, allHitsNamec, 200, 0, 180, 200, -1.1, 1.1);
}

/**
 * Method MCTruthDrawing: Draws the MCtruth hits as 2D histogram with transverse coordinate as x-axis, parallel coordinate as y-axis and hit time as color code.
 * This is done for every LAPPD individually and also for all LAPPDs in an event.
 * @param eventNumber The number of the event, needed for the naming of the histograms.
 * @param lappdmchits Map of all LAPPD hits
 */

void LAPPDDisplay::MCTruthDrawing(int eventNumber, unsigned long actualTubeNo, vector <MCLAPPDHit> mchits)
{
		std::string eventnumber = boost::lexical_cast < std::string > (eventNumber);
		//Create the name for the histogram for individual LAPPDs
		std::string lappdnumber = boost::lexical_cast < std::string > (actualTubeNo);
		string individualName = "event" + eventnumber + "lappd" + lappdnumber;
		const char *individualNamec = individualName.c_str();
		//Initialisation of the histogram for individual LAPPDs
		TH2D* LAPPDMCHits = new TH2D(individualNamec, individualNamec, 60, -0.12, 0.12, 60, -0.12, 0.12);
		//loop over all MC hits
		double mintime = 1000000;
		for (int k = 0; k < mchits.size(); k++)
		{
			LAPPDHit ahit = mchits.at(k);
			double atime = ahit.GetTime();
			if (mintime > atime)
			{
				mintime = atime;
			}
		}
		std::string nameTrans = "event" + eventnumber + "lappd" + lappdnumber + "trans";
		std::string namePara = "event" + eventnumber + "lappd" + lappdnumber + "para";
		const char * charNameTrans = nameTrans.c_str();
		const char * charNamePara = namePara.c_str();
		TH2D* LAPPDTrans = new TH2D(charNameTrans, charNameTrans, 256, mintime, mintime + 25.6, 60, -0.12, 0.12);
		TH2D* LAPPDPara = new TH2D(charNamePara, charNamePara, 256, mintime, mintime + 25.6, 60, -0.12, 0.12);

		for (int i = 0; i < mchits.size(); i++)
		{
			LAPPDHit ahit = mchits.at(i);
			double atime = ahit.GetTime(); //*1000.;
			vector<double> localpos = ahit.GetLocalPosition();
			double trans = localpos.at(1);
			double para = localpos.at(0);
			LAPPDMCHits->Fill(trans, para, atime);
			LAPPDTrans->Fill(atime, trans, 1);
			LAPPDPara->Fill(atime, para, 1);
			Position globalHit(ahit.GetPosition().at(0), ahit.GetPosition().at(1), ahit.GetPosition().at(2));
			double phi = TMath::ATan2(globalHit.Z(), globalHit.X()) * (180 / TMath::Pi());
			_all_hits->Fill(phi, globalHit.Y(), atime);
		}

		LAPPDMCHits->GetXaxis()->SetTitle("Transverse coordinate [m]");
		LAPPDMCHits->GetYaxis()->SetTitle("Parallel coordinate [m]");
		LAPPDMCHits->GetZaxis()->SetTitle("Arrival time [ns]");
		LAPPDMCHits->GetZaxis()->SetTitleOffset(1.4);
		LAPPDMCHits->Write();

		LAPPDTrans->GetYaxis()->SetTitle("Transverse coordinate [m]");
		LAPPDTrans->GetYaxis()->SetTitleOffset(1.4);
		LAPPDTrans->GetXaxis()->SetTitle("Arrival time [ns]");
		LAPPDTrans->GetZaxis()->SetTitle("Events");
		LAPPDTrans->Write();

		LAPPDPara->GetYaxis()->SetTitle("Parallel coordinate [m]");
		LAPPDPara->GetYaxis()->SetTitleOffset(1.4);
		LAPPDPara->GetXaxis()->SetTitle("Arrival time [ns]");
		LAPPDPara->GetZaxis()->SetTitle("Events");
		LAPPDPara->Write();
		if(_config_number == 2)
		{
		_LAPPD_sim_MC->cd();
		LAPPDMCHits->SetStats(0);
		LAPPDMCHits->Draw("COLZ");
		_LAPPD_sim_MC->Modified();
		_LAPPD_sim_MC->Update();

		_LAPPD_sim_MC_time->cd(1);
//		LAPPDTrans->SetMarkerStyle(21);
//		LAPPDTrans->SetMarkerSize(1.5);
		LAPPDTrans->SetStats(0);
		LAPPDTrans->Draw("COLZ");
		_LAPPD_sim_MC_time->Modified();
		_LAPPD_sim_MC_time->Update();

		_LAPPD_sim_MC_time->cd(2);
//		LAPPDPara->SetMarkerStyle(21);
//		LAPPDPara->SetMarkerSize(1.5);
		LAPPDPara->SetStats(0);
		LAPPDPara->Draw("COLZ");
		_LAPPD_sim_MC_time->Modified();
		_LAPPD_sim_MC_time->Update();
		}

		LAPPDMCHits->Clear();
		LAPPDTrans->Clear();
		LAPPDPara->Clear();
}


void LAPPDDisplay::FinaliseHistoAllLAPPDs(){
	_all_hits->GetXaxis()->SetTitle("Radius [m]");
	_all_hits->GetYaxis()->SetTitle("Height [m]");
	_all_hits->GetZaxis()->SetTitle("Arrival time [ns]");
	_all_hits->GetZaxis()->SetTitleOffset(1.4);
	_all_hits->Write();
	if(_config_number == 2)
	{
	_LAPPD_sim_MC_all->cd();
	_all_hits->SetStats(0);
	_all_hits->Draw("COLZ");
	_LAPPD_sim_MC_all->Modified();
	_LAPPD_sim_MC_all->Update();
	}
	_all_hits->Clear();
}


void LAPPDDisplay::RecoDrawing(int eventCounter, unsigned long tubeNumber, std::vector<Waveform<double>> waveformVector)
{
	std::string eventnumber = boost::lexical_cast < std::string > (eventCounter);
	std::string lappdnumber = boost::lexical_cast < std::string > (tubeNumber);
	string nameleft = "event" + eventnumber + "lappd" + lappdnumber + "left";
	const char *d = nameleft.c_str();
	double nbinsx = 25.6;
	_left = new TH2D(d, d, 256, 0.0, nbinsx, 30, 0.0, 30.0);
	string nameright = "event" + eventnumber + "lappd" + lappdnumber + "right";
	const char *e = nameright.c_str();
	_right = new TH2D(e, e, 256, 0.0, nbinsx, 30, 0.0, 30.0);

	for (int i = 0; i < 30; i++)
	{
		std::vector<double>* samplesleft = waveformVector[i].GetSamples();
		std::vector<double>* samplesright = waveformVector[59 - i].GetSamples();
		for (int j = 0; j < samplesleft->size(); j++)
		{
			double time = (double)j/10;
			_left->Fill(time + 0.05, i, -samplesleft->at(j));
			_right->Fill(time + 0.05, i, -samplesright->at(j));
		}
	}

	_left->GetXaxis()->SetTitle("Time [ns]");
	_left->GetYaxis()->SetTitle("Strip number");
	_left->GetZaxis()->SetTitle("Voltage [mV]");
	//_left->GetZaxis()->SetTitleOffset(1.4);
	_left->Write();

	_right->GetXaxis()->SetTitle("Time [ns]");
	_right->GetYaxis()->SetTitle("Strip number");
	_right->GetZaxis()->SetTitle("Voltage [mV]");
	//_right->GetZaxis()->SetTitleOffset(1.4);
	_right->Write();

	if(_config_number == 2)
	{
	_LAPPD_sim_MC_waveform->cd(1);
	_left->SetStats(0);
	_left->Draw("COLZ");
	_LAPPD_sim_MC_waveform->Modified();
	_LAPPD_sim_MC_waveform->Update();

	_LAPPD_sim_MC_waveform->cd(2);
	_right->SetStats(0);
	_right->Draw("COLZ");
	_LAPPD_sim_MC_waveform->Modified();
	_LAPPD_sim_MC_waveform->Update();
	}
	_left->Clear();
	_right->Clear();

}

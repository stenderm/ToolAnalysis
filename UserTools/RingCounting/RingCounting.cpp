/**
 * Class RingCounting:
 * Contains functionality to calculate the number of Cherenkov rings/disks from the MC truth information.
 * The possibility to write events to .csv files is also implemented to use the resulting files with ML clustering methods.
 */

#include "RingCounting.h"
#include <iostream>
#include <fstream>

/**
 * Method readAndPrintConfigParameters
 * Reads configuration parameters from the config file and prints them.
 * A verbose parameter and a csv parameter is used. The csv parameter
 */

void RingCounting::readAndPrintConfigParameters() {

	//Get values from config file
	m_variables.Get("verbose", verbose);
	m_variables.Get("ColumnOrRow", _csv_config);

	//Ensure that the verbose value is in the expected range
	if (verbose < 0) {
		std::cout << "Verbose under 0, set to 0" << std::endl;
		verbose = 0;
	}
	if (verbose > 3) {
		std::cout << "Verbose over 3, set to 3" << std::endl;
		verbose = 3;
	}
	std::cout << "Config parameters" << std::endl;
	std::cout << "Verbose level " << verbose << std::endl;
	std::cout << "_csv_config " << _csv_config << std::endl;
}

/**
 * Method readandPrintGeometryParameters:
 * Gets tank radius, height, center point and number of tank PMTs from the Geometry class, prints them and stores them in member variables.
 */

void RingCounting::readandPrintGeometryParameters() {
	m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", geom);
	tank_radius = geom->GetTankRadius();
	tank_height = geom->GetTankHalfheight();
	n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
	detector_center = geom->GetTankCentre();
	geom->Print();
	//Something is here still wrong in the geometry class, so the values are set here manually
	tank_height /= 2;
	//tank_height = 1.98;
	double barrel_compression = 0.82;
	tank_height *= barrel_compression;

	std::cout << "tank radius before adjusting: " << tank_radius << std::endl;
	if (tank_radius < 1.0)
		tank_radius = 1.37504; //set tank radius to the standard value of old anniev2 configuration(v4 seems to have a very different radius?)

	std::cout << "rank radius: " << tank_radius << std::endl;
	std::cout << "tank half height: " << tank_height << std::endl;

	n_lappds = 0;
}

/**
 * Method readMCdata:
 * Reads the MC data from ANNIEEvent and stores them in member variables.
 */

void RingCounting::readMCdata() {
	std::cout << "Start reading MCdata" << std::endl;
	m_data->Stores["ANNIEEvent"]->Get("MCHits", mchits);
	m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
	m_data->Stores["ANNIEEvent"]->Get("RunNumber", runnumber);
	m_data->Stores["ANNIEEvent"]->Get("SubRunNumber", subrunnumber);
	m_data->Stores["ANNIEEvent"]->Get("MCLAPPDHits", mclappdhits);
	m_data->Stores["ANNIEEvent"]->Get("MCParticles", MCParticles);
	std::cout << "Stop reading MCdata" << std::endl;
}

/**
 * Standard constructor
 *
 */
RingCounting::RingCounting() :
		Tool() {
}

/**
 * Method Initialise:
 * Reads the data from the configfile and from ANNIEEvent, stores it in member variables and prints them.
 * Starts a TApplication and opens a TFile.
 * @param  configfile Path of the configfile
 * @param  data       Pointer to the ANNIE-data.
 * @return            Boolean, which is false when the method encounters an error
 */

bool RingCounting::Initialise(std::string configfile, DataModel &data) {

	/////////////////// Usefull header ///////////////////////
	if (configfile != "")
		m_variables.Initialise(configfile); //loading config file

	m_data = &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////

	readAndPrintConfigParameters();

	readandPrintGeometryParameters();

	//Define ROOT specific parameters
	//TApplication
	int myargc = 0;
	char *myargv[] = { (const char*) "somestring" };
	ringCountingApp = new TApplication("ringCountingApp", &myargc, myargv);
	//TCanvas
	Double_t canvwidth = 700;
	Double_t canvheight = 600;
	ringCountingCanv = new TCanvas("ringCountingCanv", "ringCountingCanv",
			canvwidth, canvheight);
	//Output file
	_file = new TFile(
			"/nashome/m/mstender/WorkingToolAnalysis/Rings.root",
			"Recreate");
	//Histogram for displaying the ring count over all events
	_ring_events = new TH1D("rings", "rings", 200, 0, 20);
	//csv-file, in which the values for Lilia's clustering algorithms are written
	_csvfile.open(
			"/nashome/m/mstender/ToolAnalysis-MrdEfficiency2/Example.csv");
	_hit_number = 0;
	return true;
}


/**
 * Method Execute:
 *
 * @return Boolean, which is false when the method encounters an error.
 */

bool RingCounting::Execute() {
	readMCdata();
	Detector* thepmt;
	Position testposition;
	int tubeIdPMTHits = 0;
	std::cout << "Executing RingCounting..." << std::endl;

//Initialize charge and time arrays
	Position zero_position(0., 0., 0.);
	for (int i = 0; i < 200; i++) {
		charge_pmt[i] = 0;
		channelkey_pmt[i] = 0;
		position_pmt[i] = zero_position;
	}
//initialize vector for LAPPD hits as well
	time_lappd.clear();
	std::vector<double> temp_vec;
	time_lappd.assign(n_lappds, temp_vec);

	hits_LAPPDs.clear();
	std::vector<Position> vec_zero;

	vec_zero.push_back(zero_position);
	hits_LAPPDs.assign(n_lappds, vec_zero);

//get hits and associated charges/times
	Detector* pmt;
	total_hits_pmts = 0;
	int channelmin = 100000;
	int channelmax = 0;
	std::map<unsigned long, std::vector<MCHit> >::iterator itr;
	for (itr = mchits.begin(); itr != mchits.end(); itr++) {
		unsigned long ck = itr->first;
		std::vector<MCHit> TheHit = itr->second;
		for (int i_hit = 0; i_hit < TheHit.size(); i_hit++) {
			if (channelmin > TheHit.at(i_hit).GetTubeId()) {
				channelmin = TheHit.at(i_hit).GetTubeId();
			} else if (channelmax < TheHit.at(i_hit).GetTubeId()) {
				channelmax = TheHit.at(i_hit).GetTubeId();
			}
		}
	}

	int pmtHitCounter = 0;
	for (itr = mchits.begin(); itr != mchits.end(); itr++) {
		unsigned long ck = itr->first;
		std::vector<MCHit> TheHit = itr->second;
		for (int i_hit = 0; i_hit < TheHit.size(); i_hit++) {
			if (TheHit.at(i_hit).GetCharge() > 0.1) {
				pmt = geom->ChannelToDetector(TheHit.at(i_hit).GetTubeId());
				channelkey_pmt[TheHit.at(i_hit).GetTubeId() - channelmin] =
						TheHit.at(i_hit).GetTubeId();
				charge_pmt[TheHit.at(i_hit).GetTubeId() - channelmin] +=
						TheHit.at(i_hit).GetCharge();
				position_pmt[TheHit.at(i_hit).GetTubeId() - channelmin] =
						pmt->GetDetectorPosition() - detector_center;
				pmtHitCounter++;
			}
		}
	}
	for (int i = 0; i < 200; i++) {
		if (charge_pmt[i] > 0) {
			pmtHitCounter++;
		}
	}
	std::cout << "PMTs hit " << pmtHitCounter << std::endl;

	std::cout << "read in mc lappd hits" << std::endl;
	total_hits_lappds = 0;
	num_lappds_hit = 0;
	std::cout << "size of mclappdhits: " << mclappdhits->size() << std::endl;
	if (mclappdhits) {
		for (auto&& achannel : *mclappdhits) {

			unsigned long chankey = achannel.first;
			auto& hits = achannel.second;

			for (auto&& ahit : hits) {
				std::vector<double> temp_pos = ahit.GetPosition();
				int tubeid_lappd = ahit.GetTubeId();
				total_hits_lappds++;
				if ((lappds_selected && act_lappds[tubeid_lappd] == 1)
						|| !lappds_selected) {
					Position lappd_hit(temp_pos.at(0), temp_pos.at(1),
							temp_pos.at(2));
					hits_LAPPDs.at(tubeid_lappd - 1).push_back(lappd_hit);
					charge_lappd[tubeid_lappd - 1] += 1;
					maximum_lappds++;
					total_hits_lappds++;
					double single_time_lappd = ahit.GetTime();
					time_lappd[tubeid_lappd - 1].push_back(single_time_lappd);
				}
			}
			num_lappds_hit++;
		}
		//total_hits_lappds++;
	}
	std::cout << "mclappdhits: " << total_hits_lappds << std::endl;

	std::vector<Position> startVertices;
	std::vector<Position> ringcenters;
	std::vector<Position> vectorsToTank;
	std::vector<double> cherenkovAngles;
	std::vector<double> ringradii;
	std::vector<std::pair<double, int> > ringcharacteristics;
	int counterino = 0;
	CherenkovParameters *CherenkovParameterCalculator = new CherenkovParameters(
			tank_radius, tank_height, detector_center, 1.34);
	int ringcounter = 0;
	for (int i_mcparticle = 0; i_mcparticle < MCParticles->size();
			i_mcparticle++) {
		CherenkovParameterCalculator->setParameter(
				MCParticles->at(i_mcparticle).GetPdgCode(),
				MCParticles->at(i_mcparticle).GetStartEnergy(),
				MCParticles->at(i_mcparticle).GetStartVertex(),
				MCParticles->at(i_mcparticle).GetStartDirection(),
				MCParticles->at(i_mcparticle).GetStopVertex());

		if (CherenkovParameterCalculator->createsARing()) {
			MCParticles->at(i_mcparticle).Print();

			CherenkovParameterCalculator->numberOfPossibleHits(position_pmt,
					charge_pmt, n_tank_pmts, mclappdhits);

			if ((CherenkovParameterCalculator->numberOfPmtHit() > 10)
					|| (CherenkovParameterCalculator->numberOfLappdHit() > 20)) {

				ringcharacteristics.push_back(
						std::make_pair(
								CherenkovParameterCalculator->returnCherenkovRadius(),
								counterino));
				ringradii.push_back(
						CherenkovParameterCalculator->returnCherenkovRadius());
				ringcenters.push_back(
						CherenkovParameterCalculator->returnRingCenter());
				cherenkovAngles.push_back(
						CherenkovParameterCalculator->returnCherenkovAngle());
				vectorsToTank.push_back(
						CherenkovParameterCalculator->returnVectorAsPositionToTank());
				startVertices.push_back(
						CherenkovParameterCalculator->returnsVertexWrtTankCenter(
								MCParticles->at(i_mcparticle).GetStartVertex()));
				counterino++;
			}
		}
	}

	std::vector<int> overlappingRingsMulti;
	if (ringradii.size() > 0) {
		for (int i = 0; i < ringradii.size(); i++) {
			double maxRingRadius = 0;
			int containingRingNumber = 0;
			for (int j = 0; j < ringradii.size(); j++) {
				if (ringradii[i] != 0 && ringradii[j] != 0 && i != j) {
					Position ringCentreToVertexij = ringcenters[j]
							- startVertices[i];
					Position ringCentreToVertexji = ringcenters[i]
							- startVertices[j];
					double dotij = vectorsToTank[i] * ringCentreToVertexij;
					double dotji = vectorsToTank[j] * ringCentreToVertexji;
					double angleij = dotij
							/ ((ringCentreToVertexij.Mag())
									* (vectorsToTank[i].Mag()));
					double angleji = dotji
							/ ((ringCentreToVertexji.Mag())
									* (vectorsToTank[j].Mag()));
					if (angleij <= cherenkovAngles[i]
							|| angleji <= cherenkovAngles[j]) {
						if (maxRingRadius <= ringradii[j]) {
							maxRingRadius = ringradii[j];
							containingRingNumber = j;
						}
					}
					if (maxRingRadius <= ringradii[i]) {
						maxRingRadius = ringradii[i];
						containingRingNumber = i;
					}
				}
			}
			overlappingRingsMulti.push_back(containingRingNumber);
		}
	}

	std::vector<int> overlappingRings;

	for (int i = 0; i < overlappingRingsMulti.size(); i++) {
		bool seen = false;
		for (int j = 0; j < overlappingRingsMulti.size(); j++) {
			if (i != j
					&& overlappingRingsMulti[i] == overlappingRingsMulti[j]) {
				seen = true;
				break;
			}
		}
		if (!seen) {
			overlappingRings.push_back(i);
			ringcounter++;
		}
	}

	int PmtCounter = 0;
	int PrintedPMTCounter = 0;
	double tempDistance = 0;
	double chargeTemp[200];
	double radiusPmt = 0;
	double phiPmt = 0;
	std::copy(std::begin(charge_pmt), std::end(charge_pmt),
			std::begin(chargeTemp));
	if (ringcounter > 0) {
		_ring_events->Fill(ringcounter);
		std::cout << "Number of rings " << ringcounter << std::endl;
		if(_csv_config ==1){
		_csvfile << "Event number," << evnum << ",Run number," << runnumber << ",Subrun number," << subrunnumber << ",Ring count," << ringcounter << "\n";
		_csvfile << "Detector type,Radius in m,Angle in rad,Height in m,x in m,y in m,z in m, Ring number \n";
		}
		for (int i_pmt = 0; i_pmt < n_tank_pmts; i_pmt++) {
			radiusPmt = sqrt(
					position_pmt[i_pmt].X() * position_pmt[i_pmt].X()
							+ position_pmt[i_pmt].Z()
									* position_pmt[i_pmt].Z());
			phiPmt = TMath::ATan2(position_pmt[i_pmt].Z(),
					position_pmt[i_pmt].X());
			if (charge_pmt[i_pmt] >= 0.00001) {
				PmtCounter++;
				for (int i = 0; i < overlappingRings.size(); i++) {
					Position pmtToVertex = position_pmt[i_pmt]
							- startVertices[overlappingRings[i]];
					double dot = vectorsToTank[overlappingRings[i]]
							* pmtToVertex;
					double angle =
							dot
									/ ((pmtToVertex.Mag())
											* (vectorsToTank[overlappingRings[i]].Mag()));
					if (angle <= cherenkovAngles[overlappingRings[i]]
							&& chargeTemp[i_pmt] > 0.00001) {
						chargeTemp[i_pmt] = 0;
						// std::cout << "PMT " << radiusPmt << "  " << phiPmt << "  " <<  position_pmt[i_pmt].Y() << "  "  <<  position_pmt[i_pmt].X() << "  " <<  position_pmt[i_pmt].Y() << "  " <<  position_pmt[i_pmt].Z() << "  " << i << std::endl;
						if(_csv_config == 1){_csvfile << "PMT," << radiusPmt << "," << phiPmt << "," <<  position_pmt[i_pmt].Y() << ","  <<  position_pmt[i_pmt].X() << "," <<  position_pmt[i_pmt].Y() << "," <<  position_pmt[i_pmt].Z() << "," << i << "\n";}
						if(_csv_config == 0){_csvfile << radiusPmt << "," << phiPmt << ","
								<< position_pmt[i_pmt].Y() << ","
								<< position_pmt[i_pmt].X() << ","
								<< position_pmt[i_pmt].Y() << ","
								<< position_pmt[i_pmt].Z() << ","; }
						PrintedPMTCounter++;
					}
				}
				//  std::cout << "charge " << charge[i_pmt] << std::endl;
				//  std::cout << "Distance between PMT and Ring Center " << CherenkovParameterCalculator->distanceBetweenPositionAndRingCenter(pmtPos) << std::endl;
			}
			if (chargeTemp[i_pmt] >= 0.00001) {
				//std::cout << "PMT " << radiusPmt << "  " << phiPmt << "  " <<  position_pmt[i_pmt].Y() << "  "  <<  position_pmt[i_pmt].X() << "  " <<  position_pmt[i_pmt].Y() << "  " <<  position_pmt[i_pmt].Z() << "  " << overlappingRings.size() << std::endl;
				if(_csv_config == 1){_csvfile << "PMT," << radiusPmt << "," << phiPmt << "," <<  position_pmt[i_pmt].Y() << ","  <<  position_pmt[i_pmt].X() << "," <<  position_pmt[i_pmt].Y() << "," <<  position_pmt[i_pmt].Z() << "," << overlappingRings.size() << "\n";}
				if(_csv_config == 0){_csvfile << radiusPmt << "," << phiPmt << ","
						<< position_pmt[i_pmt].Y() << ","
						<< position_pmt[i_pmt].X() << ","
						<< position_pmt[i_pmt].Y() << ","
						<< position_pmt[i_pmt].Z() << ",";}
				//Print
				PrintedPMTCounter++;
			}
		}
		if (PmtCounter != PrintedPMTCounter) {
			std::cout << "PMTCounter " << PmtCounter << std::endl;
			std::cout << "PMTPrintCounter " << PrintedPMTCounter << std::endl;
		}
		double chargeLAPPDTemp[200];
		int counterMissingHits = 0;
		int counterhits = 0;
		int LAPPDhitcounter = 0;
		int PrintCounter = 0;
		std::copy(std::begin(charge_pmt), std::end(charge_pmt),
				std::begin(chargeLAPPDTemp));
		if (mclappdhits) {
			for (auto&& achannel : *mclappdhits) {
				unsigned long chankey = achannel.first;
				auto& hits = achannel.second;
				for (auto&& ahit : hits) {
					LAPPDhitcounter++;
					std::vector<double> temp_pos = ahit.GetPosition();
					Position pos(temp_pos.at(0), temp_pos.at(1),
							temp_pos.at(2));
					Position lappdPos =
							CherenkovParameterCalculator->returnsVertexWrtTankCenter(
									pos);
					double radiusLappd = sqrt(
							lappdPos.X() * lappdPos.X()
									+ lappdPos.Z() * lappdPos.Z());
					double phiLappd = TMath::ATan2(lappdPos.Z(), lappdPos.X());
					for (int i = 0; i < overlappingRings.size(); i++) {
						counterhits++;
						Position lappdToVertex = lappdPos
								- startVertices[overlappingRings[i]];
						double dot = vectorsToTank[overlappingRings[i]]
								* lappdToVertex;
						double angle =
								dot
										/ ((lappdToVertex.Mag())
												* (vectorsToTank[overlappingRings[i]].Mag()));
						if (angle <= cherenkovAngles[overlappingRings[i]]) {
							// std::cout << "LAPPD " << radiusLappd << "  " << phiLappd << "  " << lappdPos.Y() << "  " << lappdPos.X() << "  " << lappdPos.Y() << "  " << lappdPos.Z() << "  " << i << std::endl;
							if(_csv_config == 1){_csvfile << "LAPPD," << radiusLappd << "," << phiLappd << "," << lappdPos.Y() << "," << lappdPos.X() << "," << lappdPos.Y() << "," << lappdPos.Z() << "," << i << "\n";}
							if(_csv_config == 0){_csvfile << radiusLappd << "," << phiLappd << ","
									<< lappdPos.Y() << "," << lappdPos.X()
									<< "," << lappdPos.Y() << ","
									<< lappdPos.Z() << ",";}
							PrintCounter++;
							break;
						} else {
							counterMissingHits++;
						}
					}
					if (counterMissingHits == counterhits) {
						//std::cout << "LAPPD " << radiusLappd << "  " << phiLappd << "  " << lappdPos.Y() << "  " << lappdPos.X() << "  " << lappdPos.Y() << "  " << lappdPos.Z() << "  " << overlappingRings.size() << std::endl;
						if(_csv_config == 1){_csvfile << "LAPPD," << radiusLappd << "," << phiLappd << "," << lappdPos.Y() << "," << lappdPos.X() << "," << lappdPos.Y() << "," << lappdPos.Z() << "," << overlappingRings.size() << "\n";}
						if(_csv_config == 0){_csvfile << radiusLappd << "," << phiLappd << ","
								<< lappdPos.Y() << "," << lappdPos.X() << ","
								<< lappdPos.Y() << "," << lappdPos.Z() << ",";}
						PrintCounter++;
					}
					counterhits = 0;
					counterMissingHits = 0;
				}
			}
		}
		if (LAPPDhitcounter != PrintCounter) {
			std::cout << "ACHTUNG!" << std::endl;
			std::cout << "LAPPDhitcounter " << LAPPDhitcounter << std::endl;
			std::cout << "PrintCounter " << PrintCounter << std::endl;
		}
		if(_csv_config == 0){
		int hitcounter = ((PrintCounter + PrintedPMTCounter) - 2) * 6;
		while (hitcounter < 4248) {
			_csvfile << "0,";
			hitcounter++;
		}
		_csvfile << ringcounter << "\n";
		}
		if(_csv_config == 0){_csvfile << "\n";}
		std::cout << PrintCounter + PrintedPMTCounter << std::endl;
		if (_hit_number < (PrintCounter + PrintedPMTCounter)) {
			_hit_number = PrintCounter + PrintedPMTCounter;
		}
	}

	return true;

}

bool RingCounting::Finalise() {
	gSystem->ProcessEvents();
	_ring_events->Write();
	if (_ring_events)
		delete _ring_events;
	if (gROOT->FindObject("ringCountingCanv") != nullptr)
		delete ringCountingCanv;
	delete ringCountingApp;
	_file->Close();
	_csvfile.close();
	std::cout << "Maximum hit number " << _hit_number - 2 << std::endl;
	std::cout << "Maximum column number " << (_hit_number - 2) * 6 << std::endl;
	return true;
}

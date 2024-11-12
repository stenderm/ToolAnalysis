#include "TopoRecoConverter.h"

// ToDo: Cleanup commented code

TopoRecoConverter::TopoRecoConverter() :
		Tool(), verbosity_(0), error_verbosity_(1), posDirFile(nullptr),
		recoOutputfile(nullptr), mcTree(nullptr), write_pos_dir_file_(false), spatial_resolution_lappd_(1.0),
		start_channel_number_lappds_(1000) {
	run_mode_list_.push_back("Normal");
	run_mode_list_.push_back("Hybrid");
	run_mode_list_.push_back("LightSeparation");
}

bool TopoRecoConverter::Initialise(std::string configfile, DataModel &data) {

	/////////////////// Useful header ///////////////////////
	if (configfile != "") {
		m_variables.Initialise(configfile); // loading config file
	}
	//m_variables.Print();
	m_variables.Get("verbosity", verbosity_);
	m_variables.Get("rootFile", output_root_file_name_);
	m_variables.Get("datFile", output_dat_file_name_);
	m_variables.Get("runMode", run_mode_);
	m_variables.Get("writePosDirFile", write_pos_dir_file_);
	m_variables.Get("spatialResolutionLAPPD", spatial_resolution_lappd_);

	if(verbosity_ > 0){
		std::cout << "verbosity " << verbosity_ << std::endl;
		std::cout << "rootFile " << output_root_file_name_ << std::endl;
		std::cout << "datFile " << output_dat_file_name_ << std::endl;
		std::cout << "runMode " << run_mode_ << std::endl;
		std::cout << "writePosDirFile " << write_pos_dir_file_ << std::endl;
		std::cout << "spatialResolutionLAPPD " << spatial_resolution_lappd_ << std::endl;
	}

	//Is the specified run mode available? If no, stop the program.
	if (std::find(run_mode_list_.begin(), run_mode_list_.end(), run_mode_) == run_mode_list_.end()) {
		std::string logMessage = "TopoRecoConverter: Unexpected run mode provided! ";
		logMessage.append(run_mode_);
		logMessage.append(" was provided, candidates are ");
		for (size_t iMode = 0; iMode < run_mode_list_.size(); iMode++) {
			logMessage.append(run_mode_list_.at(iMode));
			if (iMode != run_mode_list_.size() - 1) {
				logMessage.append(", ");
			}
		}

		throw std::invalid_argument(logMessage);
	}

	m_data = &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////

	// Normal TTR mode: Hit id, charge and time with a list of PMT/LAPPD Ids
	// Ideal Normal TTR mode
	// Hybrid TTR mode: Hit id, charge and time and Cherenkov flag (all Cherenkov for now) with a list of PMT/LAPPD IDs with one hit per ID in list file
	// Ideal Hybrid mode:
	// Light Separation mode: Hit, id, charge, time and cherenkov flag with a list of PMT/LAPPD IDs with one hit per ID in list file

	//Open files
	// .dat file for positions and directions
	posDirFile = new ofstream(output_dat_file_name_, std::ios::out);
	if (!posDirFile->is_open()) {
		std::string logMessage = "File with name  ";
		logMessage.append(output_dat_file_name_);
		logMessage.append(" is not open!");

	    throw std::ios_base::failure(logMessage);
	}

	// .root file for the MC truth
	const char *outputFileNameChar = output_root_file_name_.c_str();
	recoOutputfile = new TFile(outputFileNameChar, "RECREATE");

	//Create the MC Truth tree
	mcTree = new TTree("MC", "Monte Carlo Truth");

	// Set branch addresses for the MC truth file
	mcTree->Branch("particle_IDs", &pID);
	mcTree->Branch("track_IDs", &pTrackID);
	mcTree->Branch("creator_processes", &pCreatorProcess);
	mcTree->Branch("particle_masses", &pMass);
	mcTree->Branch("particle_momenta", &pMomentum);
	mcTree->Branch("particle_energies", &pEnergy);
	mcTree->Branch("particle_end_momenta", &pEndMomentum);
	mcTree->Branch("particle_end_energies", &pEndEnergy);
	mcTree->Branch("particle_direction", &pDir);
	mcTree->Branch("particle_vertex", &pVertex);
	mcTree->Branch("particle_endtex", &pEndpoint);
	mcTree->Branch("particle_parent_type", &pParentType);
	mcTree->Branch("particle_start_time", &pStartTime);
	mcTree->Branch("particle_stop_time", &pStopTime);

	mcTree->Branch("hit_tube_IDs", &hitTubeID);
	mcTree->Branch("hit_times", &hitTime);
	mcTree->Branch("hit_charges", &hitCharge);
	if (run_mode_ != "Normal") {
		mcTree->Branch("light_flags", &hitFlag);
	}

	return true;
}

//Check if this function runs over each individual event or not?
bool TopoRecoConverter::Execute() {

	//ToDo: Throwing exceptions instead of returning false?

	//Retrieve digits
	//To my knowledge, the digits are already corrected for the fact that the detector center is not (0,0,0)
	//Also, the digits are already in cm
	int isStoreOk = m_data->Stores["RecoEvent"]->Get("RecoDigit", reco_digits_);
	if (not isStoreOk) {
		Log("TopoRecoConverter: Not able to retrieve reco digits! Exiting.", error_verbosity_, verbosity_);
		return false;
	}

	//Retrieve MC truth
	isStoreOk = m_data->Stores["ANNIEEvent"]->Get("MCParticles", mc_particles_);
	if (not isStoreOk) {
		Log("TopoRecoConverter: Not able to retrieve MCParticles! Exiting.", error_verbosity_, verbosity_);
		return false;
	}

	//Retrieve map for connecting PDG codes to particle mass
	isStoreOk = m_data->CStore.Get("PdgMassMap", pdg_code_to_mass_);
	if (not isStoreOk) {
		Log("TopoRecoConverter: Not able to retrieve PdgMassMap! Exiting.", error_verbosity_, verbosity_);
		return false;
	}

	//Retrieve the geometry information
	isStoreOk = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry", geometry_);
	if (not isStoreOk) {
		Log("TopoRecoConverter: Error retrieving Geometry from ANNIEEvent! Exiting.", error_verbosity_, verbosity_);
		return false;
	}

	m_data->CStore.Get("pmt_tubeid_to_channelkey_data", pmtid_to_channelkey_);
	m_data->CStore.Get("lappd_tubeid_to_detectorkey", lappd_tubeid_to_detectorkey_);

	// This is in meters. Have to convert it to cm for the conversion of the coordinates
	// Output is: (0, -0.144649, 1.681)
	center_ = geometry_->GetTankCentre();
//	std::cout << "Tank radius " << geometry_->GetTankRadius() << std::endl;
//	std::cout << "Tank half height " << geometry_->GetTankHalfheight() << std::endl;
	center_.UnitToCentimeter();

	//Clear hit information stores
	hitTubeID.clear();
	hitCharge.clear();
	hitTime.clear();
	hitFlag.clear();

	// reserve space
	size_t numberOfHits = reco_digits_->size();
	hitTubeID.reserve(numberOfHits);
	hitCharge.reserve(numberOfHits);
	hitTime.reserve(numberOfHits);
	hitFlag.reserve(numberOfHits);

	// First the pixels have to be calculated, before there can be a matching between hit position and pixel ID
	if (write_pos_dir_file_) {
		WritePMTPositionsDirectionsToFile();
		WriteLAPPDPositionsDirectionsToFile();
		//ToDo: This has to depend on the operation mode
		write_pos_dir_file_ = false;
	}

	std::cout << "Number of Digits " << reco_digits_->size() << std::endl;

	// loop over digits
	// for PMTs, it's easy, just write channelkey with charge and time and flag into vectors
	// For LAPPDs a segmentation has to happen for simulating the spatial resolution
	for (auto aDigit : *reco_digits_) {
		if (aDigit.GetDigitType() == aDigit.PMT8inch) {
			if (pmtid_to_channelkey_.count(aDigit.GetDetectorID()) != 0) {
				hitTubeID.push_back(pmtid_to_channelkey_.at(aDigit.GetDetectorID()));
				hitCharge.push_back(aDigit.GetCalCharge());
				hitTime.push_back(aDigit.GetCalTime());
				if (run_mode_ != "Normal") {
					hitFlag.push_back("Cherenkov");
				}
			}
		}
		if (aDigit.GetDigitType() == aDigit.lappd_v0) {
			if (lappd_tubeid_to_detectorkey_.count(aDigit.GetDetectorID()) != 0) {
				// find correct LAPPD in map with channelkey
				std::vector<std::pair<int, Position> > TempChannelPositionVector = lappd_pixel_keys_positions_.at(lappd_tubeid_to_detectorkey_.at(aDigit.GetDetectorID()));
				double minDifference = 1000000;
				int pixelKey = 0;
				//Find closest position
				for(size_t iPixel = 0; iPixel < TempChannelPositionVector.size(); iPixel++){
					Position difference = TempChannelPositionVector.at(iPixel).second-aDigit.GetPosition();
					if(difference.Mag() < minDifference){
						minDifference = difference.Mag();
						pixelKey = TempChannelPositionVector.at(iPixel).first;
					}
				}

				hitTubeID.push_back(pixelKey);
				hitCharge.push_back(aDigit.GetCalCharge());
				hitTime.push_back(aDigit.GetCalTime());
				if (run_mode_ != "Normal") {
					hitFlag.push_back("Cherenkov");
				}
			}
		}
	}



	// loop over detectors of type "tank" for the tank PMTs and the "LAPPD" for the LAPPDs

	//apart.first == "LAPPD" ||

	//////////////////////////////////////////////////////////////
	//////////////////Convert MC Truth////////////////////////////
	//////////////////////////////////////////////////////////////

	// clear MC Truth stores
	pID.clear();
	pTrackID.clear();
	pCreatorProcess.clear();
	pMass.clear();
	pMomentum.clear();
	pEnergy.clear();
	pEndMomentum.clear();
	pEndEnergy.clear();
	pDir.clear();
	pVertex.clear();
	pEndpoint.clear();
	pParentType.clear();
	pStartTime.clear();
	pStopTime.clear();

	//reserve enough memory
	int nTracksToReserve = mc_particles_->size();

	//ToDo: Delete Check comments
	pID.reserve(nTracksToReserve); //Check
	pTrackID.reserve(nTracksToReserve); //Check
	pCreatorProcess.reserve(nTracksToReserve); //Check
	pMass.reserve(nTracksToReserve); //Check
	pMomentum.reserve(nTracksToReserve); //Check
	pEnergy.reserve(nTracksToReserve); //Checl
	pEndMomentum.reserve(nTracksToReserve); //Check
	pEndEnergy.reserve(nTracksToReserve); //Check
	pDir.reserve(nTracksToReserve); //Check
	pVertex.reserve(nTracksToReserve); //Check
	pEndpoint.reserve(nTracksToReserve); //Check
	pParentType.reserve(nTracksToReserve); //Check
	pStartTime.reserve(nTracksToReserve); //Check
	pStopTime.reserve(nTracksToReserve); //Check

	//loop over particles, convert them and push them to the corresponding vectors
	for (MCParticle particle : *mc_particles_) {
//		std::cout << "Parent: " << particle.GetParentPdg() << "Particle: " << particle.GetPdgCode() << "Energy " << particle.GetStartEnergy() * conversion_factor_GeV_to_MeV << std::endl;


		//ToDo: Test these calculations against the LoadWCSim tool
		pID.push_back(particle.GetPdgCode());
		pTrackID.push_back(particle.GetParticleID());
		//ToDo: Am I missing something here or is WCSim not saving the process?
		pCreatorProcess.push_back("undefined");
		pMass.push_back(GetMassFromPDGCode(particle.GetPdgCode()));
		pMomentum.push_back(CalculateMomentum(particle.GetPdgCode(), particle.GetStartEnergy()));
		//TTR uses cm, ns and MeV as units
		pEnergy.push_back(particle.GetStartEnergy());
		pEndMomentum.push_back(CalculateMomentum(particle.GetPdgCode(), particle.GetStopEnergy()));
		pEndEnergy.push_back(particle.GetStopEnergy());
		Direction particleStartDirection = particle.GetStartDirection();
		std::vector<double> particleStartDirectionVector { particleStartDirection.X(), particleStartDirection.Y(), particleStartDirection.Z() };
		pDir.push_back(particleStartDirectionVector);
		pVertex.push_back(ConvertPosition(particle.GetStartVertex()));
		pEndpoint.push_back(ConvertPosition(particle.GetStopVertex()));

		// TTR should focus on muon, as should do the MC truth drawing
		// Therefore, the proton is set to have a ParentPDG of 2212 for now
		if(particle.GetParentPdg() == 0 && particle.GetPdgCode() == 2212){
			pParentType.push_back(2212);
		}
		else{
			pParentType.push_back(particle.GetParentPdg());
		}
		pStartTime.push_back(particle.GetStartTime());
		pStopTime.push_back(particle.GetStopTime());
	}

	mcTree->Fill();

	return true;
}

bool TopoRecoConverter::Finalise() {
	recoOutputfile->Write();
	posDirFile->close();
	recoOutputfile->Close();

	delete posDirFile;
	delete recoOutputfile;
	return true;
}

// Return a particle's mass based on its PDG code. Return a default value for particle IDs not found in table.
double TopoRecoConverter::GetMassFromPDGCode(int pdgCode) {
	if (pdg_code_to_mass_.count(pdgCode) != 0) {
		return pdg_code_to_mass_.at(pdgCode);
	} else {
		return default_mass_;
	}
}

//Calculate the momentum based on mass and energy. Mass is looked up from the pdg code to mass table.
double TopoRecoConverter::CalculateMomentum(int pdgCode, double energy) {
	// Energy is saved as a GeV value; mass in the PDG code table is given in MeV
	//->Convert the energy to MeV
	double scaledEnergy = energy * conversion_factor_GeV_to_MeV;
	double mass = GetMassFromPDGCode(pdgCode);
	if (mass == default_mass_) {
		return default_mass_;
	} else {
		return sqrt(pow(scaledEnergy, 2) - pow(mass, 2));
	}
}

// Convert the positions in ToolAnalysis to positions expected in the TTR
// The positions have to be corrected with the tank centre.
std::vector<double> TopoRecoConverter::ConvertPosition(Position position) {
	position.UnitToCentimeter();
	position = position - center_;
	std::vector<double> convertedPosition { position.X(), position.Y(), position.Z() };
	return convertedPosition;
}

// Write the PMT position and direction vectors to the .dat file
void TopoRecoConverter::WritePMTPositionsDirectionsToFile() {
	std::map<std::string, std::map<unsigned long, Detector*> > *detectorMap = geometry_->GetDetectors();
	for (auto apart : *detectorMap) {
		if (apart.first == "Tank") {
			std::map<unsigned long, Detector*> detectorMap = apart.second;
			for (auto aDetector : detectorMap) {
				Detector *aChannel = aDetector.second;
				Position positionChannel = aChannel->GetPositionInTank();
				Direction directionChannel = aChannel->GetDetectorDirection();
				positionChannel.UnitToCentimeter();

				// TTR expects direction vectors to point outwards

				*posDirFile << aDetector.first << " " << positionChannel.X()
											   << " " << positionChannel.Y()
											   << " " << positionChannel.Z()
											   << " " << -directionChannel.X()
											   << " " << -directionChannel.Y()
											   << " " << -directionChannel.Z()
											   << " " << aChannel->GetDetectorType() << "\n";

			}
		}
	}
}

// ToDo: This uses TVector3, which is a legacy class by now. This might need an update for newer ROOT versions!
// Write the LAPPD position and direction vectors to the .dat file
void TopoRecoConverter::WriteLAPPDPositionsDirectionsToFile() {
	double sizeLAPPD = 20.0; //cm
	double epsilon = 0.00001;
	//Check if the size is divisible by the spatial resolution. Otherwise the code will produce unwanted results.
	if (std::fmod(sizeLAPPD, spatial_resolution_lappd_) > epsilon) {
		throw std::invalid_argument("TopoRecoConverter: LAPPD size is not divisible by LAPPD resolution");
	}
	// Define the running variables for the for-loops later
	int numberOfSteps = static_cast<int>(sizeLAPPD / spatial_resolution_lappd_);
	int startValue = -numberOfSteps / 2;
	int endValue = numberOfSteps / 2;

	// loop over detectors
	std::map<std::string, std::map<unsigned long, Detector*> > *detectorMap = geometry_->GetDetectors();
	for (auto apart : *detectorMap) {
		// only look at LAPPDs
		if (apart.first == "LAPPD") {
			std::map<unsigned long, Detector*> detectorMap = apart.second;
			// loop over individual LAPPDs
			for (auto aDetector : detectorMap) {
				Detector *aChannel = aDetector.second;
				Position positionChannel = aChannel->GetPositionInTank();
				Direction directionChannel = aChannel->GetDetectorDirection();
				Direction directionNormalised = directionChannel.Unit();
				// ToDo: This might need a little more thought later
				// For now, it is assumed that LAPPD always stand upright in the tank without any tilt along the y-axis
				// For segmenting the LAPPD, two vectors are then used, one for the height, which can be simply chose as (0,1,0)
				// and one as cross product of the direction (or normal) vector and the height vector
				// These two vectors should be able to reach every point on the LAPPD from its position in the tank and can be used
				// to define the position of the pixels
				TVector3 normalVectorLAPPD = TVector3(directionNormalised.X(), directionNormalised.Y(), directionNormalised.Z());
				TVector3 yDirection = TVector3(0, 1.0, 0);
				TVector3 xzDirection = normalVectorLAPPD.Cross(yDirection);
				xzDirection = xzDirection.Unit();
				// All direction vectors should now have length of one, which in our coordinate system should be cm

				positionChannel.UnitToCentimeter();

				//FIXME: It appears that there might be a problem with the aChannel->GetPositionInTank(); not getting
				// corrected with the detector's center
				positionChannel = positionChannel - center_;
				TVector3 midOfLAPPD = TVector3(positionChannel.X(), positionChannel.Y(), positionChannel.Z());

				std::vector<std::pair<int, Position> > TempVectorPerLAPPD;
				TempVectorPerLAPPD.clear();
				TempVectorPerLAPPD.reserve(numberOfSteps * numberOfSteps);
				for (int iStepY = startValue; iStepY < endValue; iStepY++) {
					for (int iStepXZ = startValue; iStepXZ < endValue; iStepXZ++) {
						TVector3 pixelPosition = TVector3(0, 0, 0);
						if (numberOfSteps % 2 == 0) {
							pixelPosition = midOfLAPPD
												   + (iStepXZ + 0.5) * spatial_resolution_lappd_ * xzDirection
												   + (iStepY + 0.5) * spatial_resolution_lappd_ * yDirection;
						}
						else{
							pixelPosition = midOfLAPPD
									               + iStepXZ * spatial_resolution_lappd_ * xzDirection
												   + iStepY * spatial_resolution_lappd_ * yDirection;
						}

						TempVectorPerLAPPD.push_back(std::make_pair(start_channel_number_lappds_, Position(pixelPosition.X(), pixelPosition.Y(), pixelPosition.Z())));



						//TTR expects direction vectors to point outwards
						*posDirFile << start_channel_number_lappds_ << " " << pixelPosition.X()
																    << " " << pixelPosition.Y()
																    << " " << pixelPosition.Z()
																    << " " << -directionChannel.X()
																    << " " << -directionChannel.Y()
																    << " " << -directionChannel.Z()
																    << " " << aChannel->GetDetectorType() << "\n";
						start_channel_number_lappds_++;
					} //end for loop xy
				} //end for loop z
				lappd_pixel_keys_positions_.insert({aDetector.first, TempVectorPerLAPPD});
			} //end for loop LAPPDs
		} //end if LAPPD
	} //end for loop over all detectors
}


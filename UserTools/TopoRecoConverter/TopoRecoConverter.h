#ifndef TopoRecoConverter_H
#define TopoRecoConverter_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "Position.h"
#include "Direction.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

/**
 * \class TopoRecoConverter
 *
 * This tool is a new iteration of a converter formally written by Felix Benckwitz and is heavily based on his version.
 * Its purpose is the conversion of ANNIE MC data (later also real data) to the data format of the Topological Track Reconstruction (TTR).
 * These are the adjustments that have to be made to the ANNIE format:
 * 1. Adjust coordinate system so that (0,0,0) is center of the tank.
 * 2. Write PMT and LAPPD positions into a position and direction file.
 * 3. Match hits on LAPPDs to LAPPD pixel IDs.
 * 4. Write IDs, time and charge of hits to a file.
 *
 * Two files are therefore generated:
 * .dat file for positions and directions.
 * .root file for MCtruth and hits
*
* $Author: M.Stender $
* $Date: 29.05.2024 $
* Contact: malte.stender@desy.de
*/
class TopoRecoConverter: public Tool {


 public:

  TopoRecoConverter(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.


 private:
  	///////////Config variables///////////
	// Integer threshold for when errors are printed
	int error_verbosity_;

	// Overall verbosity integer
	int verbosity_;

	// output file name for MCtruth and hit information
	std::string output_root_file_name_;

	// output file name for the .dat pmt positon and direction file
	std::string output_dat_file_name_;

	// running mode for the converter
	std::string run_mode_;

	// expected spatial resolution for the LAPPDs in cm
	double spatial_resolution_lappd_;

	// list of possible keywords
	std::vector<std::string> run_mode_list_;

	// stream for the .dat file
	ofstream* posDirFile;


	//ToDo: Refactor to better naming
	///////////MCtruth variables///////////
	//File for MC truth
	TFile* recoOutputfile;
	//TTree for the MC truth
	TTree* mcTree;
	//particle IDs (particle data group ID)
	std::vector<int> pID;
	//The IDs of the tracks of the particles in the current event
	std::vector<int> pTrackID;
	//The creator processes of the particles
	std::vector<std::string> pCreatorProcess;
	//The masses of the particle
	std::vector<double> pMass;
	//The momenta of the particle
	std::vector<double> pMomentum;
	//The energies of the particles
	std::vector<double> pEnergy;
	//The momenta of the particles at the endtexes of the particles
	std::vector<double> pEndMomentum;
	//The energies of the particles at the endtexes of the particles
	std::vector<double> pEndEnergy;
	//The directions of the particles
	std::vector<std::vector<double> > pDir;
	//The vertices of the particles
	std::vector<std::vector<double> > pVertex;
	//The endtexes of the particles
	std::vector<std::vector<double> > pEndpoint;
	//The parent type of the particles
	std::vector<int> pParentType;
	//The timees the particles are born
	std::vector<double> pStartTime;
	//The timees the particles die :(
	std::vector<double> pStopTime;

	//The IDs of the tubes measuring the hits
	std::vector<int> hitTubeID;
	//The times the hits are measured
	std::vector<double> hitTime;
	//The charges of the hits
	std::vector<double> hitCharge;
	//The light flag or the hits (if it is a Cherenkov or scintillation photon)
	std::vector<string> hitFlag;

	//////////Helper Variables//////////
	// Depending on operation mode, the file should only be written once for all events
	bool write_pos_dir_file_;

	// Conversion factor for calculating energy/momentum/mass in MeV from a value given in GeV
	const int conversion_factor_GeV_to_MeV = 1000;

	//ToDo: Negative mass does not make sense and can be used to exclude non-common particles.
	//Not sure, whether that can break something in the later usages of the data,
	//but I also do not want to set it to 0 or some other random number to avoid calculated other variables based on wrong values.
	const double default_mass_ = -1000;

	// The start number for the first LAPPD channel for the segmentation approach
	// ToDo: This does not necessarily be a member variable
	int start_channel_number_lappds_;

	// Storage for the LAPPD pixel position to do matching between digit and detector position
	// map of channelkeys and vectors of detectorID that was written to position file and the position
	std::map<int, std::vector<std::pair<int, Position> > > lappd_pixel_keys_positions_;


	//////////Stores//////////
	// Map for finding channelkey from PMT ID
	std::map<int,unsigned long> pmtid_to_channelkey_;

	// Map for finding channelkey from LAPPD ID
	std::map<int,unsigned long> lappd_tubeid_to_detectorkey_;

	// Store for the MC particles
	std::vector<MCParticle>* mc_particles_ = nullptr;

	// Store for the Reco digits
	std::vector<RecoDigit>* reco_digits_ = nullptr;

	//ANNIE geomety
	Geometry* geometry_ = nullptr;

	//Center of the detector
	Position center_;

	// Map for connecting pdgcode to mass of particle in MeV
	std::map<int,double> pdg_code_to_mass_;

	// Return a particle's mass based on its PDG code. Return a default value for particle IDs not found in table.
	double GetMassFromPDGCode(int pdgCode);

	// Calculate the momentum based on the energy and the mass from the pdg table
	double CalculateMomentum(int pdgCode, double energy);

	// Convert the positions in ToolAnalysis to positions expected in the TTR
	// Y and Z are getting swapped and the position is scaled to cm from m
	std::vector<double> ConvertPosition(Position position);

	// Write the PMT position and direction vectors to the .dat file
	void WritePMTPositionsDirectionsToFile();

	// Write the PMT position and direction vectors to the .dat file
	void WriteLAPPDPositionsDirectionsToFile();
};


#endif

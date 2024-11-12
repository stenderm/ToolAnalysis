#include "MonitorDataQuality.h"
#include <iostream>

MonitorDataQuality::MonitorDataQuality() :Tool(){
	// Do nothing
}

bool MonitorDataQuality::Initialise(std::string configfile, DataModel &data) {
	// Load Configuration File
	if (configfile != ""){
		m_variables.Initialise(configfile);
	}

	//assigning transient data pointer
	m_data = &data;

	// Get configuration variables
	m_variables.Get("verbose", m_verbosity);
	m_variables.Get("outputFileName", m_output_file_name);
	if (m_verbosity) {
		std::cout << "Following Configuration was loaded:\n";
		m_variables.Print();
	}

	// ToDo: Maybe have a placeholder name and add naming depending on the loaded data
	// Open Output File
	m_output_file = new TFile(m_output_file_name.c_str(), "recreate");
	if (!m_output_file->IsOpen()) {
		std::string logMessage = "File with name " + m_output_file_name + " could not have been opened.";
		throw std::invalid_argument(logMessage);
	}
	if (m_verbosity) {
		std::cout << "Output file with name " << m_output_file_name << " recreated and opened.\n";
	}
	// Create Trees ToDo: Maybe add some more if necessary
	m_data_quality_tree_facc = new TTree("DataQualityMonitoringTreeFACC", "ANNIE Data Quality Monitoring Tree FACC");
	m_data_quality_tree_tank_pmt = new TTree("DataQualityMonitoringTreeTankPMT", "ANNIE Data Quality Monitoring Tree Tank PMT");
	m_data_quality_tree_tank_lappd = new TTree("DataQualityMonitoringTreeTankLAPPD", "ANNIE Data Quality Monitoring Tree Tank LAPPD");
	m_data_quality_tree_mrd = new TTree("DataQualityMonitoringTreeMRD", "ANNIE Data Quality Monitoring Tree MRD");
	m_data_quality_tree_trigger = new TTree("DataQualityMonitoringTreeTrigger", "ANNIE Data Quality Monitoring Tree Trigger");

	// Set branches for tank PMTs
	m_data_quality_tree_tank_pmt->Branch("SumOfChargeInAllClusters", &m_sum_cluster_charge_tank_pmt, "SumOfChargeInAllClusters/D");
	m_data_quality_tree_tank_pmt->Branch("SumOfChargeInAllClustersInPE", &m_sum_cluster_charge_tank_pmt_pe, "SumOfChargeInAllClustersInPE/D");
	m_data_quality_tree_tank_pmt->Branch("SumOfCharge", &m_sum_charge_tank_pmt, "SumOfCharge/D");
	m_data_quality_tree_tank_pmt->Branch("SumOfChargeInPE", &m_sum_charge_tank_pmt_pe, "SumOfChargeInPE/D");
	m_data_quality_tree_tank_pmt->Branch("RatioOfChargeInClustersVsAllCharge", &ratio_charge_cluster_vs_all, "RatioOfChargeInClustersVsAllCharge/D");
	m_data_quality_tree_tank_pmt->Branch("RatioOfChargeInClustersVsAllChargeInPE", &ratio_charge_cluster_vs_all, "RatioOfChargeInClustersVsAllChargeInPE/D");


	return true;
}

bool MonitorDataQuality::Execute() {
	//ToDo: Maybe move this stuff to own function and move the maps to the header file?

    //std::map<unsigned long, std::vector<Waveform<unsigned short>>> raw_waveform_map;
    //bool has_raw = m_data->Stores["ANNIEEvent"]->Get("RawADCData",raw_waveform_map);


	//channel key with the vector of hits
	std::map<unsigned long, std::vector<Hit>>* regularHitMap = nullptr;

	//cluster time with included cluster hits
	std::map<double,std::vector<Hit>>* clusterHitMap = nullptr;
	//The position of the clusters is in another store

	bool isStoreOk = m_data->Stores["ANNIEEvent"]->Get("Hits", regularHitMap);
	if (not isStoreOk) {
		throw std::invalid_argument("Could not load the hits map.");
	}

	isStoreOk = m_data->CStore.Get("ClusterMap",clusterHitMap);
	if (not isStoreOk) {
		throw std::invalid_argument("Could not load the cluster hits map.");
	}

	delete regularHitMap;
	delete clusterHitMap;
	return true;
}

bool MonitorDataQuality::Finalise() {
	m_output_file->cd();
	m_data_quality_tree_facc->Write();
	m_data_quality_tree_tank_pmt->Write();
	m_data_quality_tree_tank_lappd->Write();
	m_data_quality_tree_mrd->Write();
	m_data_quality_tree_trigger->Write();
	m_output_file->Close();
	return true;
}

void MonitorDataQuality::ExtractMetricsTankHits(std::map<unsigned long, std::vector<Hit>>* t_regularHitMap,
												std::map<double,std::vector<Hit>>* t_clusterHitMap){



}

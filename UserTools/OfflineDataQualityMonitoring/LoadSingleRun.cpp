/*
 * LoadSingleRun.cc
 *
 *  Created on: Sep 10, 2024
 *      Author: stenderm
 */

#include "LoadSingleRun.h"
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include <map>
#include <memory>
#include "TankInformation.h"

LoadSingleRun::LoadSingleRun() {
    // default constructor
}

LoadSingleRun::~LoadSingleRun() {
    // do nothing
}

void LoadSingleRun::extractNtupleInformation(
        int t_runNumber, const std::vector<std::string> &t_fileNames,
        const std::string& t_tankTreeName, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
    for (auto aFileName : t_fileNames) {
        std::shared_ptr<TFile> nTupleFile { std::make_shared<TFile>(aFileName.c_str(), "READ") };
        if (!nTupleFile || nTupleFile->IsZombie()) {
            if (t_verbosity) {
                std::cout << "File with name " << aFileName
                        << " could not have been opened. Skip the file and the run.";
            }
            return;
        }
        loadInformationTank(nTupleFile, t_tankTreeName,
                            t_verbosity, t_ntupleInformationOneRun);
        //Define trees load some stuff at first, return NTupleInformation object?
        nTupleFile->Close();
    }
}

void LoadSingleRun::loadInformationTrigger(
        std::shared_ptr<TFile> t_inputFile, std::string t_triggerTreeName, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
    TTree *triggerTree { openTree(t_inputFile, t_triggerTreeName, t_verbosity) };
    if (!triggerTree) {
        return;
    }
}

void LoadSingleRun::loadInformationMRD(
        std::shared_ptr<TFile> t_inputFile, std::string t_MRDTreeName, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
//    TTree *mrdTree { openTree(t_inputFile, t_MRDTreeName, t_verbosity) };
//    if (!mrdTree) {
//        return;
//    }
    ROOT::RDataFrame mrdDataFrame = ROOT::RDataFrame(t_MRDTreeName, t_inputFile.get());
    ROOT::RDF::RResultPtr<std::vector<int> > eventNumbers { mrdDataFrame.Take<int>("eventNumber") };
    ROOT::RDF::RResultPtr<std::vector<int> > clusterHits { mrdDataFrame.Take<double>("clusterHits") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterTimes { mrdDataFrame.Take<double>("clusterTime") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterTimesSigma { mrdDataFrame.Take<double>("clusterTimeSigma") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > hitTimes { mrdDataFrame.Take<double>("MRDhitT") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<int> > > detectorIDs { mrdDataFrame.Take<double>("MRDhitDetID") };
    ROOT::RDF::RResultPtr<std::vector<int> > numberOfClusterTracks { mrdDataFrame.Take<double>("numClusterTracks") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > trackLengths { mrdDataFrame.Take<double>("MRDTrackLength") };
}

void LoadSingleRun::loadInformationTank(
        std::shared_ptr<TFile> t_inputFile, std::string t_tankTreeName, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
    ROOT::RDataFrame tankDataframe = ROOT::RDataFrame(t_tankTreeName, t_inputFile.get());

    //ToDo: Save all of these information to the NTupleInformationObject and calculate properties then in the CalculateMetrics class
    ROOT::RDF::RResultPtr<std::vector<int> > eventNumbers { tankDataframe.Take<int>("eventNumber") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterCharges { tankDataframe.Take<double>("clusterCharge") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterTimes { tankDataframe.Take<double>("clusterTime") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterChargesPE { tankDataframe.Take<double>("clusterPE") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterChargesMaxPE { tankDataframe.Take<double>("clusterMaxPE") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterChargeBalance { tankDataframe.Take<double>("clusterChargeBalance") };



    std::map<int, TankInformation> runTankInformation;
    int globalNumberOfClusters { 0 };
    // Assign the properties to the corresponding event number
    for (size_t iEvent = 0; iEvent < eventNumbers->size(); ++iEvent) {
        int eventNumber = (*eventNumbers)[iEvent];
        runTankInformation[eventNumber].clusterCharge.push_back((*clusterCharges)[iEvent]);
        runTankInformation[eventNumber].clusterTime.push_back((*clusterTimes)[iEvent]);
        runTankInformation[eventNumber].clusterChargePE.push_back((*clusterChargesPE)[iEvent]);
        runTankInformation[eventNumber].clusterChargeMaxPE.push_back((*clusterChargesMaxPE)[iEvent]);
        runTankInformation[eventNumber].clusterChargeBalance.push_back((*clusterChargeBalance)[iEvent]);
        runTankInformation[eventNumber].numberOfClusters++;
        globalNumberOfClusters++;
    }
    t_ntupleInformationOneRun->setGlobalClusterNumber(globalNumberOfClusters);
    t_ntupleInformationOneRun->setTankInformation(runTankInformation);

}

void LoadSingleRun::addOnePropertyToVector(int t_entry, int t_cluster, double t_propertyToAdd, std::vector<std::vector<double> > & t_vectorToAddTo){
    if(!t_cluster){
        std::vector<double> tempVector { };
        if(t_entry){
            t_vectorToAddTo.push_back(tempVector);
            tempVector.clear();
        }
    }
}

//ToDo: Find run with LAPPD information to implement this!
void LoadSingleRun::loadInformationLAPPD(
        std::shared_ptr<TFile> t_inputFile, std::string t_LAPPDTreeName, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
    TTree *lappdTree { openTree(t_inputFile, t_LAPPDTreeName, t_verbosity) };
    if (!lappdTree) {
        return;
    }
}

TTree* LoadSingleRun::openTree(std::shared_ptr<TFile> t_inputFile, std::string treeName,
                               int t_verbosity) {
    TTree *theTree { dynamic_cast<TTree*>(t_inputFile->Get(treeName.c_str())) };
    if (!theTree) {
        if (t_verbosity) {
            std::cout << "Could not open TTree with name " << treeName << " in file "
                    << t_inputFile->GetName() << ".\n";
        }
        return nullptr;
    }
    return theTree;
}

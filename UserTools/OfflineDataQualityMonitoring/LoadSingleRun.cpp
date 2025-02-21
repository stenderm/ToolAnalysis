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
#include "MRDInformation.h"

LoadSingleRun::LoadSingleRun() {
    // default constructor
}

LoadSingleRun::~LoadSingleRun() {
    // do nothing
}

void LoadSingleRun::extractNtupleInformation(
        int t_runNumber, const std::vector<std::string> &t_fileNames,
        const std::string &t_tankTreeName, const std::string &t_MRDTreeName, int t_verbosity,
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

        if (aFileName.find("BeamCluster") != std::string::npos) {

            loadInformationANNIEEvent(nTupleFile, t_verbosity, t_ntupleInformationOneRun);
        } else {

            loadInformationTankPhaseIITreeMaker(nTupleFile, t_tankTreeName, t_verbosity,
                    t_ntupleInformationOneRun);
            loadInformationMRDPhaseIITreeMaker(nTupleFile, t_MRDTreeName, t_verbosity,
                    t_ntupleInformationOneRun);
        }
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

void LoadSingleRun::loadInformationANNIEEvent(
        std::shared_ptr<TFile> t_inputFile, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
    ROOT::RDataFrame eventDataFrame = ROOT::RDataFrame("Event", t_inputFile.get());
    /// General stuff
    ROOT::RDF::RResultPtr<std::vector<int> > eventNumbers { eventDataFrame.Take<int>("eventNumber") };

    /// Tank
    ROOT::RDF::RResultPtr<std::vector<int> > numberOfClusters { eventDataFrame.Take<int>(
            "numberOfClusters") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > clusterCharges { eventDataFrame.Take<
            std::vector<double> >("clusterCharge") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > clusterTimes { eventDataFrame.Take<
            std::vector<double> >("clusterTime") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > clusterChargesPE {
            eventDataFrame.Take<std::vector<double> >("clusterPE") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > clusterChargesMaxPE {
            eventDataFrame.Take<std::vector<double> >("clusterMaxPE") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > clusterChargeBalance {
            eventDataFrame.Take<std::vector<double> >("clusterChargeBalance") };

    /// MRD
    ROOT::RDF::RResultPtr<std::vector<int> > numberOfClustersMRD { eventDataFrame.Take<int>(
            "MRDClusterNumber") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > clusterTimesMRD { eventDataFrame.Take<
            std::vector<double> >("MRDClusterTime") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > hitTimesMRD { eventDataFrame.Take<
            std::vector<double> >("MRDhitT") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<int> > > detectorIDsMRD { eventDataFrame.Take<
            std::vector<int> >("MRDhitDetID") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<int> > > clusterIDMRD { eventDataFrame.Take<
            std::vector<int> >("MRDHitClusterIndex") };

    std::map<int, TankInformation> runTankInformation;
    std::map<int, MRDInformation> runMRDInformation;
    int globalNumberOfClusters { 0 };

    for (size_t iEvent = 0; iEvent < eventNumbers->size(); ++iEvent) {
        int eventNumber = (*eventNumbers)[iEvent];
        int numberOfClustersForTriggerTank = (*numberOfClusters)[iEvent];
        runTankInformation[iEvent].numberOfClusters = numberOfClustersForTriggerTank;

        /// Tank
        for (int iCluster = 0; iCluster < numberOfClustersForTriggerTank; iCluster++) {
            runTankInformation[iEvent].clusterCharge.push_back((*clusterCharges)[iEvent][iCluster]);
            runTankInformation[iEvent].clusterTime.push_back((*clusterTimes)[iEvent][iCluster]);
            runTankInformation[iEvent].clusterChargePE.push_back(
                    (*clusterChargesPE)[iEvent][iCluster]);
            runTankInformation[iEvent].clusterChargeMaxPE.push_back(
                    (*clusterChargesMaxPE)[iEvent][iCluster]);
            runTankInformation[iEvent].clusterChargeBalance.push_back(
                    (*clusterChargeBalance)[iEvent][iCluster]);
            globalNumberOfClusters++;
        }

        int numberOfClustersForTriggerMRD = (*numberOfClustersMRD)[iEvent];
        runMRDInformation[iEvent].numberOfClusters = numberOfClustersForTriggerMRD;

        std::vector<double> tempHitTimes { };
        std::vector<int> tempIDs { };
        int clusterID { 0 };
        for (size_t iHit = 0; iHit < (*clusterIDMRD)[iEvent].size(); iHit++) {
            if (clusterID == (*clusterIDMRD)[iEvent][iHit]) {
                tempHitTimes.push_back((*hitTimesMRD)[iEvent][iHit]);
                tempIDs.push_back((*detectorIDsMRD)[iEvent][iHit]);
            } else {
                clusterID = (*clusterIDMRD)[iEvent][iHit];
                runMRDInformation[iEvent].hitTimes.push_back(tempHitTimes);
                runMRDInformation[iEvent].detectorIDs.push_back(tempIDs);
                tempHitTimes.clear();
                tempIDs.clear();
                tempHitTimes.push_back((*hitTimesMRD)[iEvent][iHit]);
                tempIDs.push_back((*detectorIDsMRD)[iEvent][iHit]);
            }
        }
        runMRDInformation[iEvent].hitTimes.push_back(tempHitTimes);
        runMRDInformation[iEvent].detectorIDs.push_back(tempIDs);

        for (int iCluster = 0; iCluster < numberOfClustersForTriggerMRD; iCluster++) {
            runMRDInformation[iEvent].clusterTimes.push_back((*clusterTimesMRD)[iEvent][iCluster]);
        }

    }

    t_ntupleInformationOneRun->setGlobalClusterNumber(globalNumberOfClusters);
    t_ntupleInformationOneRun->setTankInformation(runTankInformation);
    t_ntupleInformationOneRun->setMRDInformation(runMRDInformation);
    t_ntupleInformationOneRun->setTreeMakerVersion(NTupleInformation::ANNIEEvent);
}

void LoadSingleRun::loadInformationMRDPhaseIITreeMaker(
        std::shared_ptr<TFile> t_inputFile, std::string t_MRDTreeName, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
//    TTree *mrdTree { openTree(t_inputFile, t_MRDTreeName, t_verbosity) };
//    if (!mrdTree) {
//        return;
//    }
    ROOT::RDataFrame mrdDataFrame = ROOT::RDataFrame(t_MRDTreeName, t_inputFile.get());
    ROOT::RDF::RResultPtr<std::vector<int> > eventNumbers { mrdDataFrame.Take<int>("eventNumber") };
//    ROOT::RDF::RResultPtr<std::vector<int> > clusterHits { mrdDataFrame.Take<double>("clusterHits") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterTimes { mrdDataFrame.Take<double>(
            "clusterTime") };
//    ROOT::RDF::RResultPtr<std::vector<double> > clusterTimesSigma { mrdDataFrame.Take<double>("clusterTimeSigma") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > hitTimes { mrdDataFrame.Take<
            std::vector<double> >("MRDhitT") };
    ROOT::RDF::RResultPtr<std::vector<std::vector<int> > > detectorIDs { mrdDataFrame.Take<
            std::vector<int> >("MRDhitDetID") };
//    ROOT::RDF::RResultPtr<std::vector<int> > numberOfClusterTracks { mrdDataFrame.Take<double>("numClusterTracks") };
//    ROOT::RDF::RResultPtr<std::vector<std::vector<double> > > trackLengths { mrdDataFrame.Take<double>("MRDTrackLength") };

    std::map<int, MRDInformation> runMRDInformation;

    for (size_t iCluster = 0; iCluster < eventNumbers->size(); ++iCluster) {
        int eventNumber = (*eventNumbers)[iCluster];
        runMRDInformation[eventNumber].clusterTimes.push_back((*clusterTimes)[iCluster]);
        runMRDInformation[eventNumber].detectorIDs.push_back((*detectorIDs)[iCluster]);
        runMRDInformation[eventNumber].hitTimes.push_back((*hitTimes)[iCluster]);
        runMRDInformation[eventNumber].numberOfClusters++;
    }

    t_ntupleInformationOneRun->setMRDInformation(runMRDInformation);
    t_ntupleInformationOneRun->setTreeMakerVersion(NTupleInformation::PhaseII);
}

void LoadSingleRun::loadInformationTankPhaseIITreeMaker(
        std::shared_ptr<TFile> t_inputFile, std::string t_tankTreeName, int t_verbosity,
        std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun) {
    ROOT::RDataFrame tankDataframe = ROOT::RDataFrame(t_tankTreeName, t_inputFile.get());

    //ToDo: Save all of these information to the NTupleInformationObject and calculate properties then in the CalculateMetrics class
    ROOT::RDF::RResultPtr<std::vector<int> > eventNumbers { tankDataframe.Take<int>("eventNumber") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterCharges { tankDataframe.Take<double>(
            "clusterCharge") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterTimes { tankDataframe.Take<double>(
            "clusterTime") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterChargesPE { tankDataframe.Take<double>(
            "clusterPE") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterChargesMaxPE { tankDataframe.Take<double>(
            "clusterMaxPE") };
    ROOT::RDF::RResultPtr<std::vector<double> > clusterChargeBalance { tankDataframe.Take<double>(
            "clusterChargeBalance") };

    std::map<int, TankInformation> runTankInformation;
    int globalNumberOfClusters { 0 };
    // Assign the properties to the corresponding event number
    for (size_t iCluster = 0; iCluster < eventNumbers->size(); ++iCluster) {
        int eventNumber = (*eventNumbers)[iCluster];
        runTankInformation[eventNumber].clusterCharge.push_back((*clusterCharges)[iCluster]);
        runTankInformation[eventNumber].clusterTime.push_back((*clusterTimes)[iCluster]);
        runTankInformation[eventNumber].clusterChargePE.push_back((*clusterChargesPE)[iCluster]);
        runTankInformation[eventNumber].clusterChargeMaxPE.push_back(
                (*clusterChargesMaxPE)[iCluster]);
        runTankInformation[eventNumber].clusterChargeBalance.push_back(
                (*clusterChargeBalance)[iCluster]);
        runTankInformation[eventNumber].numberOfClusters++;
        globalNumberOfClusters++;
    }
    t_ntupleInformationOneRun->setGlobalClusterNumber(globalNumberOfClusters);
    t_ntupleInformationOneRun->setTankInformation(runTankInformation);

}

void LoadSingleRun::addOnePropertyToVector(int t_entry, int t_cluster, double t_propertyToAdd,
                                           std::vector<std::vector<double> > &t_vectorToAddTo) {
    if (!t_cluster) {
        std::vector<double> tempVector { };
        if (t_entry) {
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

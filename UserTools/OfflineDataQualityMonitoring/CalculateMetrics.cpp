/*
 * CalculateMetrics.cc
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#include "CalculateMetrics.h"
#include <iostream>
#include "TMath.h"

CalculateMetrics::CalculateMetrics() {
    // do nothing
}

CalculateMetrics::~CalculateMetrics() {
    // do nothing
}
//ToDo: Calculate Mean Charge per clusters, mean charge per event and mean charge per clusters per event
void CalculateMetrics::calculateTankCharge(
        const std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun,
        std::unique_ptr<RunMetrics> &t_runMetrics) {

    std::map<int, TankInformation> tankInformationOneRun { t_ntupleInformationOneRun->getTankInformation() };

    int numberToReserve = t_ntupleInformationOneRun->getGlobalClusterNumber();

    //int runNumber = t_ntupleInformationOneRun->getRunNumber();
    std::vector<int> numberOfClusterPerEvent { };
    std::vector<double> chargePerCluster { };
    std::vector<double> chargePerEventInClusters { };
    std::vector<double> meanChargePerCluster { };

    std::vector<double> chargePerClusterPE { };
    std::vector<double> chargePerEventInClustersPE { };
    std::vector<double> meanChargePerClusterPE { };

    std::vector<double> maxPEPerClusters { };
    std::vector<double> chargeBalancePerCluster { };

    std::vector<double> timePerClusters { };
    std::vector<double> meanClusterTime { };


    numberOfClusterPerEvent.reserve(numberToReserve);
    chargePerCluster.reserve(numberToReserve);
    chargePerEventInClusters.reserve(numberToReserve);
    meanChargePerCluster.reserve(numberToReserve);

    chargePerClusterPE.reserve(numberToReserve);
    chargePerEventInClustersPE.reserve(numberToReserve);
    meanChargePerClusterPE.reserve(numberToReserve);

    maxPEPerClusters.reserve(numberToReserve);
    chargeBalancePerCluster.reserve(numberToReserve);

    timePerClusters.reserve(numberToReserve);
    meanClusterTime.reserve(numberToReserve);

    int numberOfClustersWithChargeBalanceOverOne { 0 };

    for (const auto& [event, tankInformationOneEvent] : tankInformationOneRun) {
        double chargePerEvent { 0.0 };
        double chargerPerEventPE { 0.0 };
        double timePerEvent { 0.0 };


        for(int iCluster = 0; iCluster < tankInformationOneEvent.numberOfClusters; iCluster++){
            chargePerCluster.push_back(tankInformationOneEvent.clusterCharge.at(iCluster));
            chargePerClusterPE.push_back(tankInformationOneEvent.clusterChargePE.at(iCluster));
            maxPEPerClusters.push_back(tankInformationOneEvent.clusterChargeMaxPE.at(iCluster));
            if(tankInformationOneEvent.clusterChargeBalance.at(iCluster) > 2.0){
                numberOfClustersWithChargeBalanceOverOne++;
            }
            if(!std::isinf(tankInformationOneEvent.clusterChargeBalance.at(iCluster))){
                chargeBalancePerCluster.push_back(tankInformationOneEvent.clusterChargeBalance.at(iCluster));
            }
            timePerClusters.push_back(tankInformationOneEvent.clusterTime.at(iCluster));
            chargePerEvent += tankInformationOneEvent.clusterCharge.at(iCluster);
            chargerPerEventPE += tankInformationOneEvent.clusterChargePE.at(iCluster);
            timePerEvent += tankInformationOneEvent.clusterTime.at(iCluster);

        }

        if(tankInformationOneEvent.numberOfClusters){
            meanChargePerCluster.push_back(chargePerEvent/tankInformationOneEvent.numberOfClusters);
            meanChargePerClusterPE.push_back(chargerPerEventPE/tankInformationOneEvent.numberOfClusters);
            meanClusterTime.push_back(timePerEvent/tankInformationOneEvent.numberOfClusters);
        }
        chargePerEventInClusters.push_back(chargePerEvent);
        chargePerEventInClustersPE.push_back(chargerPerEventPE);
        numberOfClusterPerEvent.push_back(tankInformationOneEvent.numberOfClusters);
    }

    t_runMetrics->setTankChargeValues(
            TMath::Mean(numberOfClusterPerEvent.begin(), numberOfClusterPerEvent.end()),
            TMath::StdDev(numberOfClusterPerEvent.begin(), numberOfClusterPerEvent.end()),
            TMath::Mean(chargePerCluster.begin(), chargePerCluster.end()),
            TMath::StdDev(chargePerCluster.begin(), chargePerCluster.end()),
            TMath::Mean(chargePerEventInClusters.begin(), chargePerEventInClusters.end()),
            TMath::StdDev(chargePerEventInClusters.begin(), chargePerEventInClusters.end()),
            TMath::Mean(meanChargePerCluster.begin(), meanChargePerCluster.end()),
            TMath::StdDev(meanChargePerCluster.begin(), meanChargePerCluster.end()),
            chargePerEventInClusters.size());

    t_runMetrics->setTankChargePEValues(
            TMath::Mean(chargePerClusterPE.begin(), chargePerClusterPE.end()),
            TMath::StdDev(chargePerClusterPE.begin(), chargePerClusterPE.end()),
            TMath::Mean(chargePerEventInClustersPE.begin(), chargePerEventInClustersPE.end()),
            TMath::StdDev(chargePerEventInClustersPE.begin(), chargePerEventInClustersPE.end()),
            TMath::Mean(meanChargePerClusterPE.begin(), meanChargePerClusterPE.end()),
            TMath::StdDev(meanChargePerClusterPE.begin(), meanChargePerClusterPE.end()));

    t_runMetrics->setTankChargeMaxPE(
            TMath::Mean(maxPEPerClusters.begin(), maxPEPerClusters.end()),
            TMath::StdDev(maxPEPerClusters.begin(), maxPEPerClusters.end()));

    t_runMetrics->setTankChargeBalance(
            TMath::Mean(chargeBalancePerCluster.begin(), chargeBalancePerCluster.end()),
            TMath::StdDev(chargeBalancePerCluster.begin(), chargeBalancePerCluster.end()));

    t_runMetrics->setTankTimeValues(
            TMath::Mean(timePerClusters.begin(), timePerClusters.end()),
            TMath::StdDev(timePerClusters.begin(), timePerClusters.end()),
            TMath::Mean(meanClusterTime.begin(), meanClusterTime.end()),
            TMath::StdDev(meanClusterTime.begin(), meanClusterTime.end()));

    t_runMetrics->setNumberOfClusterChargeBalanceAnomaly(numberOfClustersWithChargeBalanceOverOne);


}

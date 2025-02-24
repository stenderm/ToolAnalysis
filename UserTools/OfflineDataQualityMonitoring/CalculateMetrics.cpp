/*
 * CalculateMetrics.cc
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#include "CalculateMetrics.h"
#include "TankInformation.h"
#include <iostream>
#include "TMath.h"
#include "MRDInformation.h"

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

    std::map<int, TankInformation> tankInformationOneRun {
            t_ntupleInformationOneRun->getTankInformation() };

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

    int numberOfClusterWithMaxPEInf { 0 };
    int numberOfClusterWithPEInf { 0 };
    int numberOfClusterChargeBalanceInf { 0 };
    int numberOfClusterChargeInf { 0 };

    int numberOfClusterWithMaxPENan { 0 };
    int numberOfClusterWithPENan { 0 };
    int numberOfClusterChargeBalanceNan { 0 };
    int numberOfClusterChargeNan { 0 };

    int numberOfClusters { 0 };

    for (const auto& [event, tankInformationOneEvent] : tankInformationOneRun) {
        double chargePerEvent { 0.0 };
        double chargerPerEventPE { 0.0 };
        double timePerEvent { 0.0 };
        numberOfClusters += tankInformationOneEvent.numberOfClusters;

        for (int iCluster = 0; iCluster < tankInformationOneEvent.numberOfClusters; iCluster++) {

            if (std::isinf(tankInformationOneEvent.clusterChargePE.at(iCluster))) {
                numberOfClusterWithPEInf++;
            } else if (std::isnan(tankInformationOneEvent.clusterChargePE.at(iCluster))) {
                numberOfClusterWithPENan++;
            } else {
                chargePerClusterPE.push_back(tankInformationOneEvent.clusterChargePE.at(iCluster));
                chargerPerEventPE += tankInformationOneEvent.clusterChargePE.at(iCluster);
            }

            if (std::isinf(tankInformationOneEvent.clusterChargeMaxPE.at(iCluster))) {
                numberOfClusterWithMaxPEInf++;
            } else if (std::isnan(tankInformationOneEvent.clusterChargeMaxPE.at(iCluster))) {
                numberOfClusterWithMaxPENan++;
            } else {
                maxPEPerClusters.push_back(tankInformationOneEvent.clusterChargeMaxPE.at(iCluster));
            }

            if (tankInformationOneEvent.clusterChargeBalance.at(iCluster) > 2.0) {
                numberOfClustersWithChargeBalanceOverOne++;
            }

            if (std::isinf(tankInformationOneEvent.clusterChargeBalance.at(iCluster))) {
                numberOfClusterChargeBalanceInf++;
            } else if (std::isnan(tankInformationOneEvent.clusterChargeBalance.at(iCluster))) {
                numberOfClusterChargeBalanceNan++;
            } else {
                chargeBalancePerCluster.push_back(
                        tankInformationOneEvent.clusterChargeBalance.at(iCluster));
            }

            if (std::isinf(tankInformationOneEvent.clusterCharge.at(iCluster))) {
                numberOfClusterChargeInf++;
            } else if (std::isnan(tankInformationOneEvent.clusterCharge.at(iCluster))) {
                numberOfClusterChargeNan++;
            } else {
                chargePerCluster.push_back(tankInformationOneEvent.clusterCharge.at(iCluster));
                chargePerEvent += tankInformationOneEvent.clusterCharge.at(iCluster);
            }

            timePerEvent += tankInformationOneEvent.clusterTime.at(iCluster);
            timePerClusters.push_back(tankInformationOneEvent.clusterTime.at(iCluster));

        }

        if (tankInformationOneEvent.numberOfClusters) {
            meanChargePerCluster.push_back(
                    chargePerEvent / tankInformationOneEvent.numberOfClusters);
            meanChargePerClusterPE.push_back(
                    chargerPerEventPE / tankInformationOneEvent.numberOfClusters);
            meanClusterTime.push_back(timePerEvent / tankInformationOneEvent.numberOfClusters);
        }
        chargePerEventInClusters.push_back(chargePerEvent);
//        if(chargerPerEventPE < 0){
//            std::cout << chargerPerEventPE << "\n";
//        }
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

    t_runMetrics->setTankChargeMaxPE(TMath::Mean(maxPEPerClusters.begin(), maxPEPerClusters.end()),
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

    t_runMetrics->setInfChargeNumbers(numberOfClusterWithMaxPEInf, numberOfClusterWithPEInf,
                                      numberOfClusterChargeBalanceInf, numberOfClusterChargeInf);

    t_runMetrics->setInfChargeRatios(
            static_cast<double>(numberOfClusterWithMaxPEInf) / static_cast<double>(numberOfClusters),
            static_cast<double>(numberOfClusterWithPEInf) / static_cast<double>(numberOfClusters),
            static_cast<double>(numberOfClusterChargeBalanceInf) / static_cast<double>(numberOfClusters),
            static_cast<double>(numberOfClusterChargeInf) / static_cast<double>(numberOfClusters));

    t_runMetrics->setNanChargeNumbers(numberOfClusterWithMaxPENan, numberOfClusterWithPENan,
                                      numberOfClusterChargeBalanceNan, numberOfClusterChargeNan);

    t_runMetrics->setNanChargeRatios(
            static_cast<double>(numberOfClusterWithMaxPENan) / static_cast<double>(numberOfClusters),
            static_cast<double>(numberOfClusterWithPENan) / static_cast<double>(numberOfClusters),
            static_cast<double>(numberOfClusterChargeBalanceNan) / static_cast<double>(numberOfClusters),
            static_cast<double>(numberOfClusterChargeNan) / static_cast<double>(numberOfClusters));

}

void CalculateMetrics::calculateMRDMetrics(
        const std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun,
        std::unique_ptr<RunMetrics> &t_runMetrics) {
    std::map<int, MRDInformation> mrdInformationOneRun {
            t_ntupleInformationOneRun->getMRDInformation() };
    double late = 3800.0;
    double signalBegin = 1000.0;
    double signalEnd = 2500.0;
    int hitsInWindow { 0 };
    int hitsLate { 0 };
    int hitsAll { 0 };
    int clustersInWindow { 0 };
    int clustersLate { 0 };
    int clustersAll { 0 };

    std::map<int, double> numberOfHitsPerChannelPerEvent;
    int numberOfEvents = mrdInformationOneRun.size();

    for (const auto& [event, mrdInformationOneEvent] : mrdInformationOneRun) {
        clustersAll += mrdInformationOneEvent.numberOfClusters;

        for (int iCluster = 0; iCluster < mrdInformationOneEvent.numberOfClusters; iCluster++) {

            if (mrdInformationOneEvent.clusterTimes.at(iCluster) > late) {
                clustersLate++;
            } else if (mrdInformationOneEvent.clusterTimes.at(iCluster) > signalBegin
                    and mrdInformationOneEvent.clusterTimes.at(iCluster) < signalEnd) {
                clustersInWindow++;
            }

            hitsAll += static_cast<int>(mrdInformationOneEvent.hitTimes.at(iCluster).size());
            for (size_t iHit = 0; iHit < mrdInformationOneEvent.hitTimes.at(iCluster).size();
                    iHit++) {
                int detectorID = mrdInformationOneEvent.detectorIDs.at(iCluster).at(iHit);
                numberOfHitsPerChannelPerEvent[detectorID] += (1.0
                        / static_cast<double>(numberOfEvents));
                if (mrdInformationOneEvent.hitTimes.at(iCluster).at(iHit) > late) {
                    hitsLate++;
                } else if (mrdInformationOneEvent.hitTimes.at(iCluster).at(iHit) > signalBegin
                        and mrdInformationOneEvent.hitTimes.at(iCluster).at(iHit) < signalEnd) {
                    hitsInWindow++;
                }
            } //end for loop hits

        } //end for loop clusters

    } // end for loop events

    double hitsInWindowVsAll = static_cast<double>(hitsInWindow) / static_cast<double>(hitsAll);
    double lateHitsVsAll = static_cast<double>(hitsLate) / static_cast<double>(hitsAll);
    double clustersInWindowVsAll = static_cast<double>(clustersInWindow)
            / static_cast<double>(clustersAll);
    double lateClustersVsAll = static_cast<double>(clustersLate) / static_cast<double>(clustersAll);
    double clusterVsEvents = static_cast<double>(clustersAll) / static_cast<double>(numberOfEvents);
    double hitsPerEvent = static_cast<double>(hitsAll) / static_cast<double>(numberOfEvents);
    t_runMetrics->setMRDMetrics(hitsInWindowVsAll, lateHitsVsAll, clusterVsEvents, hitsPerEvent,
            numberOfHitsPerChannelPerEvent, clustersInWindowVsAll, lateClustersVsAll);

}

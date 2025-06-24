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

    int numberOfEvents = t_ntupleInformationOneRun->getGlobalNumberOfEvents();
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

    std::map<int, std::vector<double> > chargePerTube {};
    std::map<int, std::vector<double> > chargePEPerTube {};
    std::map<int, std::vector<double> > hitsPerTube {};


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

    int numberOfClustersWithChargeBalanceOverTwo { 0 };
    int numberOfClustersWithChargeBalanceOverTen { 0 };

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

        for(const auto& [tubeId, hitInformation] : tankInformationOneEvent.hitsPerTube){
            if(!t_runMetrics->isKeyPresentInMap(chargePerTube, tubeId)){
                chargePerTube[tubeId].reserve(numberOfEvents);
                chargePEPerTube[tubeId].reserve(numberOfEvents);
                hitsPerTube[tubeId].reserve(numberOfEvents);
            }
            double chargePerEvent { 0.0 };
            double chargePEPerEvent { 0.0 };
            int hitPerEvent { 0 };
            for(int iHit = 0; iHit < hitInformation.size(); iHit++){
                chargePerEvent =+ std::get<1>(hitInformation.at(iHit));
                chargePEPerEvent =+ std::get<2>(hitInformation.at(iHit));
                hitPerEvent =+ 1;
            }
            chargePerTube[tubeId].push_back(chargePerEvent);
            chargePEPerTube[tubeId].push_back(chargePEPerEvent);
            hitsPerTube[tubeId].push_back(hitPerEvent);
        }



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
                numberOfClustersWithChargeBalanceOverTwo++;
            }
            if(tankInformationOneEvent.clusterChargeBalance.at(iCluster) > 10.0){
                numberOfClustersWithChargeBalanceOverTen++;
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

    // add zeros to the vector of tube charges and tube hits to account for all of the events, where they did not see a hit
    for(const auto& [tubeId, hitCharge] : chargePerTube){
        int numberOfEventsWithZeroHits = numberOfEvents - chargePerTube[tubeId].size();
        chargePerTube[tubeId].insert(chargePerTube[tubeId].end(), numberOfEventsWithZeroHits, 0.0);
        chargePEPerTube[tubeId].insert(chargePEPerTube[tubeId].end(), numberOfEventsWithZeroHits, 0.0);
        hitsPerTube[tubeId].insert(hitsPerTube[tubeId].end(), numberOfEventsWithZeroHits, 0);
        t_runMetrics->setTankTubeValues(tubeId,
                TMath::Mean(chargePerTube[tubeId].begin(), chargePerTube[tubeId].end()),
                TMath::StdDev(chargePerTube[tubeId].begin(), chargePerTube[tubeId].end()),
                TMath::Mean(chargePEPerTube[tubeId].begin(), chargePEPerTube[tubeId].end()),
                TMath::StdDev(chargePEPerTube[tubeId].begin(), chargePEPerTube[tubeId].end()),
                TMath::Mean(hitsPerTube[tubeId].begin(), hitsPerTube[tubeId].end()),
                TMath::StdDev(hitsPerTube[tubeId].begin(), hitsPerTube[tubeId].end()));
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

    t_runMetrics->setNumberOfClusterChargeBalanceAnomaly(numberOfClustersWithChargeBalanceOverTwo, numberOfClustersWithChargeBalanceOverTen,
            static_cast<double>(numberOfClustersWithChargeBalanceOverTwo) / static_cast<double>(numberOfClusters),
            static_cast<double>(numberOfClustersWithChargeBalanceOverTen) / static_cast<double>(numberOfClusters));

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
    // ToDo: Move this to config variables?
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
    int numberOfMRDTracks = t_ntupleInformationOneRun->getGlobalNumberOfMRDTracks();
    double numberOfMRDTracksPerEvent { 0.0 };
    if(numberOfEvents){
        numberOfMRDTracksPerEvent = static_cast<double>(numberOfMRDTracks) / static_cast<double>(numberOfEvents);
    }
    double numberOfMRDTracksPerCluster { 0.0 };
    if(clustersAll){
        numberOfMRDTracksPerCluster = static_cast<double>(numberOfMRDTracks) / static_cast<double>(clustersAll);
    }

    double hitsInWindowVsAll = static_cast<double>(hitsInWindow) / static_cast<double>(hitsAll);
    double lateHitsVsAll = static_cast<double>(hitsLate) / static_cast<double>(hitsAll);
    double clustersInWindowVsAll = static_cast<double>(clustersInWindow)
            / static_cast<double>(clustersAll);
    double lateClustersVsAll = static_cast<double>(clustersLate) / static_cast<double>(clustersAll);
    double clusterVsEvents = static_cast<double>(clustersAll) / static_cast<double>(numberOfEvents);
    double hitsPerEvent = static_cast<double>(hitsAll) / static_cast<double>(numberOfEvents);
    t_runMetrics->setMRDMetrics(hitsInWindowVsAll, lateHitsVsAll, clusterVsEvents, hitsPerEvent,
            numberOfHitsPerChannelPerEvent, clustersInWindowVsAll, lateClustersVsAll, numberOfMRDTracksPerEvent, numberOfMRDTracksPerCluster);


}

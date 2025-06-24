/*
 * RunMetrics.cc
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#include "RunMetrics.h"
#include <iostream>

RunMetrics::RunMetrics(int t_runNumber, RunMode t_runMode) :
        m_run_number(t_runNumber), m_run_mode(t_runMode) {
    // do nothing

}

RunMetrics::~RunMetrics() {
    // do nothing
}

void RunMetrics::setTankChargeValues(double t_meanClusterNumber, double t_stdClusterNumber,
                                     double t_meanChargePerCluster, double t_stdChargePerCluster,
                                     double t_meanChargePerEventInClusters,
                                     double t_stdChargePerEventInClusters,
                                     double t_meanOfMeanChargePerEvent,
                                     double t_stdOfMeanChargePerEvent, size_t t_numberOfEvents) {
    m_mean_cluster_number = t_meanClusterNumber;
    m_std_cluster_number = t_stdClusterNumber;
    m_mean_charge_per_cluster = t_meanChargePerCluster;
    m_std_charge_per_cluster = t_stdChargePerCluster;
    m_mean_charge_per_event_in_clusters = t_meanChargePerEventInClusters;
    m_std_charge_per_event_in_clusters = t_stdChargePerEventInClusters;
    m_mean_mean_charge_per_event = t_meanOfMeanChargePerEvent;
    m_std_mean_charge_per_event = t_stdOfMeanChargePerEvent;
    m_number_of_events = t_numberOfEvents;
}

void RunMetrics::setTankChargePEValues(double t_meanChargePerClusterPE,
                                       double t_stdChargePerClusterPE,
                                       double t_meanChargePerEventInClustersPE,
                                       double t_stdChargePerEventInClustersPE,
                                       double t_meanOfMeanChargePerEventPE,
                                       double t_stdOfMeanChargePerEventPE) {

    m_mean_charge_per_cluster_PE = t_meanChargePerClusterPE;
    m_std_charge_per_cluster_PE = t_stdChargePerClusterPE;
    m_mean_charge_per_event_in_clusters_PE = t_meanChargePerEventInClustersPE;
    m_std_charge_per_event_in_clusters_PE = t_stdChargePerEventInClustersPE;

    m_mean_mean_charge_per_event_PE = t_meanOfMeanChargePerEventPE;
    m_std_mean_charge_per_event_PE = t_stdOfMeanChargePerEventPE;

}
void RunMetrics::setTankChargeMaxPE(double t_meanMaxPEPerCluster, double t_stdMaxPEPerCluster) {

    m_mean_max_PE_per_cluster = t_meanMaxPEPerCluster;
    m_std_max_PE_per_cluster = t_stdMaxPEPerCluster;

}
void RunMetrics::setTankChargeBalance(double t_meanChargeBalancePerCluster,
                                      double t_stdChargeBalancePerCluster) {

    m_mean_charge_balance_per_cluster = t_meanChargeBalancePerCluster;
    m_std_charge_balance_per_cluster = t_stdChargeBalancePerCluster;

}
void RunMetrics::setTankTimeValues(double t_meanTimePerCluster, double t_stdTimePerCluster,
                                   double t_meanOfMeanTimePerEvent,
                                   double t_stdOfMeanTimePerEvent) {
    m_mean_time_per_cluster = t_meanTimePerCluster;
    m_std_time_per_cluster = t_stdTimePerCluster;

    m_mean_mean_time_per_event = t_meanOfMeanTimePerEvent;
    m_std_mean_time_per_event = t_stdOfMeanTimePerEvent;
}

void RunMetrics::setMRDMetrics(double t_hitsInWindowVsAllHits, double t_lateHitsVsAllHits,
                               double t_clustersPerEvent, double t_HitsPerEvent,
                               std::map<int, double> t_hitsPerEventPerChannel,
                               double t_clustersInWindowVsAllClusters,
                               double t_lateClustersVsAllClusters, double t_tracksPerEvent, double t_tracksPerCluster) {
    m_hits_in_signal_window_vs_all_hits = t_hitsInWindowVsAllHits;
    m_late_hits_vs_all_hits = t_lateHitsVsAllHits;
    m_clusters_per_event = t_clustersPerEvent;
    m_hits_per_event = t_HitsPerEvent;
    m_hits_per_event_per_channel = t_hitsPerEventPerChannel;
    m_clusters_in_signal_window_vs_all_clusters = t_clustersInWindowVsAllClusters;
    m_late_clusters_vs_all_clusters = t_lateClustersVsAllClusters;
    m_tracks_per_event = t_tracksPerEvent;
    m_tracks_per_cluster = t_tracksPerCluster;
}

void RunMetrics::printTankCharge() {
    std::cout << "\n";
    std::cout << "Run Number: " << m_run_number << "\n";
    std::cout << "MeanClusterNumber " << m_mean_cluster_number << " +- " << m_std_cluster_number
            << "\n";
    std::cout << "MeanChargePerCluster " << m_mean_charge_per_cluster << " +- "
            << m_std_charge_per_cluster << "\n";
    std::cout << "MeanChargePerEvent " << m_mean_charge_per_event_in_clusters << " +- "
            << m_std_charge_per_event_in_clusters << "\n";
    std::cout << "MeanClusterBalance " << m_mean_charge_balance_per_cluster << "+- "
            << m_std_charge_balance_per_cluster << "\n";
    std::cout << "Number of Clusters with Charge Balance Anomaly "
            << m_number_of_charge_balances_above_two << "\n";
    std::cout << "Ratio of Charge Balances Above Two "
           << m_ratio_of_charge_balances_above_two << "\n";
    std::cout << "Ratio of Charge Balances Above Ten "
           << m_ratio_of_charge_balances_above_ten << "\n";

}

void RunMetrics::setInfChargeNumbers(int t_numberOfInfMaxPE, int t_numberOfInfPE,
                                     int t_numberOfInfChargeBalance, int t_numberOfInfCharge) {
    m_number_of_inf_max_pe = t_numberOfInfMaxPE;
    m_number_of_inf_pe = t_numberOfInfPE;
    m_number_of_inf_charge_balance = t_numberOfInfChargeBalance;
    m_number_of_inf_charge = t_numberOfInfCharge;
}

void RunMetrics::setInfChargeRatios(double t_ratioOfInfMaxPE, double t_ratioOfInfPE,
                                    double t_ratioOfInfChargeBalance, double t_ratioOfInfCharge) {
    m_ratio_of_inf_max_pe = t_ratioOfInfMaxPE;
    m_ratio_of_inf_pe = t_ratioOfInfPE;
    m_ratio_of_inf_charge_balance = t_ratioOfInfChargeBalance;
    m_ratio_of_inf_charge = t_ratioOfInfCharge;
}

void RunMetrics::setNanChargeNumbers(int t_numberOfNanMaxPE, int t_numberOfNanPE,
                                     int t_numberOfNanChargeBalance, int t_numberOfNanCharge) {
    m_number_of_nan_max_pe = t_numberOfNanMaxPE;
    m_number_of_nan_pe = t_numberOfNanPE;
    m_number_of_nan_charge_balance = t_numberOfNanChargeBalance;
    m_number_of_nan_charge = t_numberOfNanCharge;
}

void RunMetrics::setNanChargeRatios(double t_ratioOfNanMaxPE, double t_ratioOfNanPE,
                                    double t_ratioOfNanChargeBalance, double t_ratioOfNanCharge) {
    m_ratio_of_nan_max_pe = t_ratioOfNanMaxPE;
    m_ratio_of_nan_pe = t_ratioOfNanPE;
    m_ratio_of_nan_charge_balance = t_ratioOfNanChargeBalance;
    m_ratio_of_nan_charge = t_ratioOfNanCharge;
}

void RunMetrics::setTankTubeValues(int t_tubeID, double t_meanNumberOfHitsPerEvent, double t_stdNumberOfHitsPerEvent,
                                   double t_meanChargePerEvent, double t_stdChargePerEvent,
                                   double t_meanChargePEPerEvent, double t_stdChargePEPerEvent){
    m_mean_hits_per_event_per_channel_tank[t_tubeID] = t_meanNumberOfHitsPerEvent;
    m_std_hits_per_event_per_channel_tank[t_tubeID] = t_stdNumberOfHitsPerEvent;
    m_mean_charge_per_event_per_channel_tank[t_tubeID] = t_meanChargePerEvent;
    m_std_charge_per_event_per_channel_tank[t_tubeID] = t_stdChargePerEvent;
    m_mean_charge_pe_per_event_per_channel_tank[t_tubeID] = t_meanChargePEPerEvent;
    m_std_charge_pe_per_event_per_channel_tank[t_tubeID] = t_stdChargePEPerEvent;
}




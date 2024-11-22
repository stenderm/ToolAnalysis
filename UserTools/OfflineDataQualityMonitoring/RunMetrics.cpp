/*
 * RunMetrics.cc
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#include "RunMetrics.h"
#include <iostream>

RunMetrics::RunMetrics(int t_runNumber, RunMode t_runMode): m_run_number(t_runNumber), m_run_mode(t_runMode) {
    // do nothing

}

RunMetrics::~RunMetrics() {
    // do nothing
}

void RunMetrics::setTankChargeValues(double t_meanClusterNumber, double t_stdClusterNumber,
                         double t_meanChargePerCluster, double t_stdChargePerCluster,
                         double t_meanChargePerEventInClusters,
                         double t_stdChargePerEventInClusters, double t_meanOfMeanChargePerEvent, double t_stdOfMeanChargePerEvent, size_t t_numberOfEvents){
  m_mean_cluster_number = t_meanClusterNumber;
  m_std_cluster_number = t_stdClusterNumber;
  m_mean_charge_per_cluster =  t_meanChargePerCluster;
  m_std_charge_per_cluster = t_stdChargePerCluster;
  m_mean_charge_per_event_in_clusters = t_meanChargePerEventInClusters;
  m_std_charge_per_event_in_clusters = t_stdChargePerEventInClusters;
  m_mean_mean_charge_per_event = t_meanOfMeanChargePerEvent;
  m_std_mean_charge_per_event = t_stdOfMeanChargePerEvent;
  m_number_of_events = t_numberOfEvents;
}

void RunMetrics::setTankChargePEValues(double t_meanChargePerClusterPE, double t_stdChargePerClusterPE,
                           double t_meanChargePerEventInClustersPE,
                           double t_stdChargePerEventInClustersPE,
                           double t_meanOfMeanChargePerEventPE,
                           double t_stdOfMeanChargePerEventPE){

    m_mean_charge_per_cluster_PE = t_meanChargePerClusterPE;
    m_std_charge_per_cluster_PE = t_stdChargePerClusterPE;

    m_mean_charge_per_event_in_clusters_PE = t_meanChargePerEventInClustersPE;
    m_std_charge_per_event_in_clusters_PE = t_stdChargePerEventInClustersPE;

    m_mean_mean_charge_per_event_PE = t_meanOfMeanChargePerEventPE;
    m_std_mean_charge_per_event_PE = t_stdOfMeanChargePerEventPE;

}
void RunMetrics::setTankChargeMaxPE(double t_meanMaxPEPerCluster, double t_stdMaxPEPerCluster){

    m_mean_max_PE_per_cluster = t_meanMaxPEPerCluster;
    m_std_max_PE_per_cluster = t_stdMaxPEPerCluster;

}
void RunMetrics::setTankChargeBalance(double t_meanChargeBalancePerCluster,
                          double t_stdChargeBalancePerCluster){

    m_mean_charge_balance_per_cluster = t_meanChargeBalancePerCluster;
    m_std_charge_balance_per_cluster = t_stdChargeBalancePerCluster;


}
void RunMetrics::setTankTimeValues(double t_meanTimePerCluster, double t_stdTimePerCluster,
                       double t_meanOfMeanTimePerEvent, double t_stdOfMeanTimePerEvent){
    m_mean_time_per_cluster = t_meanTimePerCluster;
    m_std_time_per_cluster = t_stdTimePerCluster;

    m_mean_mean_time_per_event = t_meanOfMeanTimePerEvent;
    m_std_mean_time_per_event = t_stdOfMeanTimePerEvent;
}

void RunMetrics::setMRDMetrics(double t_hitsInWindowVsAllHits, double t_lateHitsVsAllHits,
                   double t_clustersPerEvent, double t_HitsPerEvent,
                   std::map<int, double> t_hitsPerEventPerChannel,double t_clustersInWindowVsAllClusters,
                   double t_lateClustersVsAllClusters){
    m_hits_in_signal_window_vs_all_hits = t_hitsInWindowVsAllHits;
    m_late_hits_vs_all_hits = t_lateHitsVsAllHits;
    m_clusters_per_event = t_clustersPerEvent;
    m_hits_per_event = t_HitsPerEvent;
    m_hits_per_event_per_channel = t_hitsPerEventPerChannel;
    m_clusters_in_signal_window_vs_all_clusters = t_clustersInWindowVsAllClusters;
    m_late_clusters_vs_all_clusters = t_lateClustersVsAllClusters;
}


void RunMetrics::printTankCharge(){
    std::cout << "\n";
    std::cout << "Run Number: " << m_run_number << "\n";
    std::cout << "MeanClusterNumber " << m_mean_cluster_number << " +- " << m_std_cluster_number << "\n";
    std::cout << "MeanChargePerCluster " << m_mean_charge_per_cluster << " +- " << m_std_charge_per_cluster << "\n";
    std::cout << "MeanChargePerEvent " << m_mean_charge_per_event_in_clusters << " +- " << m_std_charge_per_event_in_clusters << "\n";
    std::cout << "MeanClusterBalance " << m_mean_charge_balance_per_cluster << "+- " << m_std_charge_balance_per_cluster << "\n";
    std::cout << "Number of Clusters with Charge Balance Anomaly " << m_number_of_charge_balance_anomalies << "\n";
}

/*
 * RunMetrics.hh
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_RUNMETRICS_HH_
#define INCLUDE_RUNMETRICS_HH_

#include <vector>
#include <string>
#include <map>

#include "RunMode.h"
/*
 *
 */
class RunMetrics {
    //ToDo: Refactor this to use a map with string keys?
public:
    RunMetrics(int t_runNumber, RunMode t_runMode);
    virtual ~RunMetrics();
    void setTankChargeValues(double t_meanClusterNumber, double t_stdClusterNumber,
                             double t_meanChargePerCluster, double t_stdChargePerCluster,
                             double t_meanChargePerEventInClusters,
                             double t_stdChargePerEventInClusters,
                             double t_meanOfMeanChargePerEvent, double t_stdOfMeanChargePerEvent,
                             size_t t_numberOfEvents);
    void setTankChargePEValues(double t_meanChargePerClusterPE, double t_stdChargePerClusterPE,
                               double t_meanChargePerEventInClustersPE,
                               double t_stdChargePerEventInClustersPE,
                               double t_meanOfMeanChargePerEventPE,
                               double t_stdOfMeanChargePerEventPE);
    void setTankChargeMaxPE(double t_meanMaxPEPerCluster, double t_stdMaxPEPerCluster);
    void setTankChargeBalance(double t_meanChargeBalancePerCluster,
                              double t_stdChargeBalancePerCluster);
    void setTankTimeValues(double t_meanTimePerCluster, double t_stdTimePerCluster,
                           double t_meanOfMeanTimePerEvent, double t_stdOfMeanTimePerEvent);
    void setNumberOfClusterChargeBalanceAnomaly(int t_numberOfClusterChargeBalanceAnomaly) {
        m_number_of_charge_balance_anomalies = t_numberOfClusterChargeBalanceAnomaly;
    }
    void setInfChargeNumbers(int t_numberOfInfMaxPE, int t_numberOfInfPE,
                             int t_numberOfInfChargeBalance, int t_numberOfInfCharge);

    void setInfChargeRatios(double t_ratioOfInfMaxPE, double t_ratioOfInfPE,
                            double t_ratioOfInfChargeBalance, double t_ratioOfInfCharge);

    void setNanChargeNumbers(int t_numberOfNanMaxPE, int t_numberOfNanPE,
                             int t_numberOfNanChargeBalance, int t_numberOfNanCharge);

    void setNanChargeRatios(double t_ratioOfNanMaxPE, double t_ratioOfNanPE,
                            double t_ratioOfNanChargeBalance, double t_ratioOfNanCharge);
    void setMRDMetrics(double t_hitsInWindowVsAllHits, double t_lateHitsVsAllHits,
                       double t_clustersPerEvent, double t_HitsPerEvent,
                       std::map<int, double> t_hitsPerEventPerChannel,
                       double t_clustersInWindowVsAllClusters, double t_lateClustersVsAllClusters);

    void printTankCharge();
    int getRunNumber() const {
        return m_run_number;
    }
    double getMeanClusterNumber() const {
        return m_mean_cluster_number;
    }
    double getStdClusterNumber() const {
        return m_std_cluster_number;
    }
    double getMeanChargePerCluster() const {
        return m_mean_charge_per_cluster;
    }
    double getStdChargePerCluster() const {
        return m_std_charge_per_cluster;
    }
    double getMeanChargePerEvent() const {
        return m_mean_charge_per_event_in_clusters;
    }
    double getStdChargePerEvent() const {
        return m_std_charge_per_event_in_clusters;
    }
    double getMeanMeanChargePerEvent() const {
        return m_mean_mean_charge_per_event;
    }
    double getStdMeanChargePerEvent() const {
        return m_std_mean_charge_per_event;
    }
    double getMeanChargePerClusterPE() const {
        return m_mean_charge_per_cluster_PE;
    }
    double getStdChargePerClusterPE() const {
        return m_std_charge_per_cluster_PE;
    }
    double getMeanChargePerEventPE() const {
        return m_mean_charge_per_event_in_clusters_PE;
    }
    double getStdChargePerEventPE() const {
        return m_std_charge_per_event_in_clusters_PE;
    }
    double getMeanMeanChargePerEventPE() const {
        return m_mean_mean_charge_per_event_PE;
    }
    double getStdMeanChargePerEventPE() const {
        return m_std_mean_charge_per_event_PE;
    }
    double getMeanMaxPEPerCluster() const {
        return m_mean_max_PE_per_cluster;
    }
    double getStdMaxPEPerCluster() const {
        return m_std_max_PE_per_cluster;
    }
    double getMeanChargeBalancePerCluster() const {
        return m_mean_charge_balance_per_cluster;
    }
    double getStdChargeBalancePerCluster() const {
        return m_std_charge_balance_per_cluster;
    }
    double getMeanTimePerCluster() const {
        return m_mean_time_per_cluster;
    }
    double getStdTimePerCluster() const {
        return m_std_time_per_cluster;
    }
    double getMeanMeanTimePerEvent() const {
        return m_mean_mean_time_per_event;
    }
    double getStdMeanTimePerEvent() const {
        return m_std_mean_time_per_event;
    }
    size_t getNumberOfEvents() const {
        return m_number_of_events;
    }
    RunMode getRunMode() const {
        return m_run_mode;
    }

    double getClustersPerEvent() const {
        return m_clusters_per_event;
    }

    double getHitsInSignalWindowVsAllHits() const {
        return m_hits_in_signal_window_vs_all_hits;
    }

    double getHitsPerEvent() const {
        return m_hits_per_event;
    }

    const std::map<int, double>& getHitsPerEventPerChannel() const {
        return m_hits_per_event_per_channel;
    }

    double getLateHitsVsAllHits() const {
        return m_late_hits_vs_all_hits;
    }

    double getLateClustersVsAllClusters() const {
        return m_late_clusters_vs_all_clusters;
    }

    double getClustersInSignalWindowVsAllClusters() const {
        return m_clusters_in_signal_window_vs_all_clusters;
    }

    int getNumberOfChargeBalanceAnomalies() const {
        return m_number_of_charge_balance_anomalies;
    }

    int getNumberOfInfChargeBalance() const {
        return m_number_of_inf_charge_balance;
    }

    int getNumberOfInfMaxPe() const {
        return m_number_of_inf_max_pe;
    }

    int getNumberOfInfPe() const {
        return m_number_of_inf_pe;
    }

    double getRatioOfInfChargeBalance() const {
        return m_ratio_of_inf_charge_balance;
    }

    double getRatioOfInfMaxPe() const {
        return m_ratio_of_inf_max_pe;
    }

    double getRatioOfInfPe() const {
        return m_ratio_of_inf_pe;
    }

    int getNumberOfInfCharge() const {
        return m_number_of_inf_charge;
    }

    int getNumberOfNanCharge() const {
        return m_number_of_nan_charge;
    }

    int getNumberOfNanChargeBalance() const {
        return m_number_of_nan_charge_balance;
    }

    int getNumberOfNanMaxPe() const {
        return m_number_of_nan_max_pe;
    }

    int getNumberOfNanPe() const {
        return m_number_of_nan_pe;
    }

    double getRatioOfInfCharge() const {
        return m_ratio_of_inf_charge;
    }

    double getRatioOfNanCharge() const {
        return m_ratio_of_nan_charge;
    }

    double getRatioOfNanChargeBalance() const {
        return m_ratio_of_nan_charge_balance;
    }

    double getRatioOfNanMaxPe() const {
        return m_ratio_of_nan_max_pe;
    }

    double getRatioOfNanPe() const {
        return m_ratio_of_nan_pe;
    }

private:
    int m_run_number { };
    size_t m_number_of_events { };
    double m_mean_cluster_number { };
    double m_std_cluster_number { };

    double m_mean_charge_per_cluster { };
    double m_std_charge_per_cluster { };

    double m_mean_charge_per_event_in_clusters { };
    double m_std_charge_per_event_in_clusters { };

    double m_mean_mean_charge_per_event { };
    double m_std_mean_charge_per_event { };

    double m_mean_charge_per_cluster_PE { };
    double m_std_charge_per_cluster_PE { };

    double m_mean_charge_per_event_in_clusters_PE { };
    double m_std_charge_per_event_in_clusters_PE { };

    double m_mean_mean_charge_per_event_PE { };
    double m_std_mean_charge_per_event_PE { };

    double m_mean_max_PE_per_cluster { };
    double m_std_max_PE_per_cluster { };

    double m_mean_charge_balance_per_cluster { };
    double m_std_charge_balance_per_cluster { };

    double m_mean_time_per_cluster { };
    double m_std_time_per_cluster { };

    double m_mean_mean_time_per_event { };
    double m_std_mean_time_per_event { };

    int m_number_of_charge_balance_anomalies { };

    int m_number_of_inf_charge_balance { };
    int m_number_of_inf_max_pe { };
    int m_number_of_inf_pe { };
    int m_number_of_inf_charge { };

    double m_ratio_of_inf_charge_balance { };
    double m_ratio_of_inf_max_pe { };
    double m_ratio_of_inf_pe { };
    double m_ratio_of_inf_charge { };

    int m_number_of_nan_charge_balance { };
    int m_number_of_nan_max_pe { };
    int m_number_of_nan_pe { };
    int m_number_of_nan_charge { };

    double m_ratio_of_nan_charge_balance { };
    double m_ratio_of_nan_max_pe { };
    double m_ratio_of_nan_pe { };
    double m_ratio_of_nan_charge { };

    RunMode m_run_mode { };

    ///MRD Values
    // ToDo: These could also all be produced with mean and std.
    double m_hits_in_signal_window_vs_all_hits { };
    double m_late_hits_vs_all_hits { };
    double m_clusters_in_signal_window_vs_all_clusters { };
    double m_late_clusters_vs_all_clusters { };
    double m_clusters_per_event { };
    double m_hits_per_event { };
    std::map<int, double> m_hits_per_event_per_channel { };

};

#endif /* INCLUDE_RUNMETRICS_HH_ */

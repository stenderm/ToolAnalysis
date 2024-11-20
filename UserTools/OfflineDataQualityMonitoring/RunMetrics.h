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
    void setNumberOfClusterChargeBalanceAnomaly(int t_numberOfClusterChargeBalanceAnomaly){
        m_number_of_charge_balance_anomalies = t_numberOfClusterChargeBalanceAnomaly;
    }

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

    RunMode m_run_mode { };

};

#endif /* INCLUDE_RUNMETRICS_HH_ */

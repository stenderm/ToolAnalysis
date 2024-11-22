/*
 * NTupleInformation.hh
 *
 *  Created on: Sep 10, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_NTUPLEINFORMATION_HH_
#define INCLUDE_NTUPLEINFORMATION_HH_

#include <vector>
#include <map>
#include "TankInformation.h"
#include "MRDInformation.h"

/*
 *  This class holds the information of a run based on the loaded ntuples.
 */
class NTupleInformation {
public:
    NTupleInformation(int t_runNumber);
    virtual ~NTupleInformation();
    void setTankInformation(std::map<int, TankInformation> t_tankInformation);
    void setMRDInformation(std::map<int, MRDInformation> t_mrdInformation);
    void setGlobalClusterNumber(int t_globalNumberOfClusters);
    int getGlobalClusterNumber() const{
        return m_global_number_of_clusters;
    }
    std::map<int, TankInformation> getTankInformation() const{
        return m_tank_information;
    }
    std::map<int, MRDInformation> getMRDInformation() const{
        return m_mrd_information;
    }
    int getRunNumber() const{
        return m_run_number;
    }

private:
    // Map that holds all of the tank information with the event number as the key and a struct as value
    std::map<int, TankInformation> m_tank_information { };
    // Map that holds all of the MRD information with the event number as key and the MRDInformation struct as value
    std::map<int, MRDInformation> m_mrd_information { };
    // The run number of the current NTupleInformation
    int m_run_number { };
    int m_global_number_of_clusters { };
};

#endif /* INCLUDE_NTUPLEINFORMATION_HH_ */

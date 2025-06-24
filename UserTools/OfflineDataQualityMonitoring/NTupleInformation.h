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
#include <string>

/*
 *  This class holds the information of a run based on the loaded ntuples.
 */
class NTupleInformation {
public:
    //ToDo: Add more if there will be a new version, but also add the specific use cases to loading of the runs and calculation of the metrics
    enum TreeMakerVersion{
        PhaseII = 0,
        ANNIEEvent = 1
    };
    NTupleInformation(int t_runNumber);
    virtual ~NTupleInformation();
    void setTankInformation(std::map<int, TankInformation> t_tankInformation);
    void setMRDInformation(std::map<int, MRDInformation> t_mrdInformation);
    void setGlobalClusterNumber(int t_globalNumberOfClusters);
    void setGlobalMRDTrackNumber(int t_globalNumberMRDTracks){m_global_number_of_mrd_tracks = t_globalNumberMRDTracks;}
    void setGlobalNumberOfEvents(int t_globalNumberOfEvents){m_global_number_of_events = t_globalNumberOfEvents;}
    int getGlobalClusterNumber() const {
        return m_global_number_of_clusters;
    }
    std::map<int, TankInformation> getTankInformation() const {
        return m_tank_information;
    }
    std::map<int, MRDInformation> getMRDInformation() const {
        return m_mrd_information;
    }
    int getRunNumber() const {
        return m_run_number;
    }

    const TreeMakerVersion& getTreeMakerVersion() const {
        return m_tree_maker_version;
    }

    void setTreeMakerVersion(TreeMakerVersion t_treeMakerVersion) {
        m_tree_maker_version = t_treeMakerVersion;
    }

    int getGlobalNumberOfMRDTracks() const{
        return m_global_number_of_mrd_tracks;
    }

    int getGlobalNumberOfEvents() const{
        return m_global_number_of_events;
    }

private:
    // Map that holds all of the tank information with the event number as the key and a struct as value
    std::map<int, TankInformation> m_tank_information { };
    // Map that holds all of the MRD information with the event number as key and the MRDInformation struct as value
    std::map<int, MRDInformation> m_mrd_information { };
    // The run number of the current NTupleInformation
    int m_run_number { };
    int m_global_number_of_clusters { };
    TreeMakerVersion m_tree_maker_version { };
    int m_global_number_of_mrd_tracks { };
    int m_global_number_of_events { };
};

#endif /* INCLUDE_NTUPLEINFORMATION_HH_ */

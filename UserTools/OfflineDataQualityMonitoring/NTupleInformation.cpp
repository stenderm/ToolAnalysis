/*
 * NTupleInformation.cc
 *
 *  Created on: Sep 10, 2024
 *      Author: stenderm
 */

#include "NTupleInformation.h"

NTupleInformation::NTupleInformation(int t_runNumber) :
        m_run_number(t_runNumber) {

}

NTupleInformation::~NTupleInformation() {
    // do nothing
}

void NTupleInformation::setTankInformation(
        std::map<int, TankInformation > t_tankInformation) {
    if(m_tank_information.empty()){
        m_tank_information = t_tankInformation;
    }
    else{
        m_tank_information.insert(std::begin(t_tankInformation), std::end(t_tankInformation));
    }
}

void NTupleInformation::setGlobalClusterNumber(int t_globalNumberOfClusters){
    m_global_number_of_clusters = t_globalNumberOfClusters;
}

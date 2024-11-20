//struct that holds the information of the tank from the NTuples for better organisation
//This is meant to be loaded into a map with an int as key telling the number of the event
#ifndef INCLUDE_TANKINFORMATION_HH_
#define INCLUDE_TANKINFORMATION_HH_

#include <vector>

struct TankInformation{
    int numberOfClusters;
    std::vector<double> clusterChargePE;
    std::vector<double> clusterChargeMaxPE;
    std::vector<double> clusterCharge;
    std::vector<double> clusterTime;
    std::vector<double> clusterChargeBalance;
};

#endif /* INCLUDE_TANKINFORMATION_HH_ */
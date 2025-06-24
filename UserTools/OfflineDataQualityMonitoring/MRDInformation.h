#ifndef INCLUDE_MRDINFORMATION_HH_
#define INCLUDE_MRDINFORMATION_HH_

#include <vector>

struct MRDInformation{
    int numberOfClusters;
    std::vector<double> clusterTimes;
    // ToDo: Refactor to map maybe?
    std::vector<std::vector<int> > detectorIDs;
    std::vector<std::vector<double> > hitTimes;
};

#endif /* INCLUDE_MRDINFORMATION_HH_ */

/*
 * CalculateMetrics.hh
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_CALCULATEMETRICS_HH_
#define INCLUDE_CALCULATEMETRICS_HH_

#include "NTupleInformation.h"
#include "RunMetrics.h"
#include <memory>

/*
 *
 */
class CalculateMetrics {
public:
    CalculateMetrics();
    virtual ~CalculateMetrics();
    void calculateTankCharge(const std::unique_ptr<NTupleInformation>& t_ntupleInformationOneRun,
                             std::unique_ptr<RunMetrics>& t_runMetrics);
    void calculateMRDMetrics(
            const std::unique_ptr<NTupleInformation> &t_ntupleInformationOneRun,
            std::unique_ptr<RunMetrics> &t_runMetrics);
};

#endif /* INCLUDE_CALCULATEMETRICS_HH_ */

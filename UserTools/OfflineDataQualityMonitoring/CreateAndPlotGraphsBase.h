/*
 * CreateAndPlotGraphsBase.hh
 *
 *  Created on: Nov 7, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_CREATEANDPLOTGRAPHSBASE_HH_
#define INCLUDE_CREATEANDPLOTGRAPHSBASE_HH_

#include "RunMetrics.h"
#include <memory>

/*
 *
 */
class CreateAndPlotGraphsBase {
public:
    CreateAndPlotGraphsBase() = default;
    virtual ~CreateAndPlotGraphsBase() = default;
    virtual void setTankCharge(const std::unique_ptr<RunMetrics>& t_runMetrics) = 0;
    virtual void drawAndSave(bool t_saveHistogramsAsPictures, std::string t_saveDirectory, std::string t_fileName, bool t_saveAsRoot) = 0;

};

#endif /* INCLUDE_CREATEANDPLOTGRAPHSBASE_HH_ */

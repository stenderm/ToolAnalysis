/*
 * CreateAndPlotGraphsMultipleRunModes.hh
 *
 *  Created on: Oct 14, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_CREATEANDPLOTGRAPHSMULTIPLERUNMODES_HH_
#define INCLUDE_CREATEANDPLOTGRAPHSMULTIPLERUNMODES_HH_

#include "CanvasAndMultipleGraphs.h"
#include "CreateAndPlotGraphsBase.h"
#include "RunMetrics.h"
/*
 *
 */
class CreateAndPlotGraphsMultipleRunModes : public CreateAndPlotGraphsBase{
public:
    CreateAndPlotGraphsMultipleRunModes();
    virtual ~CreateAndPlotGraphsMultipleRunModes();
    void setTankCharge(const std::unique_ptr<RunMetrics>& t_runMetrics) override;
    void setMRDMetrics(const std::unique_ptr<RunMetrics> &t_runMetrics) override;
    void drawAndSave(bool t_saveHistogramsAsPictures,std::string t_saveDirectory, std::string t_fileName, bool t_saveAsRoot) override;

private:
    std::vector<CanvasAndMultipleGraphs> m_canvases_and_graphs_tank;
    std::vector<CanvasAndMultipleGraphs> m_canvases_and_graphs_mrd;
    void matchFillValueToGraph(const CanvasAndMultipleGraphs& t_graph, const std::unique_ptr<RunMetrics>& t_runMetrics);
    void matchFillValueToGraphMRD(const CanvasAndMultipleGraphs &t_graph, const std::unique_ptr<RunMetrics> &t_runMetrics);
};

#endif /* INCLUDE_CREATEANDPLOTGRAPHSMULTIPLERUNMODES_HH_ */

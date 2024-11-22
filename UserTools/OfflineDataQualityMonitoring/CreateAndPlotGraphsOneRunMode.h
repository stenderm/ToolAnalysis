/*
 * CreateAndPlotGraphs.hh
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_CREATEANDPLOTGRAPHSONERUNMODE_HH_
#define INCLUDE_CREATEANDPLOTGRAPHSONERUNMODE_HH_

#include <memory>
#include "RunMetrics.h"
#include "CanvasAndGraph.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "CreateAndPlotGraphsBase.h"

/*
 *
 */
class CreateAndPlotGraphsOneRunMode : public CreateAndPlotGraphsBase {
public:
    CreateAndPlotGraphsOneRunMode();
    virtual ~CreateAndPlotGraphsOneRunMode();
    void setMRDMetrics(const std::unique_ptr<RunMetrics> &t_runMetrics) override;
    void setTankCharge(const std::unique_ptr<RunMetrics>& t_runMetrics) override;
    void drawAndSave(bool t_saveHistogramsAsPictures, std::string t_saveDirectory, std::string t_fileName, bool t_saveAsRoot) override;
private:
    std::vector<CanvasAndGraph> m_canvases_and_graphs_tank;
    std::vector<CanvasAndGraph> m_canvases_and_graphs_mrd;
    void matchFillValueToGraph(const CanvasAndGraph& t_graph, const std::unique_ptr<RunMetrics>& t_runMetrics);
    void matchFillValueToGraphMRD(const CanvasAndGraph &t_graph, const std::unique_ptr<RunMetrics> &t_runMetrics);
    bool m_draw_error {false};
    bool m_draw_horizontal_error { false };
};

#endif /* INCLUDE_CREATEANDPLOTGRAPHSONERUNMODE_HH_ */

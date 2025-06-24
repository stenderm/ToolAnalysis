//struct that groups canvases and graphs, together with the titles and the axis titles.
#ifndef INCLUDE_CANVASANDGRAPH_HH_
#define INCLUDE_CANVASANDGRAPH_HH_

#include "TCanvas.h"
#include "TGraphErrors.h"
#include <algorithm>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "TAxis.h"
//ToDo: Add a run modifier, that can be used to change the color and to fill the right information.
struct CanvasAndGraph{
    std::unique_ptr<TCanvas> canvas { };
    std::unique_ptr<TGraphErrors> graph { };
    //Most of the cases this should be the run number, 
    //but there can be in the future compilation of runs.
    //It is also imaginable that there could be graphs with other titles.
    std::string xAxisTitle { };
    std::string yAxisTitle { };
    void removeSpaces(std::string& t_string){
        boost::erase_all(t_string, " ");
    }
    CanvasAndGraph(std::string t_xAxisTitle, std::string t_yAxisTitle, std::string t_subsystemIdentifier, int t_width, int t_height):
                   xAxisTitle(t_xAxisTitle), yAxisTitle(t_yAxisTitle){
        std::string canvasTitle { t_yAxisTitle + " In " + t_subsystemIdentifier };
        graph = std::make_unique<TGraphErrors>();
        graph->SetMarkerStyle(21);
        std::string graphTitle { canvasTitle + ";" + t_xAxisTitle + ";" + t_yAxisTitle};
        graph->SetTitle(graphTitle.c_str());
        std::string yAxisStripped = t_yAxisTitle;
        std::string xAxisStripped = t_xAxisTitle;
        graph->GetYaxis()->SetTitleOffset(1.5);
        removeSpaces(yAxisStripped);
        removeSpaces(xAxisStripped);
        std::string canvasName { yAxisStripped + "Per" + xAxisStripped + "In" + t_subsystemIdentifier};
        canvas = std::make_unique<TCanvas>(canvasName.c_str(), canvasTitle.c_str(), t_width, t_height);
    }
};

#endif /* INCLUDE_TANKINFORMATION_HH_ */

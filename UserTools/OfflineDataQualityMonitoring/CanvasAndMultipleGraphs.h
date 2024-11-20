//struct that groups one canvas with multiple graphs in order to plot different run types at once
#ifndef INCLUDE_CANVASANDMULTIPLEGRAPHs_HH_
#define INCLUDE_CANVASANDMULTIPLEGRAPHs_HH_

#include "TCanvas.h"
#include "TGraphErrors.h"
#include <algorithm>
#include <iostream>
#include "RunMode.h"
#include <boost/algorithm/string.hpp>


struct CanvasAndMultipleGraphs {
    std::unique_ptr<TCanvas> canvas { };
    // ToDo: Add run modes later
    std::map<RunMode, std::shared_ptr<TGraphErrors> > graphsAllModes { };
    //Most of the cases this should be the run number, 
    //but there can be in the future compilation of runs.
    //It is also imaginable that there could be graphs with other titles.
    std::string xAxisTitle { };
    std::string yAxisTitle { };
    std::string genericGraphTitle { };
    void removeSpaces(std::string &t_string) {
        boost::erase_all(t_string, " ");
    }
    CanvasAndMultipleGraphs(std::string t_xAxisTitle, std::string t_yAxisTitle,
                            std::string t_subsystemIdentifier, int t_width, int t_height) :
            xAxisTitle(t_xAxisTitle), yAxisTitle(t_yAxisTitle) {
        std::string canvasTitle { t_yAxisTitle + " In " + t_subsystemIdentifier };
        std::string yAxisStripped = t_yAxisTitle;
        std::string xAxisStripped = t_xAxisTitle;
        removeSpaces(yAxisStripped);
        removeSpaces(xAxisStripped);
        std::string canvasName { yAxisStripped + "Per" + xAxisStripped + "In"
                + t_subsystemIdentifier };
        canvas = std::make_unique<TCanvas>(canvasName.c_str(), canvasTitle.c_str(), t_width,
                t_height);
        for (const auto& oneMode : runModeToString) {
            RunMode currentRunMode = oneMode.first;
            std::shared_ptr<TGraphErrors> graph { std::make_shared<TGraphErrors>() };
            switch (currentRunMode) {
                case RunMode::Beam:
                    graph->SetMarkerColor(kBlack);
                    graph->SetLineColor(kBlack);
                    break;
                case RunMode::Cosmic:
                    graph->SetMarkerColor(kBlue);
                    graph->SetLineColor(kBlue);
                    break;
                case RunMode::LED:
                    graph->SetMarkerColor(kRed);
                    graph->SetLineColor(kRed);
                    break;
                case RunMode::AmBe:
                    graph->SetMarkerColor(kGreen);
                    graph->SetLineColor(kGreen);
                    break;
                case RunMode::Laser:
                    graph->SetMarkerColor(kYellow);
                    graph->SetLineColor(kYellow);
                    break;
                default:
                    throw std::invalid_argument("Specified run mode does not exist in enum!");
            }
            graph->SetMarkerStyle(21);

            std::string graphTitle { canvasTitle + oneMode.second + ";" + t_xAxisTitle + ";"
                    + t_yAxisTitle };
            genericGraphTitle = canvasTitle + ";" + t_xAxisTitle + ";" + t_yAxisTitle;
            graph->SetTitle(graphTitle.c_str());
            // This vector is now exactly ordered like the enum
            graphsAllModes.insert({currentRunMode, std::move(graph)});
        }
    }//end construct

    // Delete the copy constructor
    CanvasAndMultipleGraphs(const CanvasAndMultipleGraphs&) = delete;
    // Delete the copy assignment operator
    CanvasAndMultipleGraphs& operator=(const CanvasAndMultipleGraphs&) = delete;

    // Optional: Implement move constructor and move assignment if needed
    CanvasAndMultipleGraphs(CanvasAndMultipleGraphs&&) noexcept = default;
    CanvasAndMultipleGraphs& operator=(CanvasAndMultipleGraphs&&) noexcept = default;

};

#endif /* INCLUDE_CANVASANDMULTIPLEGRAPHs_HH_ */

/*
 * CreateAndPlotGraphs.cc
 *
 *  Created on: Sep 11, 2024
 *      Author: stenderm
 */

#include "CreateAndPlotGraphsOneRunMode.h"
#include "TFile.h"

CreateAndPlotGraphsOneRunMode::CreateAndPlotGraphsOneRunMode() {
    m_draw_error = true;
    m_draw_horizontal_error = true;
    const int width { 1920 };
    const int height { 1080 };
    const std::string regularXAxis { "Run Number" };
    std::string subSystemIdentifier { "Tank" };

    //ToDo: Automate the title and name generation
    //For name just put the axis titles together and get rid of the spaces

    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Number of Clusters", subSystemIdentifier, width, height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Charge per Cluster", subSystemIdentifier, width, height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Charge per Event", subSystemIdentifier, width, height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Number of Events", subSystemIdentifier, width, height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Mean Charge per Event", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Charge per Cluster in PE", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Charge per Event in PE", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Mean Charge per Event in PE", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Max Charge per Cluster in PE", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Charge Balance per Cluster", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Time per Cluster in ns", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndGraph(regularXAxis, "Mean Time per Event in ns", subSystemIdentifier, width,
                    height));

    subSystemIdentifier = "MRD";
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndGraph(regularXAxis, "Ratio Signal Hits to All Hits", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndGraph(regularXAxis, "Ratio Late Hits to All Hits", subSystemIdentifier, width, height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndGraph(regularXAxis, "Ratio Signal Clusters to All Clusters", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndGraph(regularXAxis, "Ratio Late Clusters to All Clusters", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndGraph(regularXAxis, "Number of Clusters per Event", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndGraph(regularXAxis, "Number of Hits per Event", subSystemIdentifier, width,
                    height));
    int numberOfMaxID = 350;
    for (int iID = 0; iID < numberOfMaxID; iID++) {
        std::string title = "Hits per Event for ID " + std::to_string(iID);
        m_canvases_and_graphs_mrd.push_back(
                CanvasAndGraph(regularXAxis, title, subSystemIdentifier, width, height));
    }
}

CreateAndPlotGraphsOneRunMode::~CreateAndPlotGraphsOneRunMode() {
    // do nothing
}

void CreateAndPlotGraphsOneRunMode::setMRDMetrics(const std::unique_ptr<RunMetrics> &t_runMetrics) {
    for (const CanvasAndGraph &oneGraph : m_canvases_and_graphs_mrd) {
        //ToDo: Function to sort the to fill variable to the corresponding graph container
        matchFillValueToGraphMRD(oneGraph, t_runMetrics);
    }
}

void CreateAndPlotGraphsOneRunMode::matchFillValueToGraphMRD(
        const CanvasAndGraph &t_graph, const std::unique_ptr<RunMetrics> &t_runMetrics) {
    double errorToFill { 0.0 };
    double pointToFill { 0.0 };

    if (t_graph.yAxisTitle == "Ratio Signal Hits to All Hits") {
        pointToFill = t_runMetrics->getHitsInSignalWindowVsAllHits();
        errorToFill = pointToFill/100.0;
    } else if (t_graph.yAxisTitle == "Ratio Late Hits to All Hits") {
        pointToFill = t_runMetrics->getLateHitsVsAllHits();
        errorToFill = pointToFill/100.0;
    } else if (t_graph.yAxisTitle == "Ratio Signal Clusters to All Clusters") {
        pointToFill = t_runMetrics->getClustersInSignalWindowVsAllClusters();
        errorToFill = pointToFill/100.0;
    } else if (t_graph.yAxisTitle == "Ratio Late Clusters to All Clusters") {
        pointToFill = t_runMetrics->getLateClustersVsAllClusters();
        errorToFill = pointToFill/100.0;
    } else if (t_graph.yAxisTitle == "Number of Clusters per Event") {
        pointToFill = t_runMetrics->getClustersPerEvent();
        errorToFill = pointToFill/100.0;
    } else if (t_graph.yAxisTitle == "Number of Hits per Event") {
        pointToFill = t_runMetrics->getHitsPerEvent();
        errorToFill = pointToFill/100.0;
    } else {
        for (const auto &IdValuePair : t_runMetrics->getHitsPerEventPerChannel()) {
            std::string title = "Hits per Event for ID " + std::to_string(IdValuePair.first);
            if (t_graph.yAxisTitle == title) {
                pointToFill = IdValuePair.second;
                errorToFill = pointToFill/100.0;
            }
        }
    }
    if(!m_draw_error){
        errorToFill = 0.0;
    }
    double horizontalError = 0.0;
    if(m_draw_horizontal_error){
        horizontalError = 0.25;
    }

    t_graph.graph->AddPoint(t_runMetrics->getRunNumber(), pointToFill);
    t_graph.graph->SetPointError(t_graph.graph->GetN() - 1, horizontalError, errorToFill);

}

void CreateAndPlotGraphsOneRunMode::setTankCharge(const std::unique_ptr<RunMetrics> &t_runMetrics) {
    for (const CanvasAndGraph &oneGraph : m_canvases_and_graphs_tank) {
        //ToDo: Function to sort the to fill variable to the corresponding graph container
        matchFillValueToGraph(oneGraph, t_runMetrics);
    }
}

//ToDo: Currently everything expects to be in terms of run number; add other cases?
void CreateAndPlotGraphsOneRunMode::matchFillValueToGraph(
        const CanvasAndGraph &t_graph, const std::unique_ptr<RunMetrics> &t_runMetrics) {
    double errorToFill { 0.0 };
    double pointToFill { 0.0 };

    // Lots of if cases, it's sad, but I know no better solution T.T
    if (t_graph.yAxisTitle == "Number of Clusters") {
        pointToFill = t_runMetrics->getMeanClusterNumber();
        errorToFill = t_runMetrics->getStdClusterNumber();
    } else if (t_graph.yAxisTitle == "Charge per Cluster") {
        pointToFill = t_runMetrics->getMeanChargePerCluster();
        errorToFill = t_runMetrics->getStdChargePerCluster();
    } else if (t_graph.yAxisTitle == "Charge per Event") {
        pointToFill = t_runMetrics->getMeanChargePerEvent();
        errorToFill = t_runMetrics->getStdChargePerEvent();
    } else if (t_graph.yAxisTitle == "Number of Events") {
        pointToFill = t_runMetrics->getNumberOfEvents();
        errorToFill = 0.0;
    } else if (t_graph.yAxisTitle == "Mean Charge per Event") {
        pointToFill = t_runMetrics->getMeanMeanChargePerEvent();
        errorToFill = t_runMetrics->getStdMeanChargePerEvent();
    } else if (t_graph.yAxisTitle == "Charge per Cluster in PE") {
        pointToFill = t_runMetrics->getMeanChargePerClusterPE();
        errorToFill = t_runMetrics->getStdChargePerClusterPE();
    } else if (t_graph.yAxisTitle == "Charge per Event in PE") {
        pointToFill = t_runMetrics->getMeanChargePerEventPE();
        errorToFill = t_runMetrics->getStdChargePerEventPE();
    } else if (t_graph.yAxisTitle == "Mean Charge per Event in PE") {
        pointToFill = t_runMetrics->getMeanMeanChargePerEventPE();
        errorToFill = t_runMetrics->getStdMeanChargePerEventPE();
    } else if (t_graph.yAxisTitle == "Max Charge per Cluster in PE") {
        pointToFill = t_runMetrics->getMeanMaxPEPerCluster();
        errorToFill = t_runMetrics->getStdMaxPEPerCluster();
    } else if (t_graph.yAxisTitle == "Charge Balance per Cluster") {
        pointToFill = t_runMetrics->getMeanChargeBalancePerCluster();
        errorToFill = t_runMetrics->getStdChargeBalancePerCluster();
    } else if (t_graph.yAxisTitle == "Time per Cluster in ns") {
        pointToFill = t_runMetrics->getMeanTimePerCluster();
        errorToFill = t_runMetrics->getStdTimePerCluster();
    } else if (t_graph.yAxisTitle == "Mean Time per Event in ns") {
        pointToFill = t_runMetrics->getMeanMeanTimePerEvent();
        errorToFill = t_runMetrics->getStdMeanTimePerEvent();
    }
    double horizontalError = 0.0;
    if(m_draw_horizontal_error){
        horizontalError = 0.25;
    }
    t_graph.graph->AddPoint(t_runMetrics->getRunNumber(), pointToFill);
    t_graph.graph->SetPointError(t_graph.graph->GetN() - 1, horizontalError, errorToFill);

}

void CreateAndPlotGraphsOneRunMode::drawAndSave(bool t_saveHistogramsAsPictures,
                                                std::string t_saveDirectory, std::string t_fileName,
                                                bool t_saveAsRoot) {

    std::unique_ptr<TFile> outputROOTFile;
    if (t_saveAsRoot) {
        outputROOTFile = std::make_unique < TFile > (t_fileName.c_str(), "RECREATE");
        if (!outputROOTFile->IsOpen()) {
            std::cout << "Output file " << t_fileName
                    << " couldn't be opened. Pictures are not saved!\n";
            outputROOTFile->Close();
            return;
        }
    }
    // loop over all histogram containers
    for (size_t iGraph { 0 }; iGraph < m_canvases_and_graphs_tank.size(); iGraph++) {
        // set to the correct canvasS
        m_canvases_and_graphs_tank.at(iGraph).canvas->cd();
        //draw
        m_canvases_and_graphs_tank.at(iGraph).graph->Draw("AP");
        if (t_saveAsRoot and outputROOTFile) {
            outputROOTFile->cd();
            m_canvases_and_graphs_tank.at(iGraph).canvas->Write(
                    m_canvases_and_graphs_tank.at(iGraph).canvas->GetName());
        }
        // save histograms as pictures if specified
        if (t_saveHistogramsAsPictures) {
            std::string saveFileName = t_saveDirectory
                    + std::string(m_canvases_and_graphs_tank.at(iGraph).canvas->GetName()) + ".pdf";
            m_canvases_and_graphs_tank.at(iGraph).canvas->Print(saveFileName.c_str());
        }
    }

    // loop over all histogram containers
    for (size_t iGraph { 0 }; iGraph < m_canvases_and_graphs_mrd.size(); iGraph++) {
        if(!m_canvases_and_graphs_mrd.at(iGraph).graph->GetN()){
            continue;
        }
        // set to the correct canvasS
        m_canvases_and_graphs_mrd.at(iGraph).canvas->cd();
        //draw
        m_canvases_and_graphs_mrd.at(iGraph).graph->Draw("AP");
        if (t_saveAsRoot and outputROOTFile) {
            outputROOTFile->cd();
            m_canvases_and_graphs_mrd.at(iGraph).canvas->Write(
                    m_canvases_and_graphs_mrd.at(iGraph).canvas->GetName());
        }
        // save histograms as pictures if specified
        if (t_saveHistogramsAsPictures) {
            std::string saveFileName = t_saveDirectory
                    + std::string(m_canvases_and_graphs_mrd.at(iGraph).canvas->GetName()) + ".pdf";
            m_canvases_and_graphs_mrd.at(iGraph).canvas->Print(saveFileName.c_str());
        }
    }



    if (t_saveAsRoot and outputROOTFile) {
        outputROOTFile->Close();
    }
}

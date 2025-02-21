/*
 * CreateAndPlotGraphsMultipleRunModes.cc
 *
 *  Created on: Oct 14, 2024
 *      Author: stenderm
 */

#include "CreateAndPlotGraphsMultipleRunModes.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TFile.h"

CreateAndPlotGraphsMultipleRunModes::CreateAndPlotGraphsMultipleRunModes() {
    const int width { 1920 };
    const int height { 1080 };
    const std::string regularXAxis { "Run Number" };
    std::string subSystemIdentifier { "Tank" };

    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Number of Clusters", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Cluster", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Event", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Number of Events", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Mean Charge per Event", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Cluster in PE", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Event in PE", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Mean Charge per Event in PE", subSystemIdentifier,
                    width, height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Max Charge per Cluster in PE", subSystemIdentifier,
                    width, height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge Balance per Cluster", subSystemIdentifier,
                    width, height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Time per Cluster in ns", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_tank.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Mean Time per Event in ns", subSystemIdentifier,
                    width, height));
    subSystemIdentifier = "MRD";
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Ratio Signal Hits to All Hits", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Ratio Late Hits to All Hits", subSystemIdentifier, width, height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Ratio Signal Clusters to All Clusters", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_mrd.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Ratio Late Clusters to All Clusters", subSystemIdentifier, width,
                    height));
    m_canvases_and_graphs_mrd.push_back(CanvasAndMultipleGraphs(regularXAxis, "Number of Clusters per Event", subSystemIdentifier, width, height));
    m_canvases_and_graphs_mrd.push_back(CanvasAndMultipleGraphs(regularXAxis, "Number of Hits per Event", subSystemIdentifier, width, height));
    int numberOfMaxID  = 350;
    for ( int iID = 0; iID < numberOfMaxID; iID++) {
        std::string title = "Hits per Event for ID " + std::to_string(iID);
        m_canvases_and_graphs_mrd.push_back(CanvasAndMultipleGraphs(regularXAxis, title, subSystemIdentifier, width, height));
    }
}

CreateAndPlotGraphsMultipleRunModes::~CreateAndPlotGraphsMultipleRunModes() {
    // do nothing
}

void CreateAndPlotGraphsMultipleRunModes::setTankCharge(
        const std::unique_ptr<RunMetrics> &t_runMetrics) {
    for (const CanvasAndMultipleGraphs &oneGraph : m_canvases_and_graphs_tank) {
        //ToDo: Function to sort the to fill variable to the corresponding graph container
        matchFillValueToGraph(oneGraph, t_runMetrics);
    }
}

void CreateAndPlotGraphsMultipleRunModes::setMRDMetrics(const std::unique_ptr<RunMetrics> &t_runMetrics){
    for (const CanvasAndMultipleGraphs &oneGraph : m_canvases_and_graphs_mrd) {
        //ToDo: Function to sort the to fill variable to the corresponding graph container
        matchFillValueToGraphMRD(oneGraph, t_runMetrics);
    }
}
void CreateAndPlotGraphsMultipleRunModes::matchFillValueToGraphMRD(const CanvasAndMultipleGraphs &t_graph, const std::unique_ptr<RunMetrics> &t_runMetrics){
    double errorToFill { 0.0 };
    double pointToFill { 0.0 };

    if (t_graph.yAxisTitle == "Ratio Signal Hits to All Hits") {
        pointToFill = t_runMetrics->getHitsInSignalWindowVsAllHits();
    } else if (t_graph.yAxisTitle == "Ratio Late Hits to All Hits") {
        pointToFill = t_runMetrics->getLateHitsVsAllHits();
    } else if (t_graph.yAxisTitle == "Ratio Signal Clusters to All Clusters") {
        pointToFill = t_runMetrics->getClustersInSignalWindowVsAllClusters();
    } else if (t_graph.yAxisTitle == "Ratio Late Clusters to All Clusters") {
        pointToFill = t_runMetrics->getLateClustersVsAllClusters();
    } else if (t_graph.yAxisTitle == "Number of Clusters per Event") {
        pointToFill = t_runMetrics->getClustersPerEvent();
    } else if (t_graph.yAxisTitle == "Number of Hits per Event") {
        pointToFill = t_runMetrics->getHitsPerEvent();
    } else {
        for (const auto &IdValuePair : t_runMetrics->getHitsPerEventPerChannel()) {
            std::string title = "Hits per Event for ID " + std::to_string(IdValuePair.first);
            if (t_graph.yAxisTitle == title) {
                pointToFill = IdValuePair.second;
            }
        }
    }

    t_graph.graphsAllModes.at(t_runMetrics->getRunMode())->AddPoint(t_runMetrics->getRunNumber(),
            pointToFill);
    t_graph.graphsAllModes.at(t_runMetrics->getRunMode())->SetPointError(
            t_graph.graphsAllModes.at(t_runMetrics->getRunMode())->GetN() - 1, 0.0, errorToFill);

}

//ToDo: Currently everything expects to be in terms of run number; add other cases?
void CreateAndPlotGraphsMultipleRunModes::matchFillValueToGraph(
        const CanvasAndMultipleGraphs &t_graph, const std::unique_ptr<RunMetrics> &t_runMetrics) {
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

    t_graph.graphsAllModes.at(t_runMetrics->getRunMode())->AddPoint(t_runMetrics->getRunNumber(),
            pointToFill);
    t_graph.graphsAllModes.at(t_runMetrics->getRunMode())->SetPointError(
            t_graph.graphsAllModes.at(t_runMetrics->getRunMode())->GetN() - 1, 0.0, errorToFill);

}

void CreateAndPlotGraphsMultipleRunModes::drawAndSave(bool t_saveHistogramsAsPictures,
                                                      const std::string& t_saveDirectory,
                                                      const std::string& t_fileName, bool t_saveAsRoot) {
    std::unique_ptr<TFile> outputROOTFile;
    if (t_saveAsRoot) {
        outputROOTFile = std::make_unique<TFile>(t_fileName.c_str(), "RECREATE");
        if (!outputROOTFile->IsOpen()) {
            std::cout << "Output file " << t_fileName
                    << " couldn't be opened. Pictures are not saved!\n";
            outputROOTFile->Close();
            return;
        }
    }
    // loop over all histogram containers
    for (size_t iGraph { 0 }; iGraph < m_canvases_and_graphs_tank.size(); iGraph++) {
        // set to the correct canvases
        m_canvases_and_graphs_tank.at(iGraph).canvas->cd();
        std::string graphTitle = m_canvases_and_graphs_tank.at(iGraph).genericGraphTitle;
        std::unique_ptr<TMultiGraph> theMultiGraph = std::make_unique<TMultiGraph>();
        std::unique_ptr<TLegend> theLegend = std::make_unique<TLegend>(0.8, 0.7, 0.9, 0.9);
        for (auto &aGraph : m_canvases_and_graphs_tank.at(iGraph).graphsAllModes) {
            theMultiGraph->Add(aGraph.second.get(), "AP");
            // ToDo: Add RunMode as specifier
            theLegend->AddEntry(aGraph.second.get(), getRunModeName(aGraph.first).c_str(), "p");
        }

        theMultiGraph->SetTitle(graphTitle.c_str());
        theMultiGraph->Draw("AP");

        theLegend->Draw();
        if (t_saveAsRoot and outputROOTFile) {
            outputROOTFile->cd();
            m_canvases_and_graphs_tank.at(iGraph).canvas->Write(
                    m_canvases_and_graphs_tank.at(iGraph).genericGraphTitle.c_str());
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
        // set to the correct canvases
        m_canvases_and_graphs_mrd.at(iGraph).canvas->cd();
        std::string graphTitle = m_canvases_and_graphs_mrd.at(iGraph).genericGraphTitle;
        std::unique_ptr<TMultiGraph> theMultiGraph = std::make_unique<TMultiGraph>();
        std::unique_ptr<TLegend> theLegend = std::make_unique<TLegend>(0.8, 0.7, 0.9, 0.9);
        for (auto &aGraph : m_canvases_and_graphs_mrd.at(iGraph).graphsAllModes) {
            theMultiGraph->Add(aGraph.second.get(), "AP");
            // ToDo: Add RunMode as specifier
            theLegend->AddEntry(aGraph.second.get(), getRunModeName(aGraph.first).c_str(), "p");
        }

        theMultiGraph->SetTitle(graphTitle.c_str());
        theMultiGraph->Draw("AP");

        theLegend->Draw();
        if (t_saveAsRoot and outputROOTFile) {
            outputROOTFile->cd();
            m_canvases_and_graphs_mrd.at(iGraph).canvas->Write(
                    m_canvases_and_graphs_mrd.at(iGraph).genericGraphTitle.c_str());
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


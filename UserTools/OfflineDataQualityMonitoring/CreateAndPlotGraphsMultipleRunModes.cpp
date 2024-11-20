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
    const std::string tankIdentifier { "Tank" };

    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Number of Clusters", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Cluster", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Event", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Number of Events", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Mean Charge per Event", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Cluster in PE", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge per Event in PE", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Mean Charge per Event in PE", tankIdentifier,
                    width, height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Max Charge per Cluster in PE", tankIdentifier,
                    width, height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Charge Balance per Cluster", tankIdentifier,
                    width, height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Time per Cluster in ns", tankIdentifier, width,
                    height));
    m_canvases_and_graphs.push_back(
            CanvasAndMultipleGraphs(regularXAxis, "Mean Time per Event in ns", tankIdentifier,
                    width, height));

}

CreateAndPlotGraphsMultipleRunModes::~CreateAndPlotGraphsMultipleRunModes() {
    // do nothing
}

void CreateAndPlotGraphsMultipleRunModes::setTankCharge(
        const std::unique_ptr<RunMetrics> &t_runMetrics) {
    for (const CanvasAndMultipleGraphs &oneGraph : m_canvases_and_graphs) {
        //ToDo: Function to sort the to fill variable to the corresponding graph container
        matchFillValueToGraph(oneGraph, t_runMetrics);
    }
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
                                                      std::string t_saveDirectory,
                                                      std::string t_fileName, bool t_saveAsRoot) {
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
    for (size_t iGraph { 0 }; iGraph < m_canvases_and_graphs.size(); iGraph++) {
        // set to the correct canvases
        m_canvases_and_graphs.at(iGraph).canvas->cd();
        std::string graphTitle = m_canvases_and_graphs.at(iGraph).genericGraphTitle;
        std::unique_ptr<TMultiGraph> theMultiGraph = std::make_unique<TMultiGraph>();
        std::unique_ptr<TLegend> theLegend = std::make_unique<TLegend>(0.8, 0.7, 0.9, 0.9);
        for (auto &aGraph : m_canvases_and_graphs.at(iGraph).graphsAllModes) {
            theMultiGraph->Add(aGraph.second.get(), "AP");
            // ToDo: Add RunMode as specifier
            theLegend->AddEntry(aGraph.second.get(), getRunModeName(aGraph.first).c_str(), "p");
        }

        theMultiGraph->SetTitle(graphTitle.c_str());
        theMultiGraph->Draw("AP");

        theLegend->Draw();
        if (t_saveAsRoot and outputROOTFile) {
            outputROOTFile->cd();
            m_canvases_and_graphs.at(iGraph).canvas->Write(
                    m_canvases_and_graphs.at(iGraph).genericGraphTitle.c_str());
        }

        // save histograms as pictures if specified
        if (t_saveHistogramsAsPictures) {
            std::string saveFileName = t_saveDirectory
                    + std::string(m_canvases_and_graphs.at(iGraph).canvas->GetName()) + ".pdf";
            m_canvases_and_graphs.at(iGraph).canvas->Print(saveFileName.c_str());
        }
    }
    if (t_saveAsRoot and outputROOTFile) {
        outputROOTFile->Close();
    }
}


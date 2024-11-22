#include "OfflineDataQualityMonitoring.h"
#include "LoadFileList.h"
#include "LoadSingleRun.h"
#include "NTupleInformation.h"
#include "CalculateMetrics.h"
#include "RunMetrics.h"
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <memory>
#include "CreateAndPlotGraphsOneRunMode.h"
#include "CreateAndPlotGraphsMultipleRunModes.h"
#include "CreateAndPlotGraphsBase.h"

OfflineDataQualityMonitoring::OfflineDataQualityMonitoring() :
        Tool() {
}

bool OfflineDataQualityMonitoring::Initialise(std::string configfile, DataModel &data) {

    m_run_modes = { "Beam", "Cosmic", "LED", "AmBe", "Laser", "Any" };

    loadAndCheckConfiguration(configfile);

    return true;
}

bool OfflineDataQualityMonitoring::Execute() {
    std::vector<std::tuple<int, RunMode, std::vector<std::string> > > fileNamesAsPairs =
            loadListFile();
    LoadSingleRun loadSingleRunObject = LoadSingleRun();
    CalculateMetrics calculateMetricsObject = CalculateMetrics();

    std::unique_ptr<CreateAndPlotGraphsBase> graphs;
    if (m_config_run_mode == "Any") {
        graphs = std::make_unique<CreateAndPlotGraphsMultipleRunModes>();
    } else {
        graphs = std::make_unique<CreateAndPlotGraphsOneRunMode>();
    }
    // main loop over all runs
    for (auto aRun : fileNamesAsPairs) {
        std::unique_ptr<NTupleInformation> nTupleInformationOneRun = std::make_unique
                < NTupleInformation > (std::get<0>(aRun));
        loadSingleRunObject.extractNtupleInformation(std::get<0>(aRun), std::get<2>(aRun),
                m_config_tank_tree_name, m_config_MRD_tree_name, m_config_verbosity, nTupleInformationOneRun);
        std::unique_ptr<RunMetrics> runMetricsOneRun = std::make_unique < RunMetrics
                > (std::get<0>(aRun), std::get<1>(aRun));
        calculateMetricsObject.calculateTankCharge(nTupleInformationOneRun, runMetricsOneRun);
        calculateMetricsObject.calculateMRDMetrics(nTupleInformationOneRun, runMetricsOneRun);
        runMetricsOneRun->printTankCharge();
        graphs->setTankCharge(runMetricsOneRun);
        graphs->setMRDMetrics(runMetricsOneRun);
    }
    graphs->drawAndSave(m_config_save_plots_as_pictures, m_config_picture_save_directory,
            m_config_root_file_name, m_config_save_plots_as_file);

    return true;
}

bool OfflineDataQualityMonitoring::Finalise() {

    return true;
}

std::vector<std::tuple<int, RunMode, std::vector<std::string> > > OfflineDataQualityMonitoring::loadListFile() {
    LoadFileList loadFileListObject = LoadFileList();
    if (m_config_run_mode == "Any") {
        loadFileListObject.retrieveFileNamesFromListFileForAllRunModes(
                m_config_list_file_name,
                m_config_verbosity);
    } else {
        loadFileListObject.retrieveFileNamesFromListFileWithMatchingRunMode(
                m_config_list_file_name,
                m_config_verbosity, m_config_run_mode);
    }
    loadFileListObject.assignFilesToRun(m_config_verbosity,
            m_config_run_number_prefix, m_config_run_number_suffix);
    if (m_config_verbosity > 1) {
        loadFileListObject.printFileNamesAssignedToRunNumber();
    }
    return loadFileListObject.getFileNamesAsTuples();
}

void OfflineDataQualityMonitoring::loadAndCheckConfiguration(const std::string &t_configfileName) {
    if (t_configfileName != "") {
        m_variables.Initialise(t_configfileName); // loading config file
    }
    m_variables.Print();
    m_variables.Get("filename", m_config_list_file_name);
    m_variables.Get("prefix", m_config_run_number_prefix);
    m_variables.Get("suffix", m_config_run_number_suffix);
    m_variables.Get("runMode", m_config_run_mode);
    m_variables.Get("triggerTreeName", m_config_trigger_tree_name);
    m_variables.Get("tankTreeName", m_config_tank_tree_name);
    m_variables.Get("mrdTreeName", m_config_MRD_tree_name);
    m_variables.Get("lappdTreeName", m_config_LAPPD_tree_name);
    m_variables.Get("verbosity", m_config_verbosity);
    m_variables.Get("numberOfRunsPerPoint", m_config_number_of_runs_per_point);
    m_variables.Get("saveHistogramsAsPictures", m_config_save_plots_as_pictures);
    m_variables.Get("pictureDirectory", m_config_picture_save_directory);
    m_variables.Get("saveHistogramsInROOTFile", m_config_save_plots_as_file);
    m_variables.Get("outputROOTFileName", m_config_root_file_name);
    checkConfigurationVariables();
    if (m_config_verbosity) {
        m_variables.Print();
    }
}

void OfflineDataQualityMonitoring::checkConfigurationVariables() {
    checkRunMode();
    checkNumberOfRunsPerPoint();
    checkAndAdjustVerbosity();
    checkSaveConfig();
}

void OfflineDataQualityMonitoring::checkRunMode() const {
    if (std::find(m_run_modes.begin(), m_run_modes.end(), m_config_run_mode) != m_run_modes.end()) {
        return;
    }
    //If no, print an error message and stop the program
    std::string logMessage { "Unexpected run mode provided! " };
    logMessage.append(m_config_run_mode);
    logMessage.append(" was provided, candidates are ");
    for (size_t iMode = 0; iMode < m_run_modes.size(); iMode++) {
        logMessage.append(m_run_modes.at(iMode));
        if (iMode != m_run_modes.size() - 1) {
            logMessage.append(", ");
        }
    }
    throw std::invalid_argument(logMessage);
}

void OfflineDataQualityMonitoring::checkNumberOfRunsPerPoint() const {
    if (m_config_number_of_runs_per_point < 1) {
        throw std::invalid_argument("Number of runs per point is smaller than one. Abort.");
    }
}

void OfflineDataQualityMonitoring::checkAndAdjustVerbosity() {
    if (m_config_verbosity < 0) {
        m_config_verbosity = 0;
        std::cout << "Verbosity was negative and is set to 0.\n";
        return;
    }
    if (m_config_verbosity > 2) {
        m_config_verbosity = 2;
        std::cout << "Verbosity was above 2 and is set to 2.n";
    }
}

void OfflineDataQualityMonitoring::checkSaveConfig() const {
    if (!m_config_save_plots_as_pictures and !m_config_save_plots_as_file) {
        throw std::invalid_argument("No output mode was specified as true.");
    }
    // ToDo: Does not work under C++17
//    std::filesystem::directory_entry entry { m_config_picture_save_directory };
//    if (!entry.exists()) {
//        throw std::invalid_argument("Specified output directory for pictures does not exist.");
//    }
    TFile outputROOTFile(m_config_root_file_name.c_str(), "RECREATE");
    if (!outputROOTFile.IsOpen()) {
        std::string logMessage { "Output file " + m_config_root_file_name
                + " couldn't be opened. Pictures are not saved!" };
        throw std::invalid_argument(logMessage);
    }
    outputROOTFile.Close();

}


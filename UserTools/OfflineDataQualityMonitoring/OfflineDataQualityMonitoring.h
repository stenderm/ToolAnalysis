#ifndef OfflineDataQualityMonitoring_H
#define OfflineDataQualityMonitoring_H

#include <string>
#include <iostream>
#include "RunMode.h"

#include "Tool.h"

/**
 * \class OfflineDataQualityMonitoring
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
 *
 * $Author: B.Richards $
 * $Date: 2019/05/28 10:44:00 $
 * Contact: b.richards@qmul.ac.uk
 */
class OfflineDataQualityMonitoring: public Tool {

public:

    OfflineDataQualityMonitoring(); ///< Simple constructor
    bool Initialise(std::string configfile, DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
    bool Execute(); ///< Execute function used to perform Tool purpose.
    bool Finalise(); ///< Finalise function used to clean up resources.

private:
    ///Configuration Variables
    // Input Mode, if files created with PhaseIITreeMaker (Old) or files created with the ANNIEEventTreeMaker (New) is expected
    std::string m_config_input_mode { };
    // Name of the list file to load for PhaseIITreeMaker files (old)
    std::string m_config_list_file_name_old { };
    // Name of the list file to load for ANNIEEventTreeMaker files (new)
    std::string m_config_list_file_name_new { };
    // Prefix of the run number ToDo: Needs adjustment probably for newer files
    std::string m_config_run_number_prefix { };
    // Suffix of the run number
    std::string m_config_run_number_suffix { };
    // Run mode of the runs that are to be analysed
    std::string m_config_run_mode { };
    // Trigger data tree name
    std::string m_config_trigger_tree_name { };
    // Tank data tree name
    std::string m_config_tank_tree_name { };
    // MRD data tree name
    std::string m_config_MRD_tree_name { };
    // LAPPD data tree name
    std::string m_config_LAPPD_tree_name { };
    // Verbosity
    int m_config_verbosity { };
    // Number of runs that are grouped together as data point
    int m_config_number_of_runs_per_point { };
    // Should the plots be saved as picture (.pdf)
    bool m_config_save_plots_as_pictures { };
    // Name of the directory, to which the pictures are saved
    std::string m_config_picture_save_directory { };
    // Should the plots be saved as a root file
    bool m_config_save_plots_as_file { };
    // Name of output root file
    std::string m_config_root_file_name { };


    // vector of available run modes, that is filled in constructor
    std::vector<std::string> m_run_modes { };

    void loadAndCheckConfiguration(const std::string& t_configfileName);
    // checks validity of configuration variables, exceptions are thrown in the individual methods if parameter is invalid
    void checkConfigurationVariables();
    // checks if the run mode from the configuration file exists in the vector of run modes
    void checkRunMode() const;
    // check if the input mode from the configuration file is one of Old, New or Both
    void checkInputMode() const;
    // checks if the number of runs per point from the configuration file is bigger or equal to one, throws exception if not
    void checkNumberOfRunsPerPoint() const;
    // checks if the verbosity is 0, 1 or 2, if it is higher, set to 2, if it is lower to 0
    void checkAndAdjustVerbosity();
    // checks the validity of the save paths and modes
    void checkSaveConfig() const;

    std::vector<std::tuple<int, RunMode, std::vector<std::string> > > loadListFileOld();

    std::vector<std::tuple<int, RunMode, std::vector<std::string> > > loadListFileNew();
};

#endif

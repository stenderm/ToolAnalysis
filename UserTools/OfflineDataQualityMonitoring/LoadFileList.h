/*
 * LoadFileList.hh
 *
 *  Created on: Sep 9, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_LOADFILELIST_HH_
#define INCLUDE_LOADFILELIST_HH_

#include <string>
#include <vector>

#include "RunMode.h"

/*
 * This class loads the name of the files that are contained in a list file, and assigns the file names to the according run number.
 */
class LoadFileList {
public:
    LoadFileList();
    virtual ~LoadFileList();
    // Retrieve the .root file names from the list file. This list will be unordered.
    void retrieveFileNamesFromListFileWithMatchingRunMode(const std::string& t_fileListName, int t_verbosity, const std::string& t_runMode);
    // Retrieve the .root file names from the list file for all possible modes. This list will be unordered.
    void retrieveFileNamesFromListFileForAllRunModes(const std::string& t_fileListName, int t_verbosity);
    // Assigns the file names from the list of file names to the corresponding runs, since some runs are in multiple files
    void assignFilesToRun(int t_verbosity, const std::string &t_prefix,
                          const std::string &t_suffix);
    // Prints the run number and then the files that are associated with that number
    void printFileNamesAssignedToRunNumber();
    // Getter for the vector of pairs<run number, vector of file names corresponding to run number>
    std::vector<std::tuple<int, RunMode, std::vector<std::string> > > getFileNamesAsTuples(){return m_file_names_assigned_run_number;}
    // Retrieves the files for the newer format from the file list
    void retrieveFileNamesFromNewListFile(const std::string& t_fileListName, int t_verbosity);
private:
    // A vector of pairs of filename and run mode
    std::vector<std::pair<std::string, RunMode> > m_list_of_file_names { };
    // Vector of file tuples consisting of run number, run mode and the file names that correspond to said run number
    std::vector<std::tuple<int, RunMode, std::vector<std::string> > > m_file_names_assigned_run_number { };
    // Extracts the run number from the file name
    int extractRunNumberFromFileName(const std::string& t_fileName, const std::string& t_prefix, const std::string& t_suffix, int t_verbosity);
};

#endif /* INCLUDE_LOADFILELIST_HH_ */

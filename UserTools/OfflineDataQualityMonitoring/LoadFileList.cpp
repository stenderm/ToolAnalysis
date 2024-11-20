/*
 * LoadFileList.cc
 *
 *  Created on: Sep 9, 2024
 *      Author: stenderm
 */

#include "LoadFileList.h"
#include <fstream>
#include <iostream>
#include <algorithm>

LoadFileList::LoadFileList() {
    // do nothing

}

LoadFileList::~LoadFileList() {
    // do nothing
}

//ToDo: Here is some duplicate code!
void LoadFileList::retrieveFileNamesFromListFileWithMatchingRunMode(std::string t_fileListName,
                                                                    int t_verbosity,
                                                                    std::string t_runMode) {
    m_list_of_file_names.clear();
    std::fstream listFile;
    listFile.open(t_fileListName, std::ios::in);
    if (!listFile.is_open()) {
        throw std::invalid_argument("Cannot open list file!");
    }
    std::string aLine;
    if (t_verbosity > 1) {
        std::cout << "\n";
        std::cout << "Loaded file names: \n";
    }
    while (getline(listFile, aLine)) {
        if (t_verbosity > 1) {
            std::cout << aLine << "\n";
        }
        // check if specified run mode is contained in file name
        // ensure that the path is ignored and only the file name is looked at
        size_t positionInStringSeparator { aLine.find_last_of("/") };
        size_t positionInStringKeyword { aLine.find(t_runMode, positionInStringSeparator) };
        if (positionInStringKeyword != std::string::npos) {
            // put only file names into the list that contain the run mode keyword together with the run mode, which is redundant in this case
            m_list_of_file_names.push_back(std::make_pair(aLine, getRunModeFromString(t_runMode)));
        }
    }
    listFile.close();
}

void LoadFileList::retrieveFileNamesFromListFileForAllRunModes(
        const std::string &t_fileListName, int t_verbosity) {
    m_list_of_file_names.clear();
    std::fstream listFile;
    listFile.open(t_fileListName, std::ios::in);
    if (!listFile.is_open()) {
        throw std::invalid_argument("Cannot open list file!");
    }
    std::string aLine;
    if (t_verbosity > 1) {
        std::cout << "\n";
        std::cout << "Loaded file names: \n";
    }
    while (getline(listFile, aLine)) {
        if (t_verbosity > 1) {
            std::cout << aLine << "\n";
        }
        // check if specified run mode is contained in file name
        // ensure that the path is ignored and only the file name is looked at
        size_t positionInStringSeparator { aLine.find_last_of("/") };
        for (const auto &aModeName : runModeToString) {
            // Obviously it is necessary that no run keyword is part of another keyword for this to work
            size_t positionInStringKeyword { aLine.find(aModeName.second, positionInStringSeparator) };
            if (positionInStringKeyword != std::string::npos) {
                // put only file names into the list that contain the run mode keyword
                m_list_of_file_names.push_back(std::make_pair(aLine, aModeName.first));
            }
        }
    }
    listFile.close();
}

void LoadFileList::assignFilesToRun(int t_verbosity, const std::string &t_prefix,
                                    const std::string &t_suffix) {
    // By default sort works on the first element of the pair, which is ideal here
    std::sort(m_list_of_file_names.begin(), m_list_of_file_names.end());
    if (t_verbosity > 1) {
        std::cout << "\n";
        std::cout << "Sorted file Names: \n";
    }
    for (const std::pair<std::string, RunMode>& oneFile : m_list_of_file_names) {
        int runNumber = extractRunNumberFromFileName(oneFile.first, t_prefix, t_suffix, t_verbosity);

        // see if run number already exists in list
        auto it = std::find_if(m_file_names_assigned_run_number.begin(),
                m_file_names_assigned_run_number.end(),
                [&runNumber](const std::tuple<int, RunMode, std::vector<std::string> > & element) {
                    return std::get<0>(element) == runNumber;
                });
        // element is not found, new pair has to be created and pushed into the vector
        if (it == m_file_names_assigned_run_number.end()) {
            std::vector<std::string> tempFileName { };
            tempFileName.push_back(oneFile.first);
            m_file_names_assigned_run_number.push_back( std::make_tuple(runNumber, oneFile.second, tempFileName) ); // @suppress("Function cannot be instantiated") // @suppress("Invalid arguments")
        } else {
            // element is found and the file name can be pushed to the vector inside the pair
            std::get<2>(*it).push_back(oneFile.first);
        }

        if (t_verbosity > 1) {
            std::cout << oneFile.first << " with extracted run number " << runNumber << "\n";
        }
    }
}

int LoadFileList::extractRunNumberFromFileName(const std::string &t_fileName,
                                               const std::string &t_prefix,
                                               const std::string &t_suffix, int t_verbosity) {
    int sizePrefix = t_prefix.size();
    size_t positionPrefix { t_fileName.find(t_prefix) + sizePrefix };
    // check if prefix is found at all
    if (positionPrefix == std::string::npos) {
        if (t_verbosity) {
            std::cout << "Prefix " << t_prefix
                    << " not found to extract the run number. Skip the file.\n";
            return -1;
        }
    }
    //check if prefix is not unique in path
    size_t secondPositionPrefix { t_fileName.find(t_prefix, positionPrefix + sizePrefix) };
    if (secondPositionPrefix != std::string::npos) {
        if (t_verbosity) {
            std::cout << "Prefix " << t_prefix << " was not unique in file name. Skip the file.\n";
            return -1;
        }
    }

    size_t positionSuffix { t_fileName.find(t_suffix, positionPrefix) };
    // check if suffix is found at all
    if (positionSuffix == std::string::npos) {
        if (t_verbosity) {
            std::cout << "Suffix " << t_suffix
                    << " not found to extract the run number. Skip the file.\n";
            return -1;
        }
    }
    std::string runNumberString = t_fileName.substr(positionPrefix,
            positionSuffix - positionPrefix);
    int runNumber = std::stoi(runNumberString);
    return runNumber;
}

void LoadFileList::printFileNamesAssignedToRunNumber() {
    std::cout << "\n";
    for (auto aRun : m_file_names_assigned_run_number) {
        std::cout << "Run number >>>>> " << std::get<0>(aRun) << " <<<<< consists of the files:\n";
        for (auto aFileName : std::get<2>(aRun)) {
            std::cout << aFileName << "\n";
        }
        std::cout << "\n";
    }

}


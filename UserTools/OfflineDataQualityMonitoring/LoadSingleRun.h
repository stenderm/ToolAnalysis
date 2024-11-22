/*
 * LoadSingleRun.hh
 *
 *  Created on: Sep 10, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_LOADSINGLERUN_HH_
#define INCLUDE_LOADSINGLERUN_HH_

#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "NTupleInformation.h"

/*
 *  This class loads the ntuple data of a single run, that can then be processed further
 */
class LoadSingleRun {
public:
    LoadSingleRun();
    virtual ~LoadSingleRun();
    void extractNtupleInformation(int t_runNumber, const std::vector<std::string> &t_fileNames,
                                  const std::string& t_tankTreeName, const std::string& t_MRDTreeName, int t_verbosity,
                                  std::unique_ptr<NTupleInformation>& t_ntupleInformationOneRun);
private:
    // checks if trigger tree exists and loads the trigger ntuple information into a NTupleInformation object
    void loadInformationTrigger(std::shared_ptr<TFile> t_inputFile, std::string t_triggerTreeName,
                                int t_verbosity, std::unique_ptr<NTupleInformation>& t_ntupleInformationOneRun);
    // checks if MRD tree exists and loads the MRD ntuple information into a NTupleInformation object
    void loadInformationMRD(std::shared_ptr<TFile> t_inputFile, std::string t_MRDTreeName,
                            int t_verbosity, std::unique_ptr<NTupleInformation>& t_ntupleInformationOneRun);
    // checks if tank tree exists and loads the tank ntuple information into a NTupleInformation object
    void loadInformationTank(std::shared_ptr<TFile> t_inputFile, std::string t_tankTreeName,
                             int t_verbosity, std::unique_ptr<NTupleInformation>& t_ntupleInformationOneRun);
    // checks if LAPPD tree exists and loads the LAPPD ntuple information into a NTupleInformation object
    //ToDo: Find run with LAPPD information to implement this!
    void loadInformationLAPPD(std::shared_ptr<TFile> t_inputFile, std::string t_LAPPDTreeName,
                              int t_verbosity, std::unique_ptr<NTupleInformation>& t_ntupleInformationOneRun);
    // open tree
    TTree* openTree(std::shared_ptr<TFile> t_inputFile, std::string treeName, int t_verbosity);
    //
    void addOnePropertyToVector(int t_entry, int t_cluster, double t_propertyToAdd, std::vector<std::vector<double> > & t_vectorToAddTo);

};

#endif /* INCLUDE_LOADSINGLERUN_HH_ */

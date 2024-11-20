/*
 * RunModes.hh
 *
 *  Created on: Oct 14, 2024
 *      Author: stenderm
 */

#ifndef INCLUDE_RUNMODE_HH_
#define INCLUDE_RUNMODE_HH_

#include <string>
#include <unordered_map>
#include <algorithm>
#include <optional>

enum class RunMode {
    Beam, Cosmic, LED, AmBe, Laser
};

const std::unordered_map<RunMode, std::string> runModeToString { { RunMode::Beam, "Beam" }, {
        RunMode::Cosmic, "Cosmic" }, { RunMode::LED, "LED" }, { RunMode::AmBe, "AmBe" }, {
        RunMode::Laser, "Laser" } };

inline std::string getRunModeName(RunMode t_runMode) {
    auto it = runModeToString.find(t_runMode);
    return (it != runModeToString.end()) ? it->second : "Unknown";
}

inline RunMode getRunModeFromString(const std::string &t_runMode) {
    auto it = std::find_if(runModeToString.begin(), runModeToString.end(),
            [&t_runMode](const std::pair<RunMode, std::string> &element) {
                return element.second == t_runMode;
            });
    if (it == runModeToString.end()) {
       throw std::invalid_argument("Color not found!");
    }
    return it->first;  // Return the found color
}

#endif /* INCLUDE_RUNMODE_HH_ */

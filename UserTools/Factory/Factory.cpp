#include "Factory.h"

Tool* Factory(std::string tool) {
    Tool* ret = 0;
    if (tool=="PlotWaveforms") ret=new PlotWaveforms;
    if (tool=="OfflineDataQualityMonitoring") ret=new OfflineDataQualityMonitoring;

    return ret;
}

# TopoRecoConverter

TopoRecoConverter
This tool is used to convert ANNIE data to the format for the Topological Track Reconstruction.
Tools to have upstream from this tool for it to work properly:
LoadWCSim (Loads the regular PMT data)
LoadWCSimLAPPD (Loads the LAPPD data)
MCParticleProperties (Provides map for pdg to mass)
Maybe also MCRecoEventLoader
myDigitBuilder

myDigitBuilder
## Data

Describe any data formats TopoRecoConverter creates, destroys, changes, or analyzes. E.G.

**RawLAPPDData** `map<Geometry, vector<Waveform<double>>>`
* Takes this data from the `ANNIEEvent` store and finds the number of peaks


## Configuration

Describe any configuration variables for TopoRecoConverter.

```
param1 value1
param2 value2
```

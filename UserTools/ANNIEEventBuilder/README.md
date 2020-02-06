# ANNIEEventBuilder

ANNIEEventBuilder

The ANNIEEventBuilder takes parsed PMT data from the DataDecoder tool,
matches it to trigger information in the RawData TriggerData store, and
constructs an ANNIEEvent BoostStore.

## Data

Describe any data formats ANNIEEventBuilder creates, destroys, changes, or analyzes. E.G.

**RawLAPPDData** `map<Geometry, vector<Waveform<double>>>`
* Takes this data from the `ANNIEEvent` store and finds the number of peaks


## Configuration

Describe any configuration variables for ANNIEEventBuilder.

```
param1 value1
param2 value2
```
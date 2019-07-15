#include "LoadGeometry.h"

LoadGeometry::LoadGeometry():Tool(),adet(nullptr){}


bool LoadGeometry::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Make the ANNIEEvent Store if it doesn't exist
  int recoeventexists = m_data->Stores.count("ANNIEEvent");
  if(recoeventexists==0) m_data->Stores["ANNIEEvent"] = new BoostStore(false,2);

  m_variables.Get("verbosity", verbosity);
  m_variables.Get("FACCMRDGeoFile", fFACCMRDGeoFile);
  m_variables.Get("TankPMTGeoFile", fTankPMTGeoFile);
  m_variables.Get("LAPPDGeoFile", fLAPPDGeoFile);
  m_variables.Get("DetectorGeoFile", fDetectorGeoFile);

  //Check files exist
  if(!this->FileExists(fDetectorGeoFile)){
		Log("LoadGeometry Tool: File for Detector Geometry does not exist!",v_error,verbosity);
        std::cout << "Filepath was... " << fDetectorGeoFile << std::endl;
		return false;
  }
  if(!this->FileExists(fFACCMRDGeoFile)){
		Log("LoadGeometry Tool: File for FACC/MRD Geometry does not exist!",v_error,verbosity);
        std::cout << "Filepath was... " << fFACCMRDGeoFile << std::endl;
		return false;
  }
  if(!this->FileExists(fLAPPDGeoFile)){
    Log("LoadGeometry Tool: File for the LAPPDs does not exist!",v_error,verbosity);
        std::cout << "Filepath was... " << fDetectorGeoFile << std::endl;
    return false;
  }

  //Initialize the geometry using the geometry CSV file entries
  this->InitializeGeometry();

  //Load MRD Geometry Detector/Channel Information
  this->LoadFACCMRDDetectors();

  this->LoadLAPPDs();

  m_data->Stores.at("ANNIEEvent")->Header->Set("AnnieGeometry",AnnieGeometry,true);

  //AnnieGeometry->GetChannel(0); // trigger InitChannelMap

  return true;
}


bool LoadGeometry::Execute(){
  return true;
}


bool LoadGeometry::Finalise(){
  std::cout << "LoadGeometry tool exitting" << std::endl;
  return true;
}


void LoadGeometry::InitializeGeometry(){
  Log("LoadGeometry tool: Now loading DetectorGeoFile and initializing geometry",v_message,verbosity);
  //Get the Detector file data key
  std::string DetectorLegend = this->GetLegendLine(fDetectorGeoFile);
  std::vector<std::string> DetectorLegendEntries;
  boost::split(DetectorLegendEntries,DetectorLegend, boost::is_any_of(","), boost::token_compress_on);

  //Initialize at zero; will be set later after channels are loaded
  int numtankpmts = 0;
  int numlappds = 0;
  int nummrdpmts = 0;
  int numvetopmts = 0;

  //Initialize data that will be fed to Geometry (units in meters)
  int geometry_version;
  double tank_xcenter,tank_ycenter,tank_zcenter;
  double tank_radius,tank_halfheight, pmt_enclosed_radius, pmt_enclosed_halfheight;
  double mrd_width,mrd_height,mrd_depth,mrd_start;

  std::string line;
  ifstream myfile(fDetectorGeoFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#") != std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector data (should only be one line here)
    while(getline(myfile,line)){
      std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> DataEntries;
      boost::split(DataEntries,line, boost::is_any_of(","), boost::token_compress_on);
      for (int i=0; i<DataEntries.size(); i++){
        //Check Legend at i, load correct data type
        int ivalue;
        double dvalue;
        if(DetectorLegendEntries.at(i) == "geometry_version") ivalue = std::stoi(DataEntries.at(i));
        else dvalue = std::stod(DataEntries.at(i));
        if (DetectorLegendEntries.at(i) == "geometry_version") geometry_version = ivalue;
        if (DetectorLegendEntries.at(i) == "tank_xcenter") tank_xcenter = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_ycenter") tank_ycenter = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_zcenter") tank_zcenter = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_radius") tank_radius = dvalue;
        if (DetectorLegendEntries.at(i) == "tank_halfheight") tank_halfheight = dvalue;
        if (DetectorLegendEntries.at(i) == "pmt_enclosed_radius") pmt_enclosed_radius = dvalue;
        if (DetectorLegendEntries.at(i) == "pmt_enclosed_halfheight") pmt_enclosed_halfheight = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_width") mrd_width = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_height") mrd_height = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_depth") mrd_depth = dvalue;
        if (DetectorLegendEntries.at(i) == "mrd_start") mrd_start = dvalue;
      }
    }
    Position tank_center(tank_xcenter, tank_ycenter, tank_zcenter);
    // Initialize the Geometry
    AnnieGeometry = new Geometry(geometry_version,
                                 tank_center,
                                 tank_radius,
                                 tank_halfheight,
                                 pmt_enclosed_radius,
                                 pmt_enclosed_halfheight,
                                 mrd_width,
                                 mrd_height,
                                 mrd_depth,
                                 mrd_start,
                                 numtankpmts,
                                 nummrdpmts,
                                 numvetopmts,
                                 numlappds,
                                 geostatus::FULLY_OPERATIONAL);
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  myfile.close();
}

void LoadGeometry::LoadFACCMRDDetectors(){
  //First, get the MRD file data key
  Log("LoadGeometry tool: Now loading FACC/MRD detectors",v_message,verbosity);
  std::string MRDLegend = this->GetLegendLine(fFACCMRDGeoFile);
  std::vector<std::string> MRDLegendEntries;
  boost::split(MRDLegendEntries,MRDLegend, boost::is_any_of(","), boost::token_compress_on);

  std::string line;
  ifstream myfile(fFACCMRDGeoFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector specs
    while(getline(myfile,line)){
      std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> SpecLine;
      boost::split(SpecLine,line, boost::is_any_of(","), boost::token_compress_on);
      if(verbosity>4) std::cout << "This line of data: " << line << std::endl;
      //Parse data line, make corresponding detector/channel
      bool add_ok = this->ParseMRDDataEntry(SpecLine,MRDLegendEntries);
      if(not add_ok){
        std::cerr<<"Faild to add Detector to Geometry!"<<std::endl;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  if(myfile.is_open()) myfile.close();
    Log("LoadGeometry tool: FACC/MRD Detector/Channel loading complete",v_message,verbosity);
}

bool LoadGeometry::ParseMRDDataEntry(std::vector<std::string> SpecLine,
        std::vector<std::string> MRDLegendEntries){
  //Parse the line for information needed to fill the detector & channel classes
  int detector_num,channel_num,detector_system,orientation,layer,side,num,
      rack,TDC_slot,TDC_channel,discrim_slot,discrim_ch,
      patch_panel_row,patch_panel_col,amp_slot,amp_channel,
      hv_crate,hv_slot,hv_channel,nominal_HV,polarity;
  double x_center,y_center,z_center,x_width,y_width,z_width;
  std::string PMT_type,cable_label,paddle_label;

  //Search for Legend entry.  Fill value type if found.
  Log("LoadGeometry tool: parsing data line into variables",v_debug,verbosity);
  for (int i=0; i<SpecLine.size(); i++){
    int ivalue;
    double dvalue;
    std::string svalue;
    for (int j=0; j<MRDIntegerValues.size(); j++){
      if(MRDLegendEntries.at(i) == MRDIntegerValues.at(j)){
        ivalue = std::stoi(SpecLine.at(i));
        break;
      }
    }
    for (int j=0; j<MRDStringValues.size(); j++){
      if(MRDLegendEntries.at(i) == MRDStringValues.at(j)){
        svalue = SpecLine.at(i);
        break;
      }
    }
    for (int j=0; j<MRDDoubleValues.size(); j++){
      if(MRDLegendEntries.at(i) == MRDDoubleValues.at(j)){
        dvalue = std::stod(SpecLine.at(i));
        break;
      }
    }
    //Integers
    if (MRDLegendEntries.at(i) == "detector_num") detector_num = ivalue;
    if (MRDLegendEntries.at(i) == "channel_num") channel_num = ivalue;
    if (MRDLegendEntries.at(i) == "detector_system") detector_system = ivalue;
    if (MRDLegendEntries.at(i) == "orientation") orientation = ivalue;
    if (MRDLegendEntries.at(i) == "layer") layer = ivalue;
    if (MRDLegendEntries.at(i) == "side") side = ivalue;
    if (MRDLegendEntries.at(i) == "num") num = ivalue;
    if (MRDLegendEntries.at(i) == "rack") rack = ivalue;
    if (MRDLegendEntries.at(i) == "TDC_slot") TDC_slot = ivalue;
    if (MRDLegendEntries.at(i) == "TDC_channel") TDC_channel = ivalue;
    if (MRDLegendEntries.at(i) == "discrim_slot") discrim_slot = ivalue;
    if (MRDLegendEntries.at(i) == "discrim_ch") discrim_ch = ivalue;
    if (MRDLegendEntries.at(i) == "patch_panel_row") patch_panel_row = ivalue;
    if (MRDLegendEntries.at(i) == "patch_panel_col") patch_panel_col = ivalue;
    if (MRDLegendEntries.at(i) == "amp_slot") amp_slot = ivalue;
    if (MRDLegendEntries.at(i) == "amp_channel") amp_channel = ivalue;
    if (MRDLegendEntries.at(i) == "hv_crate") hv_crate = ivalue;
    if (MRDLegendEntries.at(i) == "hv_slot") hv_slot = ivalue;
    if (MRDLegendEntries.at(i) == "hv_channel") hv_channel = ivalue;
    if (MRDLegendEntries.at(i) == "nominal_HV") nominal_HV = ivalue;
    if (MRDLegendEntries.at(i) == "polarity") polarity = ivalue;
    //Doubles
    if (MRDLegendEntries.at(i) == "x_center") x_center = dvalue;
    if (MRDLegendEntries.at(i) == "y_center") y_center = dvalue;
    if (MRDLegendEntries.at(i) == "z_center") z_center = dvalue;
    if (MRDLegendEntries.at(i) == "x_width") x_width = dvalue;
    if (MRDLegendEntries.at(i) == "y_width") y_width = dvalue;
    if (MRDLegendEntries.at(i) == "z_width") z_width = dvalue;
    //Strings
    if (MRDLegendEntries.at(i) == "PMT_type") PMT_type = svalue;
    if (MRDLegendEntries.at(i) == "paddle_label") paddle_label = svalue;
    if (MRDLegendEntries.at(i) == "cable_label") cable_label = svalue;
  }

  //FIXME Need the direction of the MRD PMT
  //FIXME Do we want the Paddle's center position?  Or PMT?
  //FIXME: things that are not loaded in with the default det/channel format:
  //  - detector_system (fed as "MRD"), orientation, layer, side, num
  //  - discrim_slot, discrim_ch
  //  - patch_panel_row, patch_panel_col, amp_slot, amp_channel
  //  - nominal_HV, polarity
  //  - cable_label, paddle_label
  //  - x_width, y_width, z_width
  //
  if(verbosity>4) std::cout << "Filling a FACC/MRD data line into Detector/Channel classes" << std::endl;
  Detector adet(detector_num,
                "MRD",
                "MRD", //Change to orientation for PaddleDetector class?
                Position( x_center/100.,
                          y_center/100.,
                          z_center/100.),
                Direction(0.,
                          0.,
                          0.),
                PMT_type,
                detectorstatus::ON,
                0.);

  int MRD_x, MRD_y, MRD_z;
  // orientation 0=horizontal, 1=vertical
  MRD_x = (orientation) ? num  : side;
  MRD_y = (orientation) ? side : num;
  // veto layers are both cabled as z=0, with the layers differentiated by x=0, x=1
  // in practice of course, both span the same x, but are offset in z.
  if(layer>0) MRD_z = layer;
  else        MRD_z = side;

  Paddle apad( MRD_x,
               MRD_y,
               MRD_z,
               orientation,
               Position( x_center/100.,
                         y_center/100.,
                         z_center/100.),
               std::pair<double,double>{x_center-(x_width/200.), x_center+(x_width/200.)},
               std::pair<double,double>{y_center-(y_width/200.), y_center+(y_width/200.)},
               std::pair<double,double>{z_center-(z_width/200.), z_center+(z_width/200.)});

  Channel pmtchannel( channel_num,
                      Position(0,0,0.),
                      -1, // stripside
                      -1, // stripnum
                      rack,
                      TDC_slot,
                      TDC_channel,
                      -1,                 // TDC has no level 2 signal handling
                      -1,
                      -1,
                      hv_crate,
                      hv_slot,
                      hv_channel,
                      channelstatus::ON);

  // Add this channel to the geometry
  if(verbosity>4) cout<<"Adding channel "<<channel_num<<" to detector "<<detector_num<<endl;
  adet.AddChannel(pmtchannel);
  if(verbosity>5) cout<<"Adding detector to Geometry"<<endl;
  AnnieGeometry->AddDetector(adet);
  if(verbosity>4) cout<<"Adding paddle to Geometry"<<endl;
  AnnieGeometry->SetDetectorPaddle(detector_num, apad);
  return true;
}

void LoadGeometry::LoadLAPPDs(){
  //First, get the LAPPD file data key
  Log("LoadGeometry tool: Now loading LAPPDs",v_message,verbosity);
  std::string LAPPDLegend = this->GetLegendLine(fLAPPDGeoFile);
  std::vector<std::string> LAPPDLegendEntries;
  boost::split(LAPPDLegendEntries,LAPPDLegend, boost::is_any_of(","), boost::token_compress_on);

  std::string line;
  ifstream myfile(fLAPPDGeoFile.c_str());
  if (myfile.is_open()){
    //First, get to where data starts
    while(getline(myfile,line)){
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataStartLineLabel)!=std::string::npos) break;
    }
    //Loop over lines, collect all detector specs
    detector_num_store = 100000;
    counter = 0;
    while(getline(myfile,line)){
      std::cout << line << std::endl; //has our stuff;
      if(line.find("#")!=std::string::npos) continue;
      if(line.find(DataEndLineLabel)!=std::string::npos) break;
      std::vector<std::string> SpecLine;
      boost::split(SpecLine,line, boost::is_any_of(","), boost::token_compress_on);
      if(verbosity>4) std::cout << "This line of data: " << line << std::endl;
      //Parse data line, make corresponding detector/channel
      bool add_ok = this->ParseLAPPDDataEntry(SpecLine,LAPPDLegendEntries);
      if(not add_ok){
        std::cerr<<"Faild to add Detector to Geometry!"<<std::endl;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  if(myfile.is_open()) myfile.close();
    Log("LoadGeometry tool: LAPPD Detector/Channel loading complete",v_message,verbosity);
}


bool LoadGeometry::ParseLAPPDDataEntry(std::vector<std::string> SpecLine,
        std::vector<std::string> LAPPDLegendEntries){
  //Parse the line for information needed to fill the detector & channel classes
   int detector_num,channel_strip_side,channel_strip_num;
   unsigned int channel_signal_crate,channel_signal_card,channel_signal_channel,channel_level2_crate,channel_level2_card,channel_level2_channel,channel_hv_crate,channel_hv_card,channel_hv_channel,channel_num;
   double detector_position_x,detector_position_y,detector_position_z,detector_direction_x,detector_direction_y,detector_direction_z,channel_position_x,channel_position_y,channel_position_z;
   std::string detector_type,detector_status,channel_status;
  //Search for Legend entry.  Fill value type if found.
  Log("LoadGeometry tool: parsing data line into variables",v_debug,verbosity);
  for (int i=0; i<SpecLine.size(); i++){
    int ivalue;
    unsigned int uivalue;
    double dvalue;
    std::string svalue;
    for (int j=0; j<LAPPDIntegerValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDIntegerValues.at(j)){
        ivalue = std::stoi(SpecLine.at(i));
        break;
      }
    }
    for (int j=0; j<LAPPDStringValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDStringValues.at(j)){
        svalue = SpecLine.at(i);
        break;
      }
    }
    for (int j=0; j<LAPPDDoubleValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDDoubleValues.at(j)){
        dvalue = std::stod(SpecLine.at(i));
        break;
      }
    }
    for (int j=0; j<LAPPDUnIntValues.size(); j++){
      if(LAPPDLegendEntries.at(i) == LAPPDUnIntValues.at(j)){
        uivalue = std::stoul(SpecLine.at(i));
        break;
      }
    }
    //Integers
    if (LAPPDLegendEntries.at(i) == "detector_num") detector_num = ivalue;
    if (LAPPDLegendEntries.at(i) == "channel_strip_side") channel_strip_side = ivalue;
    if (LAPPDLegendEntries.at(i) == "channel_strip_num") channel_strip_num = ivalue;

    //Unsigned Integers
    if (LAPPDLegendEntries.at(i) == "channel_signal_crate") channel_signal_crate = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_signal_card") channel_signal_card = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_signal_channel") channel_signal_channel = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_level2_crate") channel_level2_crate = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_level2_card") channel_level2_card = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_level2_channel") channel_level2_channel = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_hv_crate") channel_hv_crate = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_hv_card") channel_hv_card = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_hv_channel") channel_hv_channel = uivalue;
    if (LAPPDLegendEntries.at(i) == "channel_num") channel_num = uivalue;

    //Doubles
    if (LAPPDLegendEntries.at(i) == "detector_position_x") detector_position_x = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_position_y") detector_position_y = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_position_z") detector_position_z = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_direction_x") detector_direction_x = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_direction_y") detector_direction_y = dvalue;
    if (LAPPDLegendEntries.at(i) == "detector_direction_z") detector_direction_z = dvalue;
    if (LAPPDLegendEntries.at(i) == "channel_position_x") channel_position_x = dvalue;
    if (LAPPDLegendEntries.at(i) == "channel_position_y") channel_position_y = dvalue;
    if (LAPPDLegendEntries.at(i) == "channel_position_z") channel_position_z = dvalue;

    //Strings
    if (LAPPDLegendEntries.at(i) == "detector_type") detector_type = svalue;
    if (LAPPDLegendEntries.at(i) == "detector_status") detector_status = svalue;
    if (LAPPDLegendEntries.at(i) == "channel_status") channel_status = svalue;
  }

  if(verbosity>4) std::cout << "Filling a LAPPD data line into Detector/Channel classes" << std::endl;
  if(detector_num != detector_num_store){
  detectorstatus detstat;
  if(detector_status == "OFF"){
    detstat = detectorstatus::OFF;
    }
    else if(detector_status == "ON"){
      detstat = detectorstatus::ON;
    }
    else if(detector_status == "UNSTABLE"){
      detstat = detectorstatus::UNSTABLE;
    }
    else{
      std::cerr << "The chosen detector status isn't available!!!" << std::endl;
    }
  //TODO Somewhere it has to be stated that the units are in [m] for LAPPDs for now
  adet = new Detector(464+detector_num,
                "LAPPD",
                "Barrel",
                Position(detector_position_x,
                        detector_position_y,
                        detector_position_z),
                Direction(detector_direction_x,
                          detector_direction_y,
                          detector_direction_z),
                detector_type,
                detstat,
                0.);
  detector_num_store = detector_num;
  }

  channelstatus channelstat;
  if(channel_status == "OFF"){
      channelstat = channelstatus::OFF;
      }
  else if(channel_status == "ON"){
      channelstat = channelstatus::ON;
        }
  else if(channel_status == "UNSTABLE"){
      channelstat = channelstatus::UNSTABLE;
      }
  else{
  std::cerr << "The chosen channel status isn't available!!!" << std::endl;
      }
  Channel lappdchannel(464+channel_num,
                      Position(channel_position_x,
                               channel_position_y,
                               channel_position_z),
                      channel_strip_side,
                      channel_strip_num,
                      channel_signal_crate,
                      channel_signal_card,
                      channel_signal_channel,
                      channel_level2_crate,
                      channel_level2_card,
                      channel_level2_channel,
                      channel_hv_crate,
                      channel_hv_card,
                      channel_hv_channel,
                      channelstat);

  // Add this channel to the detector
  if(adet != nullptr){
  if(verbosity>4) cout<<"Adding channel "<<channel_num<<" to LAPPD "<<detector_num<<endl;
  adet->AddChannel(lappdchannel);
  }
  counter++;
  if(adet != nullptr && counter == 60){
  if(verbosity>5) cout<<"Adding LAPPD to Geometry"<<endl;
  AnnieGeometry->AddDetector(*adet);
  counter = 0;
  }
  return true;
}

bool LoadGeometry::FileExists(std::string name) {
  ifstream myfile(name.c_str());
  return myfile.good();
}


std::string LoadGeometry::GetLegendLine(std::string name) {
  if(verbosity>4) std::cout << "Getting legend of file: " << name << std::endl;
  std::string line;
  std::string legendline = "null";
  ifstream myfile(name.c_str());
  if (myfile.is_open()){
    while(std::getline(myfile,line)){
      if(verbosity>4) std::cout << line << std::endl;
      if(line.find("#") != std::string::npos) continue;
      if(line.find(LegendLineLabel) != std::string::npos){
        //Next line is the title line
        getline(myfile,line);
        legendline = line;
        if(verbosity>4) std::cout<<"Legend line loaded. Legend is: " << legendline << std::endl;
        break;
      }
    }
  } else {
    Log("LoadGeometry tool: Something went wrong opening a file!!!",v_error,verbosity);
  }
  if(legendline=="null"){
    Log("LoadGeometry tool: Legend line label not found!!!",v_error,verbosity);
  }
  myfile.close();
  return legendline;
}

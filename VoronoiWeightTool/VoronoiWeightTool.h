#ifndef VoronoiWeightTool_VoronoiWeightTool_H
#define VoronoiWeightTool_VoronoiWeightTool_H

#include <EventLoop/Algorithm.h>
#include "EventLoopAlgs/NTupleSvc.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "TRegexp.h"
#include <fstream>
#include <iostream>

#ifndef __CINT__
#include "xAODCore/ShallowCopy.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODJet/JetConstituentVector.h"
#endif

#include "xAODTracking/VertexContainer.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "PATInterfaces/SystematicRegistry.h"

#include <TTree.h>

//ASG:
#include "AsgTools/AsgTool.h"

//class VoronoiWeights : public EL::Algorithm
class VoronoiWeightTool : virtual public asg::AsgTool 
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
    std::string m_inputContainer = "CaloCalTopoClusters",
                m_outputContainer = "VoronoiClusters";

    bool m_debug = false;
    bool m_doLC = false;
    bool m_doSpread = true;
    int m_nSigma = 1;
    std::string m_eventInfo      = "EventInfo",
                m_jets          = "AntiKt4LCTopoJets";

   // Constructor with parameters:
   VoronoiWeightTool(const std::string& name);

   // Destructor:
   virtual ~VoronoiWeightTool();
   

   // methods used in the analysis
   virtual StatusCode MakeVoronoiClusters(std::vector< std::pair<fastjet::PseudoJet,std::vector<float> > >&);
   void SpreadPt(std::vector< std::pair< fastjet::PseudoJet,std::vector<float> > >& correctedptvec, float spreadr=0.4, float alpha=2);

private:
  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store;  //!
  struct PJcomp;
  std::vector<fastjet::PseudoJet> clusters; //!

public:
  // this is a standard constructor
  VoronoiWeightTool ();

  // these are the functions inherited from AsgTool
  virtual StatusCode initialize ();
  virtual StatusCode execute ();
  virtual StatusCode finalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(VoronoiWeightTool, 1);
};

#endif

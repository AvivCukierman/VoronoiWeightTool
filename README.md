# VoronoiWeightTool
This tool applies Voronoi weights to clusters. It makes a shallow copy container of the clusters with the pT of the clusters modified appropriately.
## Dependencies
This package makes use of [UChicago](https://github.com/UCATLAS)'s [xAODAnaHelpers](https://github.com/UCATLAS/xAODAnaHelpers) package and [Giordon Stark](https://github.com/kratsg)'s [xAODJetReclustering](https://github.com/kratsg/xAODJetReclustering) package.

## Installing
To install,
```bash
mkdir myRootCore && cd $_
rcSetup Base,2.3.41
git clone https://github.com/kratsg/xAODJetReclustering
git clone https://github.com/UCATLAS/xAODAnaHelpers
git clone https://github.com/AvivCukierman/VoronoiWeightTool
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJet/tags
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJetContrib/tags
rc find_packages
rc compile
```


## Configurations

 Property           | Type                      | Default                   | Description
:-------------------|:-------------------------:|--------------------------:|:-------------------------------------------------------------------------------------
InputContainer   | string                    |    "CaloCalTopoClusters"                       | name of the input cluster container to be modified
OutputContainer  | string                    |    "VoronoiClusters"                       | name of the output shallow copy cluster container with modified pT
doLCWeights       | bool                     | false                      | apply LC weights to clusters before Voronoi algorithms
doSpread  | bool     | true | after Voronoi subtraction, apply spreading to negative pT clusters
nSigma     | int                    | 0                      | at the end, suppress all clusters with pT<nSigma * sigmaRho * sqrt(area)

N.B. The `doLCWeights` option indicates whether to apply LC weights (if true) or EM weights (if false) to clusters. The output clusters have their pTs modified *in the given state*. I.e. with `doLCWeights` set to true, only the LC pT of the new cluster collection is modified, and with `doLCWeights` set to false, only the EM pT of the new cluster collection is modified. In order to set the state of the output (or any) cluster collection, the `CaloClusterChangeSignalStateList` tool is used, as follows:
```c++
#include "xAODCaloEvent/CaloClusterChangeSignalState.h"
```

```c++
  const xAOD::CaloClusterContainer*             new_clusters   (nullptr);
  if(evtStore()->retrieve(new_clusters,m_OutputContainer).isFailure()) Error(APP_NAME,"Could not retrieve the Voronoi cluster container");

  CaloClusterChangeSignalStateList stateHelperList;
  for(auto clust: *new_clusters){
    if(m_doLCWeights) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    //now clust->pt() will give Voronoi modified pT, and jet reconstruction will use Voronoi modified pT
  }
```

## Using
Add this package as a dependency in `cmt/Makefile.RootCore`.

Add a header
```c++
#include <VoronoiWeightTool/VoronoiWeightTool.h>
```

Set up the tool in the `initialize()` portion of your algorithm as a pointer

```c++
  m_VoronoiTool = new VoronoiWeightTool(m_name);
  m_VoronoiTool->setProperty("InputContainer", m_InputContainer);
  m_VoronoiTool->setProperty("OutputContainer", m_OutputContainer);
  m_VoronoiTool->setProperty("doSpread", m_doSpread);
  m_VoronoiTool->setProperty("doLCWeights", m_doLCWeights);
  m_VoronoiTool->setProperty("nSigma", m_nSigma);
  m_VoronoiTool->initialize();
```

And then simply call `m_jetReclusteringTool->execute()` in the `execute()` portion of your algorithm to fill the TStore with the appropriate container(s). Don't forget to delete the pointer when you're done.
```c++
if(m_VoronoiTool) delete m_VoronoiTool;
```

You can then do jet clustering as normal on the output cluster container.

Note that as it behaves like an `AsgTool`, the functions `setProperty()` and `initialize()` have a return type `StatusCode`.

## Authors
- [Aviv Cukierman](https://github.com/AvivCukierman)

## Acknowledgements
[Giordon Stark](https://github.com/kratsg), from whom I plagiarized this readme and also whose code I used for inspiration.

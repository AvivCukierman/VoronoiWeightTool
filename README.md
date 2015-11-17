# VoronoiWeightTool
This tool applies Voronoi weights to clusters. It makes a shallow copy container of the clusters with the pT of the clusters modified appropriately.
## Dependencies
This package makes use of [UChicago](https://github.com/UCATLAS)'s [xAODAnaHelpers](https://github.com/UCATLAS/xAODAnaHelpers) package and [Giordon Stark](https://github.com/kratsg)'s [xAODJetReclustering](https://github.com/kratsg/xAODJetReclustering) package.

## Installing
To install,
```bash
mkdir myRootCore && cd $_
rcSetup Base,2.3.23
git clone https://github.com/kratsg/xAODJetReclustering
git clone https://github.com/UCATLAS/xAODAnaHelpers
git clone https://github.com/AvivCukierman/VoronoiWeightTool
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJet/tags
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJetContrib/tags
rc find_packages
rc compile
```
(I think, I haven't tested this yet)

## Configurations

 Property           | Type                      | Default                   | Description
:-------------------|:-------------------------:|--------------------------:|:-------------------------------------------------------------------------------------
InputContainer   | string                    |    "CaloCalTopoClusters"                       | name of the input cluster container to be modified
OutputContainer  | string                    |    "VoronoiClusters"                       | name of the output shallow copy cluster container with modified pT
doLCWeights       | bool                     | false                      | apply LC weights to clusters before Voronoi algorithms
doSpread  | bool     | true | after Voronoi subtraction, apply spreading to negative pT clusters
nSigma     | int                    | 0                      | at the end, suppress all clusters with pT<nSigma * sigmaRho * sqrt(area)

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

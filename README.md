# Nature-Communications-Junctions
Code from McEvoy et al. (2022) Nature Communications

*****************************************************
This repository contains the following files:

1. A MATLAB (last tested on ver 2020b) file to simulate remodeling of a tri-cellular junction:
"Nat_comms_vertex.m"

2. A COMSOL Multiphysics (last tested on ver5.4a) finite element model to support supplementary analytical model of tri-cell junction stress:
"Vertex_FE.mph"

* All files and software were tested in Windows 10.
** Code commented within "Nat_comms_vertex.m" for detailed explanation of analysis.

*****************************************************
Matlab code demo:
1. Install MATLAB (https://www.mathworks.com/products/matlab.html). The installation usually takes approximately 1 hour.
2. Open "Nat_comms_vertex.m" and click "Run".
* New windows will then open with the simulated results

Note: This file is configured to reproduce simulated junction dynamics at cell vertices [Figure 4].
All other results may be reproduced my modifying parameters, as outlined in the manuscript. 

*****************************************************
COMSOL simulation demo:
1. Install COMSOL Multiphysics (https://www.comsol.com/). The installation usually takes approximately 1.5 hours.
2. Open "Vertex_FE.mph".
3. Navigate to "Settings" [middle panel] and click on "Compute".
* Model may be used to obtain FEM reults within Figure S2.

*****************************************************
* To ensure accurate reproduction of results we recommend you reopen the unedited files provided before modifying parameters and settings.
* The demos may take 10 minutes to run depending on the machine. 

*****************************************************
For any questions regarding these files, please contact Eoin McEvoy (eoin.mcevoy@universityofgalway.ie).

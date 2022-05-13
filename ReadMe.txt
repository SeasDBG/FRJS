FRJS Code GitHub
---Code Base for the Analysis done for Focused Rotary Jet Spinning ---
5/13/22


--PIV ipynotebooks are used to generate particle image velocimetry maps from ventricle ejection fractions --
How to use:
-1) .tif files should be placed in the Raw data folder, using the file structure Raw\Condition\Sample\File0001.tif
Data analysis scripts will extend to work for any number of samples and any number of conditions, as long as the file tree follows that architecture

-2) Run the sample ipynotebook labelled "1 - VentriclePIV.ipynb", this will perform PIV analysis on .tif files placed in the above raw folders, and will save that data analysis into the "Analyzed" Folder

-3) Run the sample ipynotebook labelled "2 -FrameAverage.ipynb", this will perform phase averaging over the PIV images, and will sum across the basal plan to calculate net Flux in the verticle direction
(Note, processing scripts for producing ejection fractions and cardiac output comparisons are coded for comparing 2 conditions explicity. To compare more, adjusts need to be made to the code, in the regions noted)

--Deformation Mapping -- Will produce deformation maps over the surface of the ventricle
How to use:
-1) .tif files should be placed in the Raw data folder, using the file structure Raw\Condition\Sample\File0001.tif
Data analysis scripts will extend to work for any number of samples and any number of conditions, as long as the file tree follows that architecture

-2) Run the sample ipynotebook labelled "Ventricle Deformation Mapping.ipynb", this will perform PIV analysis on .tif files placed in the above raw folders, and will save that data analysis into the "Analyzed" Folder, along with graphic representation of the deformation.


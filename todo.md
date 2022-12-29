
# TODO list

## Planned Structure

+tools 
    .load - main method, which calls the following (unless suppressed):
    .detectSpikes - support threshold and TEO techniques
    .segmentEpochs - code to parse visual stimulus info from analog-in channel
    .sort - spike sorting interface (PCA-based)

    <possibly removeCommonMode>

+analysis
    .f012
    .linear [pca]

+plots
    .epochs
    .heatmap
    .spikeMetrics




## Done

### intan RHD load code

intan RHD load code (for Intan2000 files) doesn't appear to work on the intan file saved by dario
unclear what the error is, the internal format does not appear to line up correctly

the xml file settings.xml has the following: 

FileFormat="Traditional" 
Filename.BaseFilename="050922" 
Filename.Path="C:/Users/.../MEA/bReaChes/2022/September" 

A little bit of inspection is revealing that the per-channel structre is internally different, will need to chase down more documentation once I'm back on the ground. 

... (13/09/22) have downloaded the 3.0 version of the Intan code from the Intan website, will make parallel changes. 
... (07/10/22) have finished refactor, edited the output style to be more similar to tools.readBlackrock (different repo)
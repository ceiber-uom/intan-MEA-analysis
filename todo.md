
# TODO list

[ ] implement analysis.f012 (for stim epoch data)
[ ] implement spike correlation analysis (TBD)
[ ] implement epoching code for spikes
[ ] finish implementing spike detection code
[ ] implement forChannels code for waves
[ ] confirm channel ordering

# Revised Plan (Dec 2022)

## +tools

the +tools is a collection of command-line tools for loading data and data processing, as well as some internal tools which significantly simplify the plots and analysis code. 

### Loading data

[x] `.readIntan` - load 'raw' wave data
[x] `.readPlexon` - load sorted spike data processed in Plexon and saved in a tabular format
[x] `.readExpoXML` - load stimulus data from Expo

### Processing data

[x] `.segmentEpochs` - split waves/spikes based on photocell epoch trigger
[x] `.removeCommonMode` - per EIBER 2015
[x] `.removeLineNoise` - per PIETERSEN EIBER 2015
[x] `.TEO` (Teager energy operater) - per EIBER 2015
[ ] `.spikeDetection` - detect spikes from wave (automated). This was deprioritised because this was being done in Plexon/Kilosort

### Internal Tools 

[x] `.forChannels` - boilerplate loop over channels/units
[x] `.forWaveType` - boilerplate loop over Intan AMP inputs
[x] `.peakseek` - alternative to findpeaks

## +plots

These are the core functions for visualusation of the loaded data (with minimal processing). There are also some plot-specific utilities (plots.tidy, plots.PDF_tools)

### 'top-level' plot functions

[x] `.epochs` - basic plot of recorded signal waves
[x] `.raster` - basic plot of recorded signal spike times
[x] `.heatmap` - heatmap based on channel layout. If Y not specified, shows heatmap of channel spikerate 
[x] `.ISI` - plot of inter-spike-interval (default: histograms)
[x] `.spikeshape` - plot of spike shapes 

### plot utilities 

[x] `.tidy` - implement better axis style (tidyPlotForIllustrator)
[x] `.layout` - source of truth for channel subplot layout
[x] 


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
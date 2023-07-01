# intan-MEA-analysis
Basic analysis code for Intan 64-channel MEA system

https://intantech.com/


This is an update to some of the code in the appendices of https://doi.org/10.26190/unsworks/18556

Includes an update to https://au.mathworks.com/matlabcentral/fileexchange/43662-read-intan-technologies-files 


## Troubleshooting common (expected) errors:

If you are using tools.simplify and trying to plot a single channel/unit and you get an error which does not appear when plotting multiple channel/units, try adding the '--undo' flag to the call to the plot code. This error will frequently appear as "The logical indices contain a true value outside of the array bounds."

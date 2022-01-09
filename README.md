# GW_Rank
 A GW Followup Pipeline Designed for CFHT MEGACam and WIRCam

## Setup:

You will need to build the 'gwemopt' tool, follow the readme instructions in that folder for setup.
Currently, to get all the components of the code to be located on the python path I need to add the following to the .rc file for my local shell:
    export PYTHONPATH=$PYTHONPATH:$HOME/path/to/GW_Rank
    export PYTHONPATH=$PYTHONPATH:$HOME/path/to/GW_Rank/GW_LIGO_gcn
    export PYTHONPATH=$PYTHONPATH:$HOME/path/to/GW_Rank/HOGWARTs
    export PYTHONPATH=$PYTHONPATH:$HOME/path/to/GW_Rank/gwemopt
where $HOME/path/to/GW_Rank is the path to wherever you save the GW_Rank repo locally. In the future I will fix this so you don't need to manually add the paths.

## Directories:

*'GW_LIGO_gcn' contains the LIGO_obs.py script, which will listen for a gcn alert and, once an alert is found, execute the galaxy ranking and tile ranking codes. 
NOTE: the default code is set up to take a skymap you already downloaded as input; comment out the last 3 lines of code as needed if you just want to listen for gcn alerts vs testing w/skymaps. You can see a list of previous signals w/skymaps for LIGO O3 here: https://gracedb.ligo.org/superevents/public/O3/.
NOTE: pay attention to whether you are looking at alerts with the 'role' = 'test' or 'observation' ('observation' represents real signals, 'test' alerts go out ~ every hour.)
*'gwemopt' contains the tile ranking code, adapted from https://github.com/mcoughlin/gwemopt. This is relevant if using MEGACam.
*'HOGWARTs' contains the galaxy ranking code, adapted from https://github.com/Lanasalmon/HOGWARTs/tree/master. This is relevant if using WIRCam.
*After running the code and receiving a gcn alert, the skymap corresponding to that alert will go into a created directory called 'skymaps'.
*After the information from the skymap has been processed by the galaxy and tile ranking codes, the resulting ranked lists and associated information will go into a created directory called 'outputs'.

## Notes:

The 'gwemopt' and 'HOGWARTs' codes contain their own README files, and 'gwemopt' also has installation documentation. Please refer to those for more instructions.
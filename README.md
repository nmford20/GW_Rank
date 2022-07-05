# GW_Rank
 A GW Followup Pipeline Designed for CFHT MEGACam and WIRCam

## Setup:

You will need to build the `gwemopt` tool, follow the readme instructions in that folder for setup.

NOTE: List of package dependencies will be added in the future.

## Directories:

- `GW_LIGO_gcn` contains the LIGO_obs.py script, which will listen for a gcn alert and, once an alert is found, execute the galaxy ranking and tile ranking codes. 

NOTE: the default code is set up to take a skymap you already downloaded as input; comment out the last 3 lines of code as needed if you just want to listen for gcn alerts vs testing w/skymaps. You can see a list of previous signals w/skymaps for LIGO O3 here: https://gracedb.ligo.org/superevents/public/O3/.

NOTE: pay attention to whether you are looking at gcn alerts with the `role` classified as `test` or `observation` (`observation` represents real signals, `test` alerts are generated ~ every hour.)

- `gwemopt` contains the tile ranking code, adapted from https://github.com/mcoughlin/gwemopt. This is relevant if using MEGACam.
- `HOGWARTs` contains the galaxy ranking code, adapted from https://github.com/Lanasalmon/HOGWARTs/tree/master. This is relevant if using WIRCam.
- After running the code and receiving a gcn alert, the skymap corresponding to that alert will go into a created directory called `skymaps`.
- After the information from the skymap has been processed by the galaxy and tile ranking codes, the resulting ranked lists and associated information will go into a created directory called `outputs`.

## Notes:

The `gwemopt` and `HOGWARTs` codes contain their own README files, and `gwemopt` also has installation documentation. I've also created more thorough documentation in the included pdf file. Please refer to those for more instructions.

## Acknowledgements:

This code relies heavily on the publicly available [gwemopt](https://github.com/mcoughlin/gwemopt) and [HOGWARTs](https://github.com/Lanasalmon/HOGWARTs/tree/master) codes. This code was developed under the supervision of Dr. Daryl Haggard at McGill University and Dr. John Ruan at Bishops University.

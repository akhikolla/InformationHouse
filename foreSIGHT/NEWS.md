foreSIGHT 1.0.0
----------------------------------------------------------------
- Major revision
- Restructured code to split existing functions into smaller functions and use updated function arguments
- The new functions are createExpSpace, generateScenario, generateScenarios, getSimSummary
- Added function runSystemModel to run system models and plotting functions named plot*
- Added helper functions named view* and writeControlFile
- Updated vignettes


foreSIGHT 0.9.8
----------------------------------------------------------------
Updated to work with latest version of cowplot.


foreSIGHT 0.9.6
----------------------------------------------------------------
Updated plotLayers, performanceSpaces and quickSpace to maintain compatability with update of ggplot2 (3.0.0).
New weather generator functionality added including ability for user to specify parameter bounds to scenarioGenerator function.
New features added. modSimulator and modCalibrator functions added so users can calibrate and simulate from weather generators directly.
New modelTag "Radn-har26-wgen" and relevant attributes added.
Bug fixes: attHold and attPertrurb arguments now able to be specified in any order (mixing hydroclimate variable type orders now ok).

foreSIGHT 0.9.2
----------------------------------------------------------------
Longer description provided with method references (inc. doi's).

foreSIGHT 0.9.1
----------------------------------------------------------------
Prebuilt vignette index included

foreSIGHT 0.9.0
----------------------------------------------------------------
First complete implementation ready for CRAN submission

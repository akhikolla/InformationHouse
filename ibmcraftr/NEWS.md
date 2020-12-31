ibmcraftr 1.0.0
===============

The new version 0.2.0 has C++ code integration to make the functions run faster. It has also a higher level
function to actually run the transitions for the number of timesteps that users specify. 

- New and faster state_transition function in C++, called `stRCPP`. 
- `stateT` is the C++ specific function which is called from RCPP.
- `run_state_trans` function now organizes how to run the state transtions
in the specified number of timesteps.
- Seperate chunks of functions to handle: 
    1. Rate to probability transformation
    2. Calculation of cummulative probability
        * Previously, cummulative probability calculation was written within state_trans and stRCPP.
        * Now this is replaced with cumprob function which is written specifically to do that. 

ibmcraftr 0.1.0
===============

This is the very first version of my package. It has 2 functionalities for now: 

- To populate a synthetic population with states
- To make state transitions between those states

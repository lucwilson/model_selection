# model_selection
Software respository for model selection spectral parameterization (ms-specparam) and generating analyses from the Principled model-selection algorithm for parameterizing oscillatory peaks in neural power spectra preprint.

ms-specparam is available in both MATLAB and Python programming languages, as well as in the [Brainstorm distribution](https://neuroimage.usc.edu/brainstorm/Introduction) (Tadel et al., 2011) where users can interact with the algorithm using a graphical user interface.

# Information

Code for running the ms-specparam and ms-SPRiNT algorithms can be found in ms_specparam.m and ms_SPRiNT.m, respectively. These algorithms are also available from the [Brainstorm distribution](https://neuroimage.usc.edu/brainstorm/Introduction) (Tadel et al., 2011).

Code used for simulating spectra in the preprint, as well as analyzing outputs from ms-specparam, can be found in ms_specparam_code_analysis.m.

Simulated power spectra (and their constituent parameters) and noise matrix used in the preprnt can be found in ms_specparam_data_spectra.m.

Code used for analyzing Cam-CAN data in the preprint can be found in the figures folder.

Empirical metadata and results in the preprint (Cam-CAN), where permissible to share, can be found in /empirical_data.

Please contact the authors of this package for any questions and suggestions.

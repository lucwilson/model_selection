ms-specparam
=======================

This folder contains the necessary code to perform `ms-specparam` in MATLAB.

If you wish to use `ms-specparam` in MATLAB with a graphical user interface, it is also included as part of the free [Brainstorm distribution](https://neuroimage.usc.edu/brainstorm/Introduction) (Tadel et al., 2011).


Dependencies
------------

`ms-specparam` requires the user to have the [Optimization toolbox](https://www.mathworks.com/products/optimization.html) installed.

Usage
-----

`ms-specparam` can be loaded and used in your MATLAB analysis pipeline by downloading the ms-specparam folder and adding it to your system path in MATLAB.

    addpath('path_to/ms-specparam')

For a worked example, see the tutorial.

Tutorial
--------

A full demonstration of `ms-specparam`, contrasted against `specparam`, can be found in the tutorial folder as a MATLAB script.

Additional Resources
--------------------

If you are interested in using `ms-specparam` for parsimonious time-resolved spectral modelling, see `ms_SPRiNT`.

License and Reference
---------------------

This code is distributed under an Apache 2.0 License.

If you use this code in your project, please cite the ms-specparam preprint:

Wilson, L. E., da Silva Castanheira, J., Lévesque Kinder, B., & Baillet, S. (2024). _A Bayesian Model-Selection Approach for Determining the Number of Spectral Peaks in Neural Power Spectra. bioRxiv._ https://doi.org/10.1101/2024.08.01.606216

You must also cite the original specparam manuscript:

Donoghue T, Haller M, Peterson EJ, Varma P, Sebastian P, Gao R, Noto T, Lara AH, Wallis JD,
Knight RT, Shestyuk A, & Voytek B (2020). _Parameterizing neural power spectra into periodic
and aperiodic components. Nature Neuroscience, 23_, 1655-1665.
DOI: [10.1038/s41593-020-00744-x](https://www.nature.com/articles/s41593-020-00744-x)

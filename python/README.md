FOOOF-ms (ms-specparam)
=======================

`fooof-ms` is the Python equivalent to ms-specparam, which was originally written in MATLAB. It produces near-identical results to the MATLAB version (where differences in versions are due to language-specific optimization algorithms).


Dependencies
------------

`fooof-ms` is written in Python. It is a subclass of [FOOOF](https://fooof-tools.github.io), and thus requires all the same dependencies including Python >= 3.7 to run.

Below are the required dependencies for FOOOF:

- `numpy <https://github.com/numpy/numpy>`
- `scipy <https://github.com/scipy/scipy>` >= 0.19

There are also optional dependencies, which offer additional functionality:

- `matplotlib <https://github.com/matplotlib/matplotlib>` for visualizing data and model fits
- `tqdm <https://github.com/tqdm/tqdm>` for printing progress bars when fitting many models
- `pandas <https://github.com/pandas-dev/pandas>` for exporting model fit results to dataframes
- `pytest <https://github.com/pytest-dev/pytest>` to run the test suite locally.

Usage
-----

`fooof-ms` can be loaded and used in your Python pipeline by downloading the fooof-ms folder, adding it to your system path in Python and importing it alongside FOOOF.

    sys.path.append('path_to/fooof-ms')
    from fooof import FOOOF
    from fooof_ms import FOOOF_MS

For a worked example, see the tutorial.

Tutorial
--------

A tutorial for using fooof-ms is provided in plain Python (.py) and Jupyter notebook (.ipynb) formats, which can be found in the tutorial folder.


License and Reference
---------------------

This code is distributed under an Apache 2.0 License.

If you use this code in your project, please cite the ms-specparam preprint:

Wilson, L. E., da Silva Castanheira, J., Lévesque Kinder, B., & Baillet, S. (2024). Model selection for spectral parameterization. bioRxiv. https://doi.org/10.1101/2024.08.01.606216

You must also cite the original specparam manuscript:

Donoghue T, Haller M, Peterson EJ, Varma P, Sebastian P, Gao R, Noto T, Lara AH, Wallis JD,
Knight RT, Shestyuk A, & Voytek B (2020). Parameterizing neural power spectra into periodic
and aperiodic components. Nature Neuroscience, 23, 1655-1665.
DOI: 10.1038/s41593-020-00744-x


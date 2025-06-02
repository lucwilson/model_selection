"""
FOOOF-MS: A FOOOF subclass for model selection using the Bayesian Information Criterion.
"""

from collections import namedtuple
import os
import json
import numpy as np
from numpy.linalg import LinAlgError
from scipy.optimize import curve_fit
from fooof import FOOOF
from fooof.core.funcs import get_ap_func, gaussian_function
from fooof.core.io import fpath, fname
from fooof.core.modutils import copy_doc_func_to_method
from fooof.core.strings import _format
from fooof.core.utils import dict_array_to_lst, dict_select_keys, dict_lst_to_array
from fooof.sim.gen import gen_aperiodic, gen_periodic
from fooof.utils.params import compute_gauss_std



def compute_bic(y_true, y_pred, n_params):
    """Compute Bayesian Information Criterion (BIC)."""
    n = len(y_true)
    mse = np.mean((y_true - y_pred) ** 2)
    loglik = -n / 2 * (1 + np.log(mse) + np.log(2 * np.pi))
    return n_params * np.log(n) - 2 * loglik


def model_func(freqs, *model_params):
    if len(model_params) % 3 == 0:
        ap_length = 3
        ap_fit = gen_aperiodic(freqs, model_params[0:ap_length], 'knee')
    else:
        ap_length = 2
        ap_fit = gen_aperiodic(freqs, model_params[0:ap_length], 'fixed')
    peak_fit = np.zeros_like(freqs)
    for p in range(ap_length,len(model_params),3):
        peak_fit += model_params[p+1] * np.exp(-0.5 * ((freqs - model_params[p]) / model_params[p+2]) ** 2)
    return ap_fit + peak_fit


class FOOOF_MS(FOOOF):
    """FOOOF model that optimizes to a parsimonious model via BIC."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bic_data_ = {}

    def fit(self, freqs = None, power_spectrum = None, freq_range=None, *args, **kwargs):
        # If freqs & power_spectrum provided together, add data to object.
        if freqs is not None and power_spectrum is not None:
            self.add_data(freqs, power_spectrum, freq_range)
        # If power spectrum provided alone, add to object, and use existing frequency data
        #   Note: be careful passing in power_spectrum data like this:
        #     It assumes the power_spectrum is already logged, with correct freq_range
        elif isinstance(power_spectrum, np.ndarray):
            self.power_spectrum = power_spectrum

        # Check that data is available
        if not self.has_data:
            raise NoDataError("No data available to fit, can not proceed.")

        # Check and warn about width limits (if in verbose mode)
        if self.verbose:
            self._check_width_limits()

        best_params = None
        best_bic = np.inf
        bic_values = []

        # Step 1: Robust aperiodic fit on original spectrum
        self.aperiodic_params_ = self._robust_ap_fit(self.freqs, self.power_spectrum)
        self._ap_fit = gen_aperiodic(self.freqs, self.aperiodic_params_)

        # Step 2: Flatten the power spectrum using fit aperiodic fit
        self._spectrum_flat = self.power_spectrum - self._ap_fit

        # Step 3: Iteratively fit peaks and estimate BIC
        bic_values = [];

        # Step 3a: Find peaks, and fit them with gaussians
        self.gaussian_params_ = self._est_peaks(np.copy(self._spectrum_flat))
        gauss_params_mem = self.gaussian_params_[self.gaussian_params_[:,1].argsort()[::-1]]
        for pks in range(gauss_params_mem.shape[0]+1):
            self.gaussian_params_ = self._est_fit(gauss_params_mem[0:pks])
            # Step 3b: Calculate the peak fit
            #   Note: if no peaks are found, this creates a flat (all zero) peak fit
            self._peak_fit = gen_periodic(self.freqs, np.ndarray.flatten(self.gaussian_params_))

            # Step 3c: Create peak-removed (but not flattened) power spectrum
            self._spectrum_peak_rm = self.power_spectrum - self._peak_fit

            # Step 3d: Run aperiodic fit on peak-removed power spectrum
            #   This overwrites previous aperiodic fit, and recomputes the flattened spectrum
            self.aperiodic_params_ = self._simple_ap_fit(self.freqs, self._spectrum_peak_rm)

            # Step 3e: Setup initials and bounds for full model optimization
            self.peak_params_ = self._create_peak_params(self.gaussian_params_)
            model_guess = np.concatenate((self.aperiodic_params_, self.gaussian_params_.flatten()))
            lo_bound = [[peak[0] - 2 * self._cf_bound * peak[2], 0, self._gauss_std_limits[0]]
                        for peak in self.gaussian_params_]
            hi_bound = [[peak[0] + 2 * self._cf_bound * peak[2], np.inf, self._gauss_std_limits[1]]
                        for peak in self.gaussian_params_]
            lo_bound = [bound if bound[0] > self.freq_range[0] else \
                [self.freq_range[0], *bound[1:]] for bound in lo_bound]
            hi_bound = [bound if bound[0] < self.freq_range[1] else \
                [self.freq_range[1], *bound[1:]] for bound in hi_bound]
            gaus_param_bounds = (tuple([item for sublist in lo_bound for item in sublist]),
                                 tuple([item for sublist in hi_bound for item in sublist]))
            model_param_bounds = (tuple(self._ap_bounds[0]+gaus_param_bounds[0]),tuple(self._ap_bounds[1]+gaus_param_bounds[1]))

            # Step 3f: Attempt full model optimization
            try:
                model_params, _ = curve_fit(model_func, self.freqs, self.power_spectrum,
                                               p0=model_guess, maxfev=self._maxfev, bounds=model_param_bounds)
            except RuntimeError:
                raise FitError("Model fitting failed due to not finding "
                               "parameters in the peak component fit.")
            except LinAlgError:
                raise FitError("Model fitting failed due to a LinAlgError during peak fitting. "
                               "This can happen with settings that are too liberal, leading, "
                               "to a large number of guess peaks that cannot be fit together.")

            # Step 3g: Recover model features
            if len(self.aperiodic_params_) % 3 == 0:
                self.aperiodic_params_ = np.array(model_params[:3])
                self.gaussian_params_ = np.array(model_params[3:].reshape(-1,3) if self.gaussian_params_.shape[0] > 0 else [])
            else:
                self.aperiodic_params_ = np.array(model_params[:2])
                self.gaussian_params_ = np.array(model_params[2:].reshape(-1,3) if self.gaussian_params_.shape[0] > 0 else [])
            self._ap_fit = gen_aperiodic(self.freqs, self.aperiodic_params_)
            self._peak_fit = gen_periodic(self.freqs, np.ndarray.flatten(self.gaussian_params_))
            self.fooofed_spectrum_ = self._ap_fit + self._peak_fit

            # Step 3h: Estimate BIC, update best model if appropriate
            bic = compute_bic(self.power_spectrum, self.fooofed_spectrum_, len(model_params))
            bic_values.append(bic)

            if bic < best_bic:
                best_bic = bic
                best_ap = self.aperiodic_params_
                if self.gaussian_params_.shape[0] > 0:
                    best_pk = self.gaussian_params_[self.gaussian_params_[:, 0].argsort()]
                else:
                    best_pk = self.gaussian_params_
                best_spec = self.fooofed_spectrum_

        # Step 4: Save parameters and stats from best model
        self.aperiodic_params_ = best_ap
        self.gaussian_params_ = best_pk
        self.peak_params_ = self._create_peak_params(self.gaussian_params_)
        self.fooofed_spectrum_ = best_spec
        self._calc_r_squared()
        self._calc_error()
        self.bic_data_ = {
            'bic_values': bic_values,
            'peak_counts': list(range(gauss_params_mem.shape[0]+1)),
            'best_index': int(np.argmin(bic_values)),
            'n_models_tested': gauss_params_mem.shape[0]+1
        }

    def _est_peaks(self, flat_iter):
        """Iteratively fit peaks to flattened spectrum.

        Parameters
        ----------
        flat_iter : 1d array
            Flattened power spectrum values.

        Returns
        -------
        gaussian_params : 2d array
            Parameters that define the gaussian fit(s).
            Each row is a gaussian, as [mean, height, standard deviation].
        """

        # Initialize matrix of guess parameters for gaussian fitting
        guess = np.empty([0, 3])

        # Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian
        #   Stopping procedures: limit on # of peaks, or relative or absolute height thresholds
        while len(guess) < self.max_n_peaks:

            # Find candidate peak - the maximum point of the flattened spectrum
            max_ind = np.argmax(flat_iter)
            max_height = flat_iter[max_ind]

            # Stop searching for peaks once height drops below height threshold
            if max_height <= self.peak_threshold * np.std(flat_iter):
                break

            # Set the guess parameters for gaussian fitting, specifying the mean and height
            guess_freq = self.freqs[max_ind]
            guess_height = max_height

            # Halt fitting process if candidate peak drops below minimum height
            if not guess_height > self.min_peak_height:
                break

            # Data-driven first guess at standard deviation
            #   Find half height index on each side of the center frequency
            half_height = 0.5 * max_height
            le_ind = next((val for val in range(max_ind - 1, 0, -1)
                           if flat_iter[val] <= half_height), None)
            ri_ind = next((val for val in range(max_ind + 1, len(flat_iter), 1)
                           if flat_iter[val] <= half_height), None)
            # Guess bandwidth procedure: estimate the width of the peak
            try:
                # Get an estimated width from the shortest side of the peak
                #   We grab shortest to avoid estimating very large values from overlapping peaks
                # Grab the shortest side, ignoring a side if the half max was not found
                short_side = min([abs(ind - max_ind) \
                    for ind in [le_ind, ri_ind] if ind is not None])

                # Use the shortest side to estimate full-width, half max (converted to Hz)
                #   and use this to estimate that guess for gaussian standard deviation
                fwhm = short_side * 2 * self.freq_res
                guess_std = compute_gauss_std(fwhm)

            except ValueError:
                # This procedure can fail (extremely rarely), if both le & ri ind's end up as None
                #   In this case, default the guess to the average of the peak width limits
                guess_std = np.mean(self.peak_width_limits)

            # Check that guess value isn't outside preset limits - restrict if so
            #   Note: without this, curve_fitting fails if given guess > or < bounds
            if guess_std < self._gauss_std_limits[0]:
                guess_std = self._gauss_std_limits[0]
            if guess_std > self._gauss_std_limits[1]:
                guess_std = self._gauss_std_limits[1]

            # Collect guess parameters and subtract this guess gaussian from the data
            guess = np.vstack((guess, (guess_freq, guess_height, guess_std)))
            peak_gauss = gaussian_function(self.freqs, guess_freq, guess_height, guess_std)
            flat_iter = flat_iter - peak_gauss
        # Check peaks based on edges, and on overlap, dropping any that violate requirements

        guess = self._drop_peak_cf(guess)
        guess = self._drop_peak_overlap(guess)

        return guess

    def _est_fit(self, guess):

        if len(guess) > 0:
            gaussian_params = self._fit_peak_guess(guess)
            gaussian_params = gaussian_params[gaussian_params[:, 0].argsort()]
        else:
            gaussian_params = np.empty([0, 3])
        return gaussian_params


    def report(self, freqs=None, power_spectrum=None, freq_range=None, plt_log=False):
        """Run model fit, and display a report, which includes a plot, and printed results.

        Parameters
        ----------
        freqs : 1d array, optional
            Frequency values for the power spectrum.
        power_spectrum : 1d array, optional
            Power values, which must be input in linear space.
        freq_range : list of [float, float], optional
            Desired frequency range to fit the model to.
            If not provided, fits across the entire given range.
        plt_log : bool, optional, default: False
            Whether or not to plot the frequency axis in log space.

        Notes
        -----
        Data is optional, if data has already been added to the object.
        """
        self.fit(freqs, power_spectrum, freq_range)
        self.plot(plt_log=plt_log)
        self.print_results(concise=False)

    def add_results(self, fooof_result):
        """Add results data into object from a FOOOFResults object.

        Parameters
        ----------
        fooof_result : FOOOFResults
            A data object containing the results from fitting a FOOOF model.
        """
        self.aperiodic_params_ = fooof_result.aperiodic_params
        self.gaussian_params_ = fooof_result.gaussian_params
        self.peak_params_ = fooof_result.peak_params
        self.r_squared_ = fooof_result.r_squared
        self.error_ = fooof_result.error
        self.bic_data_ = fooof_result.bic_data

        self._check_loaded_results(fooof_result._asdict())

    def get_results(self):
        """Return model fit parameters and goodness of fit metrics.

        Returns
        -------
        FOOOFResults
            Object containing the model fit results from the current object.
        """
        return FOOOFResults(**{key.strip('_') : getattr(self, key) \
            for key in get_description()['results']})

    def print_results(self, concise=False):
        """Print out model fitting results.

        Parameters
        ----------
        concise : bool, optional, default: False
            Whether to print the report in a concise mode, or not.
        """
        print(gen_results_fm_str(self, concise))


    @property
    def bic_data(self):
        return self.bic_data_

    def report_bic(self):
        if not self.bic_data_:
            print("No MS data available. Run '.fit()' first.")
            return

        print("BIC Model Selection Report")
        print("---------------------------")
        for i, bic in enumerate(self.bic_data_['bic_values']):
            selected = " <--- selected" if i == self.bic_data_['best_index'] else ""
            print(f"{i} peaks: BIC = {bic:.2f}{selected}")

    def save(self, file_name, file_path=None, append=False,
             save_results=False, save_settings=False, save_data=False):

        save_fm(self, file_name, file_path, append, save_results, save_settings, save_data)

def get_description():
    """Get dictionary specifying FOOOF attributes, and what kind of data they store.

    Returns
    -------
    attributes : dict
        Mapping of FOOOF object attributes, and what kind of data they are.

    Notes
    -----
    This function organizes public FOOOF object attributes into:

    - results : parameters for and measures of the model
    - settings : model settings
    - data : input data
    - meta_data : meta data of the inputs
    - arrays : data stored in arrays
    - model_components : component pieces of the model
    - descriptors : descriptors of the object status and model results
    """

    attributes = {'results' : ['aperiodic_params_', 'gaussian_params_', 'peak_params_',
                               'r_squared_', 'error_', 'bic_data_'],
                  'settings' : ['peak_width_limits', 'max_n_peaks',
                                'min_peak_height', 'peak_threshold',
                                'aperiodic_mode'],
                  'data' : ['power_spectrum', 'freq_range', 'freq_res'],
                  'meta_data' : ['freq_range', 'freq_res'],
                  'arrays' : ['freqs', 'power_spectrum', 'aperiodic_params_',
                              'peak_params_', 'gaussian_params_'],
                  'model_components' : ['fooofed_spectrum_', '_spectrum_flat',
                                        '_spectrum_peak_rm', '_ap_fit', '_peak_fit'],
                  'descriptors' : ['has_data', 'has_model', 'n_peaks_']
                  }

    return attributes

def gen_results_fm_str(fm, concise=False):
    """Generate a string representation of model fit results.

    Parameters
    ----------
    fm : FOOOF
        Object to access results from.
    concise : bool, optional, default: False
        Whether to print the report in concise mode.

    Returns
    -------
    output : str
        Formatted string of model results.
    """

    # Returns a null report if no results are available
    if np.all(np.isnan(fm.aperiodic_params_)):
        return _no_model_str(concise)
    # Create the formatted strings for printing
    str_lst = [

        # Header
        '=',
        '',
        ' FOOOF - POWER SPECTRUM MODEL',
        '',

        # Frequency range and resolution
        'The model was run on the frequency range {} - {} Hz'.format(
            int(np.floor(fm.freq_range[0])), int(np.ceil(fm.freq_range[1]))),
        'Frequency Resolution is {:1.2f} Hz'.format(fm.freq_res),
        '',

        # Aperiodic parameters
        ('Aperiodic Parameters (offset, ' + ('knee, ' if fm.aperiodic_mode == 'knee' else '') + \
         'exponent): '),
        ', '.join(['{:2.4f}'] * len(fm.aperiodic_params_)).format(*fm.aperiodic_params_),
        '',

        # Peak parameters
        '{} peaks were found:'.format(
            len(fm.peak_params_)),
        *['CF: {:6.2f}, PW: {:6.3f}, BW: {:5.2f}'.format(op[0], op[1], op[2]) \
          for op in fm.peak_params_],
        '',

        # Goodness if fit
        'Goodness of fit metrics:',
        'R^2 of model fit is {:5.4f}'.format(fm.r_squared_),
        'Error of the fit is {:5.4f}'.format(fm.error_),
        '',

        # Goodness if fit
        'BIC results:',
        'Number of models tested: {:1.0f}'.format(fm.bic_data_['n_models_tested']),
        'Best model contains: {:1.0f} peaks'.format(fm.bic_data_['best_index']),
        '',

        # Footer
        '='
    ]
    output = _format(str_lst, concise)

    return output

def save_fm(fm, file_name, file_path=None, append=False,
            save_results=False, save_settings=False, save_data=False):
    """Save out data, results and/or settings from a FOOOF object into a JSON file.

    Parameters
    ----------
    fm : FOOOF
        Object to save data from.
    file_name : str or FileObject
        File to save data to.
    file_path : str, optional
        Path to directory to save to. If None, saves to current directory.
    append : bool, optional, default: False
        Whether to append to an existing file, if available.
        This option is only valid (and only used) if 'file_name' is a str.
    save_results : bool, optional
        Whether to save out FOOOF model fit results.
    save_settings : bool, optional
        Whether to save out FOOOF settings.
    save_data : bool, optional
        Whether to save out input data.

    Raises
    ------
    ValueError
        If the save file is not understood.
    """

    # Convert object to dictionary & convert all arrays to lists, for JSON serializing
    obj_dict = dict_array_to_lst(fm.__dict__)

    # Set and select which variables to keep. Use a set to drop any potential overlap
    #   Note that results also saves frequency information to be able to recreate freq vector
    keep = set((get_description()['results'] + get_description()['meta_data'] if save_results else []) + \
               (get_description()['settings'] if save_settings else []) + \
               (get_decription()['data'] if save_data else []))
    obj_dict = dict_select_keys(obj_dict, keep)

    # Save out - create new file, (creates a JSON file)
    if isinstance(file_name, str) and not append:
        with open(fpath(file_path, fname(file_name, 'json')), 'w') as outfile:
            json.dump(obj_dict, outfile)

    # Save out - append to file_name (appends to a JSONlines file)
    elif isinstance(file_name, str) and append:
        with open(fpath(file_path, fname(file_name, 'json')), 'a') as outfile:
            json.dump(obj_dict, outfile)
            outfile.write('\n')

    # Save out - append to given file object (appends to a JSONlines file)
    elif isinstance(file_name, io.IOBase):
        json.dump(obj_dict, file_name)
        file_name.write('\n')

    else:
        raise ValueError("Save file not understood.")

class FOOOFResults(namedtuple('FOOOFResults', ['aperiodic_params', 'peak_params',
                                               'r_squared', 'error', 'gaussian_params', 'bic_data'])):
    """Model results from parameterizing a power spectrum.

    Parameters
    ----------
    aperiodic_params : 1d array
        Parameters that define the aperiodic fit. As [Offset, (Knee), Exponent].
        The knee parameter is only included if aperiodic is fit with knee.
    peak_params : 2d array
        Fitted parameter values for the peaks. Each row is a peak, as [CF, PW, BW].
    r_squared : float
        R-squared of the fit between the full model fit and the input data.
    error : float
        Error of the full model fit.
    gaussian_params : 2d array
        Parameters that define the gaussian fit(s).
        Each row is a gaussian, as [mean, height, standard deviation].

    Notes
    -----
    This object is a data object, based on a NamedTuple, with immutable data attributes.
    """
    __slots__ = ()

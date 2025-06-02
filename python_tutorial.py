# FOOOF_MS (Python ms-specparam)

# This is a tutorial for the ms-specparam algorithm in Python.
# If you have not already downloaded the fooof-ms folder, which contains all of
# the necessary code for the algorithm to run locally on your device, you should
# do so before continuing.

# Define Utility Functions and Import Dependencies

# In this section, we import all necessary libraries including NumPy, FOOOF
# (specparam), and the new FOOOF_MS (ms_specparam). The FOOOF_MS class extends
# FOOOF to include BIC-based model selection to determine the optimal number of
# peaks in the spectral model. Since FOOOF_MS is a subclass of FOOOF, it retains
# many of the same functionalities as a FOOOF object.

# (Note: For direct interepretability of the code below, we will refer to
# specparam and ms_specparam as FOOOF and FOOOF_MS, respectively, for the
# remainder of this tutorial.)

# Import dependencies
import sys
# Add fooof-ms (ms_specparam) folder to syspath
sys.path.append('your_path_here')
# example sys.path.append('/Users/lucwilson/Desktop/fooof-ms')

import numpy as np
from fooof.sim.gen import gen_power_spectrum
from fooof import FOOOF, FOOOFGroup # Import standard specparam
from fooof_ms import FOOOF_MS, FOOOFGroup_ms  # Ensure fooof_ms (ms_specparam) is installed and imported


# Generate Synthetic 1D Spectrum

# Here, we generate a single synthetic power spectrum using pre-defined
# aperiodic and periodic parameters. The result simulates a realistic neural
# power spectrum which we will use to compare model outputs between FOOOF with
# FOOOF_MS.

# Set seed for reproducibility
np.random.seed(42)

# Generate synthetic spectrum with 3 peaks
freq_range = [1, 40]
freq_1d, power_1d = gen_power_spectrum(
    freq_range,
    aperiodic_params=[1, 1],
    periodic_params=[
        [8, 0.4, 1],
        [14, 0.25, 1.5],
        [21, 0.35, 0.8]],
    nlv=0.05
)
freq_1d = np.array(freq_1d)
power_1d = np.array(power_1d)


## Fit FOOOF to 1D Spectrum

# In this section, we apply FOOOF to a single power spectrum.

print("FOOOF (Fixed 3 Peaks):")
fm_std = FOOOF(max_n_peaks=6)
fm_std.report(freq_1d, power_1d) # Import data to model, fit data, report results


# Fit FOOOF_MS to 1D Spectrum

# In this section, we apply FOOOF_MS to a single power spectrum.
# The result includes an automatically selected model based on BIC evaluation.

print("FOOOF_MS (Fixed 3 Peaks):")
fm_ms = FOOOF_MS(max_n_peaks=6)
fm_ms.report(freq_1d, power_1d) # Import data to model, fit data, report results
fm_ms.report_bic() # Report results from model selection
fm_ms.save('FOOOF_MS_output', 'your_path_here', False, True) # Save results


# Generate Synthetic 2D Spectra

# Here, we generate a group (N=64) of synthetic spectra, each with one peak,
# with varying parameters. These will be used to demonstrate group-level model
# fitting.

# Generate 64 random spectra
np.random.seed(42)
freq_range = [1, 40]
n_spectra = 64

# Preallocate list
power_2d = []
peak_list = []
for a in range(n_spectra):
    # Random aperiodic parameters: offset (0.5–1.5), exponent (0.5–1.5)
    ap = np.random.uniform(0.5, 1.5, size=2)

    # Randomize number of peaks (2–4 per spectrum)
    n_peaks = np.random.randint(1, 2)

    # Each peak: [center (3–35 Hz), height (0.2–0.8), width (0.8–2)]
    peaks = [
        [np.random.uniform(3, 35),  # center
         np.random.uniform(0.2, 0.8),  # height
         np.random.uniform(0.8, 2)]  # width
        for _ in range(n_peaks)
    ]

    # Generate spectrum
    freq_2d, power = gen_power_spectrum(
        freq_range,
        aperiodic_params=ap,
        periodic_params=peaks,
        nlv=0.05  # fixed noise level
    )
    peak_list.append(np.array(peaks))
    power_2d.append(power)

# Convert to NumPy array of shape (64, n_freqs)
power_2d = np.array(power_2d)
freq_2d = np.array(freq_2d)


# Fit FOOOFGroup to Multiple Spectra
# Here, we fit our 64 simulated spectra using the FOOOFGroup class, and display
# the results. We use alternative functions from the 1d version above to load
# data and report results, but either approach is appropriate.

fg_std = FOOOFGroup(max_n_peaks=6)
fg_std.add_data(freq_2d,power_2d) # Import data to model
fg_std.report() #  Fit data, report results (281 peaks)


# Fit FOOOFGroup to Multiple Spectra
# Here, we fit our 64 simulated spectra using the FOOOFGroup_ms class, and
# display the results. Note the differences in number of peaks fit between
# FOOOFGroup and FOOOFGroup_ms.

fg = FOOOFGroup_ms(max_n_peaks=6)
fg.add_data(freq_2d,power_2d) # Import data to model
fg.report() # Fit data, report results (65 peaks)
fg.save('FOOOFGroup_MS_output', 'your_path_here', False, True, True, True) # Save output


# Load Model Results
# Here, we demonstrate how to reload results from the previous FOOOFGroup_ms
# analysis, which has preserved the newly added BIC data.

fg2 = FOOOFGroup_ms()
fg2.load('FOOOFGroup_MS_output.json','your_path_here') # Load results
fg2.print_results() # Print report
fg2.plot() # Generate group plots for aperiodic parameters and peak locations

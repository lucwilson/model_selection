%% ms-specparam
% This is a tutorial for the ms-specparam algorithm.
% If you have not already downloaded the ms-specparam folder, which contains all of
% the necessary code for the algorithm to run locally on your device, you should
% do so before continuing.

%% Import dependencies
% Add ms-specparam folder to directory
addpath('your_path_here')
% example addpath('/Users/lucwilson/Desktop/ms-specparam')

%% Generate frequencies, 1D spectrum
rng(42);
Freqs = 1:0.5:40;
TF = zeros(1,1,79);
aperiodic_params = [1 1];
periodic_params = [8 0.4 1; 14 0.25 1.5; 21 0.35 0.8];
noise = randn(1,1,79);
TF(1,1,:) = build_model(Freqs, aperiodic_params, 'fixed', periodic_params, @gaussian);
TF = 10.^(TF+noise*0.05); % function takes PSD measured in units of power

%% Generate options structure from pre-sets
help build_opt
opt = build_opt();

%% Fit ms-specparam to 1D spectrum
help ms_specparam
[fs, fg] = ms_specparam(TF,Freqs, opt, opt.hOT);

%% Display results
figure('Position',[100 100 900 400])

subplot(1,2,1), hold on % ms-specparam with lowest-BIC model 
plot(Freqs,log10(squeeze(TF)),'k','LineWidth',2) % Plot PSD
plot(Freqs, gen_aperiodic(Freqs, fg.aperiodic_params, 'fixed'),'--b','LineWidth',2) % Aperiodic fit
plot(Freqs,log10(fg.fooofed_spectrum),'r','LineWidth',2) % Full model fit
text(30, 0.8,['R^2 = ' num2str(round(fg.r_squared,4))]) % Add R2 
text(30, 0.74,['BIC = ' num2str(round(fg.BIC,2))]) % Add BIC
text(30, 0.68,['Peaks = ' num2str(size(fg.peak_params,1))]) % Add number of peaks
xlabel('Frequency (Hz)'), ylabel('log_{10}(Power)'), title('Lowest-BIC model (*ms-specparam)')
% This represents the output of ms-specparam.

subplot(1,2,2), hold on % ms-specparam with highest R^2 model (similar to specparam)
plot(Freqs,log10(squeeze(TF)),'k','LineWidth',2)
plot(Freqs, gen_aperiodic(Freqs, fg.models(end).aperiodic_params, 'fixed'),'--b','LineWidth',2)
plot(Freqs,build_model(Freqs, fg.models(end).aperiodic_params, 'fixed', fg.models(end).peak_params, @gaussian),'r','LineWidth',2)
text(30, 0.8,['R^2 = ' num2str(round(fg.models(end).r_squared,4))])
text(30,0.74,['BIC = ' num2str(round(fg.models(end).BIC,2))])
text(30, 0.68,['Peaks = ' num2str(size(fg.models(end).peak_params,1))])
xlabel('Frequency (Hz)'), ylabel('log_{10}(Power)'), title('Highest R^2 model (*similar to specparam)')
% N.B. This output will differ slightly from default specparam due to added 
% whole-model optimization at the end (absent in specparam). However, it
% preserves the same number of peaks as specparam.

%% Support functions taken from ms_specparam.m
function model_fit = build_model(freqs, ap_pars, ap_type, pk_pars, peak_function)
%     Builds a full spectral model from parameters.
%
%     Parameters
%     ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       ap_pars : 1xm array
%           Parameter estimates for aperiodic fit.
%       pk_pars : kx3 array, where k = No. of peaks.
%           Guess parameters for peak fits.
%       pk_type : {'gaussian', 'cauchy', 'best'}
%           Which types of peaks are being fitted.
%
%       Returns
%       -------
%       model_fit : 1xn array
%           Model power spectrum, in log10-space

    ap_fit = gen_aperiodic(freqs, ap_pars, ap_type);
    model_fit = ap_fit;
    if length(pk_pars) > 1
        for peak = 1:size(pk_pars,1)
            model_fit = model_fit + peak_function(freqs,pk_pars(peak,1),...
                pk_pars(peak,2),pk_pars(peak,3));
        end  
    end
end

function ap_vals = gen_aperiodic(freqs,aperiodic_params,aperiodic_mode)
%       Generate aperiodic values, from parameter definition.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%       	Frequency vector to create aperiodic component for.
%       aperiodic_params : 1x3 array
%           Parameters that define the aperiodic component.
%       aperiodic_mode : {'fixed', 'knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       ap_vals : 1d array
%           Generated aperiodic values.

    switch aperiodic_mode
        case 'fixed'  % no knee
            ap_vals = expo_nk_function(freqs,aperiodic_params);
        case 'knee'
            ap_vals = expo_function(freqs,aperiodic_params);
        case 'floor'
            ap_vals = expo_fl_function(freqs,aperiodic_params);
    end
end

function ys = gaussian(freqs, mu, hgt, sigma)
%       Gaussian function to use for fitting.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create gaussian fit for.
%       mu, hgt, sigma : doubles
%           Parameters that define gaussian function (centre frequency,
%           height, and standard deviation).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for gaussian function.

    ys = hgt*exp(-(((freqs-mu)./sigma).^2) /2);

end

function ys = cauchy(freqs, ctr, hgt, gam)
%       Cauchy function to use for fitting.
% 
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create cauchy fit for.
%       ctr, hgt, gam : doubles
%           Parameters that define cauchy function (centre frequency,
%           height, and "standard deviation" [gamma]).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for cauchy function.

    ys = hgt./(1+((freqs-ctr)/gam).^2);

end

function ys = expo_function(freqs,params)
%       Exponential function to use for fitting 1/f, with a 'knee' (maximum at low frequencies).
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x3 array (offset, knee, exp)
%           Parameters (offset, knee, exp) that define Lorentzian function:
%           y = 10^offset * (1/(knee + x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential function.

    ys = params(1) - log10(abs(params(2)) +freqs.^params(3));

end

function ys = expo_nk_function(freqs, params)
%       Exponential function to use for fitting 1/f, without a 'knee'.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x2 array (offset, exp)
%           Parameters (offset, exp) that define Lorentzian function:
%           y = 10^offset * (1/(x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential (no-knee) function.

    ys = params(1) - log10(freqs.^params(2));

end

function ys = expo_fl_function(f, params)

    ys = log10(f.^(params(1)) * 10^(params(2)) + params(3));

end

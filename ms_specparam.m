%% Initialize simulation/analysis parameters
opt.freq_range          = [1 40];
opt.peak_width_limits   = [1 12];
opt.max_peaks           = 6;
opt.min_peak_height     = 1 / 10; % convert from dB to B
opt.aperiodic_mode      = 'fixed'; % 'knee'
opt.peak_threshold      = 2;   % 2 std dev: parameter for interface simplification
opt.border_threshold    = 1;   % 1 std dev: proximity to edge of spectrum, static in Python 
opt.return_spectrum     = 0;   % SPM/FT: set to 1
opt.power_line          = '-5'; % otherwise 50, 60 Hz
opt.proximity_threshold = 0.75;
opt.optim_obj           = 'negloglike'; % negloglike
opt.peak_type           = 'gaussian'; % 'cauchy', for interface simplification
opt.guess_weight        = 'none'; % otherwise 'weak' or 'strong'
opt.thresh_after        = true;   % Only when guess weight > 'None'
opt.hOT                 = 1; % 0 if no optimization toolbox
% generate spectra (same spectra, different noise levels)
rng(123);
Freqs = 0.5:0.5:100;
TF = zeros(4,1,200);
ap = [-4.3 1.2];
pp = [4 0.7 1.3; 11 0.9 2; 17 0.6 5];
noise = randn(1,1,200);
TF(1,1,:) = build_model(Freqs, ap, 'fixed', pp, @gaussian);
for k = 1:3
    TF(k+1,1,:) = TF(1,1,:)+noise*0.05.*k;
end
TF = 10.^TF; % function takes PSD measured in units of power

[fs, fg] = ms_specparam(TF,Freqs, opt, 1);


%% Main algorithm
function [fs, fg] = ms_specparam(TF, Freqs, opt, hOT)
    % Find all frequency values within user limits
    fMask = (round(Freqs.*10)./10 >= opt.freq_range(1)) & (round(Freqs.*10)./10 <= opt.freq_range(2)) & ~mod(sum(abs(round(Freqs.*10)./10-[1;2;3].*str2double(opt.power_line)) >= 2),3);
    fs = Freqs(fMask);
    spec = log10(squeeze(TF(:,1,fMask))); % extract log spectra
    nChan = size(TF,1);
    if nChan == 1, spec = spec'; end
    % Initalize FOOOF structs
    fg(nChan) = struct(...
            'aperiodic_params', [],...
            'peak_params',      [],...
            'peak_types',       '',...
            'ap_fit',           [],...
            'fooofed_spectrum', [],...
            'peak_fit',         [],...
            'error',            [],...
            'r_squared',        []);
    % Iterate across channels
    for chan = 1:nChan
        % Fit aperiodic
        aperiodic_pars = robust_ap_fit(fs, spec(chan,:), opt.aperiodic_mode);
        % Remove aperiodic
        flat_spec = flatten_spectrum(fs, spec(chan,:), aperiodic_pars, opt.aperiodic_mode);
        % estimate valid peaks (and determine max n)
        [est_pars, peak_function] = est_peaks(fs, flat_spec, opt.max_peaks, opt.peak_threshold, opt.min_peak_height, ...
            opt.peak_width_limits/2, opt.proximity_threshold, opt.border_threshold, opt.peak_type);
        model = struct();
        for pk = 0:size(est_pars,1)
            peak_pars = est_fit(est_pars(1:pk,:), fs, flat_spec, opt.peak_width_limits/2, opt.peak_type, opt.guess_weight,hOT);
            % Refit aperiodic
            aperiodic = spec(chan,:);
            for peak = 1:size(peak_pars,1)
                aperiodic = aperiodic - peak_function(fs,peak_pars(peak,1), peak_pars(peak,2), peak_pars(peak,3));
            end
            aperiodic_pars = simple_ap_fit(fs, aperiodic, opt.aperiodic_mode);
            guess = peak_pars;
            if ~isempty(guess)
                lb = [max([ones(size(guess(1:pk,:),1),1).*fs(1) guess(1:pk,1)-guess(1:pk,3)*2],[],2),zeros(size(guess(1:pk,2))),ones(size(guess(1:pk,3)))*opt.peak_width_limits(1)/2]';
                ub = [min([ones(size(guess(1:pk,:),1),1).*fs(end) guess(1:pk,1)+guess(1:pk,3)*2],[],2),inf(size(guess(1:pk,2))),ones(size(guess(1:pk,3)))*opt.peak_width_limits(2)/2]';

            else
                lb = [];
                ub = [];
            end
            switch opt.aperiodic_mode
                case 'fixed'
                    lb = [-inf; 0; lb(:)];
                    ub = [inf; inf; ub(:)];
                case 'knee'
                    lb = [-inf; 0; 0; lb(:)];
                    ub = [inf; 100; inf; ub(:)];
            end
            if opt.return_spectrum
                fg(chan).power_spectrum = spec(chan,:);
            end
            guess = guess(1:pk,:)';
            guess = [aperiodic_pars'; guess(:)];
            options = optimset('Display', 'off', 'TolX', 1e-7, 'TolFun', 1e-9, ...
                'MaxFunEvals', 5000, 'MaxIter', 5000); % Tuned options 
            try
                params = fmincon(@err_fm_constr,guess,[],[],[],[], ...
                    lb,ub,[],options,fs,spec(chan,:),opt.aperiodic_mode,opt.peak_type);
            catch
                a = 0; % for catching errors
            end
            switch opt.aperiodic_mode
                case 'fixed'
                    aperiodic_pars_tmp = params(1:2);
                    if length(params) > 3
                        peak_pars_tmp = reshape(params(3:end),[3 length(params(3:end))./3])';
                    end
                case 'knee'
                    aperiodic_pars_tmp = params(1:3);
                    if length(params) > 3
                        peak_pars_tmp = reshape(params(4:end),[3 length(params(4:end))./3])';
                    end
            end
            % Generate model fit
            ap_fit = gen_aperiodic(fs, aperiodic_pars_tmp, opt.aperiodic_mode);
            model_fit = ap_fit;
            if length(params) > 3
                for peak = 1:size(peak_pars_tmp,1)
                    model_fit = model_fit + peak_function(fs,peak_pars_tmp(peak,1),...
                        peak_pars_tmp(peak,2),peak_pars_tmp(peak,3));
                end  
            else
                peak_pars_tmp = [0 0 0];
            end
            % Calculate model error
            MSE = sum((spec(chan,:) - model_fit).^2)/length(model_fit);
            rsq_tmp = corrcoef(spec(chan,:),model_fit).^2;
            loglik = -length(model_fit)/2.*(1+log(MSE)+log(2*pi));
            AIC = 2.*(length(params)-loglik);
            BIC = length(params).*log(length(model_fit))-2.*loglik;
            model(pk+1).aperiodic_params = aperiodic_pars_tmp;
            model(pk+1).peak_params = peak_pars_tmp;
            model(pk+1).MSE = MSE;
            model(pk+1).r_squared = rsq_tmp(2);
            model(pk+1).loglik = loglik;
            model(pk+1).AIC = AIC;
            model(pk+1).BIC = BIC;
            model(pk+1).BF = exp((BIC-model(1).BIC)./2);
        end
        % insert data from best model
        
        [~,mi] = min([model.BIC]);
        
        aperiodic_pars = model(mi).aperiodic_params;
        peak_pars = model(mi).peak_params;
        % Return FOOOF results
        aperiodic_pars(2) = abs(aperiodic_pars(2));
        fg(chan).aperiodic_params   = aperiodic_pars;
        fg(chan).peak_params        = peak_pars;
        fg(chan).peak_types         = func2str(peak_function);
        fg(chan).ap_fit             = 10.^gen_aperiodic(fs, aperiodic_pars, opt.aperiodic_mode);
        fg(chan).fooofed_spectrum   = 10.^build_model(fs, aperiodic_pars, opt.aperiodic_mode, peak_pars, peak_function);
        fg(chan).peak_fit           = fg(chan).fooofed_spectrum ./ fg(chan).ap_fit;
        fg(chan).error              = model(mi).MSE;
        fg(chan).r_squared          = model(mi).r_squared;
        fg(chan).loglik             = model(mi).loglik; % log-likelihood
        fg(chan).AIC                = model(mi).AIC;
        fg(chan).BIC                = model(mi).BIC;
        fg(chan).models             = model;
        %plot(fs', [fg(chan).ap_fit', fg(chan).peak_fit', fg(chan).fooofed_spectrum'])
    end
end

%% ===== GENERATE APERIODIC =====
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


%% ===== CORE MODELS =====
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


%% ===== FITTING ALGORITHM =====
function aperiodic_params = simple_ap_fit(freqs, power_spectrum, aperiodic_mode)
%       Fit the aperiodic component of the power spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

%       Set guess params for lorentzian aperiodic fit, guess params set at init
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
        'MaxFunEvals', 10000, 'MaxIter', 10000);

    switch (aperiodic_mode)
        case 'fixed'  % no knee
            exp_guess = -(power_spectrum(end)-power_spectrum(1))./log10(freqs(end)./freqs(1));
            guess_vec = [power_spectrum(1), exp_guess];
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs, power_spectrum);
        case 'knee'
            exp_guess = -(power_spectrum(end)-power_spectrum(1))./log10(freqs(end)./freqs(1));
            guess_vec = [power_spectrum(1),0, exp_guess];
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs, power_spectrum);
    end

end

function aperiodic_params = robust_ap_fit(freqs, power_spectrum, aperiodic_mode)
%       Fit the aperiodic component of the power spectrum robustly, ignoring outliers.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

    % Do a quick, initial aperiodic fit
    popt = simple_ap_fit(freqs, power_spectrum, aperiodic_mode);
    initial_fit = gen_aperiodic(freqs, popt, aperiodic_mode);

    % Flatten power_spectrum based on initial aperiodic fit
    flatspec = power_spectrum - initial_fit;

    % Flatten outliers - any points that drop below 0
    flatspec(flatspec(:) < 0) = 0;

    % Use percential threshold, in terms of # of points, to extract and re-fit
    perc_thresh = prctile(flatspec, 0.025);
    perc_mask = flatspec <= perc_thresh;
    freqs_ignore = freqs(perc_mask);
    spectrum_ignore = power_spectrum(perc_mask);

    % Second aperiodic fit - using results of first fit as guess parameters

    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
        'MaxFunEvals', 10000, 'MaxIter', 10000);
    guess_vec = popt;

    switch (aperiodic_mode)
        case 'fixed'  % no knee
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs_ignore, spectrum_ignore);
        case 'knee'
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs_ignore, spectrum_ignore);
    end
end

function spectrum_flat = flatten_spectrum(freqs, power_spectrum, robust_aperiodic_params, aperiodic_mode)
%       Flatten the power spectrum by removing the aperiodic component.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       robust_aperiodic_params : 1x2 or 1x3 array (see aperiodic_mode)
%           Parameter estimates for aperiodic fit.
%       aperiodic_mode : 1 or 2
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       spectrum_flat : 1xn array
%           Flattened (aperiodic removed) power spectrum.


spectrum_flat = power_spectrum - gen_aperiodic(freqs,robust_aperiodic_params,aperiodic_mode);

end

function [guess_params,peak_function] = est_peaks(freqs, flat_iter, max_n_peaks, peak_threshold, min_peak_height, gauss_std_limits, proxThresh, bordThresh, peakType)
%       Iteratively fit peaks to flattened spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       max_n_peaks : double
%           Maximum number of gaussians to fit within the spectrum.
%       peak_threshold : double
%           Threshold (in standard deviations of noise floor) to detect a peak.
%       min_peak_height : double
%           Minimum height of a peak (in log10).
%       gauss_std_limits : 1x2 double
%           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
%       proxThresh : double
%           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
%       peakType : {'gaussian', 'cauchy', 'both'}
%           Which types of peaks are being fitted
%       guess_weight : {'none', 'weak', 'strong'}
%           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       guess_params : mx3 array, where m = No. of peaks.
%           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].

    switch peakType 
        case 'gaussian' % gaussian only
            peak_function = @gaussian; % Identify peaks as gaussian
            % Initialize matrix of guess parameters for gaussian fitting.
            guess_params = zeros(max_n_peaks, 3);
            % Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian.
            % Stopping procedure based on either the limit on # of peaks,
            % or the relative or absolute height thresholds.
            for guess = 1:max_n_peaks
                % Find candidate peak - the maximum point of the flattened spectrum.
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);

                % Stop searching for peaks once max_height drops below height threshold.
                if max_height <= peak_threshold * std(flat_iter)
                    break
                end

                % Set the guess parameters for gaussian fitting - mean and height.
                guess_freq = freqs(max_ind);
                guess_height = max_height;

                % Halt fitting process if candidate peak drops below minimum height.
                if guess_height <= min_peak_height
                    break
                end

                % Data-driven first guess at standard deviation
                % Find half height index on each side of the center frequency.
                half_height = 0.5 * max_height;

                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height)+1;

                % Keep bandwidth estimation from the shortest side.
                % We grab shortest to avoid estimating very large std from overalapping peaks.
                % Grab the shortest side, ignoring a side if the half max was not found.
                % Note: will fail if both le & ri ind's end up as None (probably shouldn't happen).
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate std from FWHM. Calculate FWHM, converting to Hz, get guess std from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_std = fwhm / (2 * sqrt(2 * log(2)));

                % Check that guess std isn't outside preset std limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_std < gauss_std_limits(1)
                    guess_std = gauss_std_limits(1);
                end
                if guess_std > gauss_std_limits(2)
                    guess_std = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq, guess_height, guess_std];

                % Subtract best-guess gaussian.
                peak_gauss = gaussian(freqs, guess_freq, guess_height, guess_std);
                flat_iter = flat_iter - peak_gauss;

            end
            % Remove unused guesses
            guess_params(guess_params(:,1) == 0,:) = [];

            % Check peaks based on edges, and on overlap
            % Drop any that violate requirements.
            guess_params = drop_peak_cf(guess_params, bordThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);
            
        case 'cauchy' % cauchy only
            peak_function = @cauchy; % Identify peaks as cauchy
            guess_params = zeros(max_n_peaks, 3);
            flat_spec = flat_iter;
            for guess = 1:max_n_peaks
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);
                if max_height <= peak_threshold * std(flat_iter)
                    break
                end
                guess_freq = freqs(max_ind);
                guess_height = max_height;
                if guess_height <= min_peak_height
                    break
                end
                half_height = 0.5 * max_height;
                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height);
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate gamma from FWHM. Calculate FWHM, converting to Hz, get guess gamma from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_gamma = fwhm/2;
                % Check that guess gamma isn't outside preset limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_gamma < gauss_std_limits(1)
                    guess_gamma = gauss_std_limits(1);
                end
                if guess_gamma > gauss_std_limits(2)
                    guess_gamma = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq(1), guess_height, guess_gamma];

                % Subtract best-guess cauchy.
                peak_cauchy = cauchy(freqs, guess_freq(1), guess_height, guess_gamma);
                flat_iter = flat_iter - peak_cauchy;

            end
            guess_params(guess_params(:,1) == 0,:) = [];
            guess_params = drop_peak_cf(guess_params, proxThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);
            
    end
end

function model_params = est_fit(guess_params, freqs, flat_spec, gauss_std_limits, peakType, guess_weight,hOT)
%       Iteratively fit peaks to flattened spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       max_n_peaks : double
%           Maximum number of gaussians to fit within the spectrum.
%       peak_threshold : double
%           Threshold (in standard deviations of noise floor) to detect a peak.
%       min_peak_height : double
%           Minimum height of a peak (in log10).
%       gauss_std_limits : 1x2 double
%           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
%       proxThresh : double
%           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
%       peakType : {'gaussian', 'cauchy', 'both'}
%           Which types of peaks are being fitted
%       guess_weight : {'none', 'weak', 'strong'}
%           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       guess_params : mx3 array, where m = No. of peaks.
%           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].

    switch peakType 
        case 'gaussian' % gaussian only

            % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 1, guess_weight, gauss_std_limits,hOT);
            else
                model_params = [];
            end
            
        case 'cauchy' % cauchy only

           % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 2, guess_weight, gauss_std_limits,hOT);
            else
                model_params = [];
            end
            
    end
end

    

function guess = drop_peak_cf(guess, bw_std_edge, freq_range)
%       Check whether to drop peaks based on center's proximity to the edge of the spectrum.
%
%       Parameters
%       ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%
%       Returns
%       -------
%       guess : qx3 where q <= m No. of peaks.
%           Guess parameters for peak fits.

    cf_params = guess(:,1)';
    bw_params = guess(:,3)' * bw_std_edge;

    % Check if peaks within drop threshold from the edge of the frequency range.

    keep_peak = abs(cf_params-freq_range(1)) > bw_params & ...
        abs(cf_params-freq_range(2)) > bw_params;

    % Drop peaks that fail the center frequency edge criterion
    guess = guess(keep_peak,:);

end

function guess = drop_peak_overlap(guess, proxThresh)
%       Checks whether to drop gaussians based on amount of overlap.
%
%       Parameters
%       ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%       proxThresh: double
%           Proximity threshold (in st. dev. or gamma) between two peaks.
%
%       Returns
%       -------
%       guess : qx3 where q <= m No. of peaks.
%           Guess parameters for peak fits.
%
%       Note
%       -----
%       For any gaussians with an overlap that crosses the threshold,
%       the lowest height guess guassian is dropped.

    % Sort the peak guesses, so can check overlap of adjacent peaks
    guess = sortrows(guess);

    % Calculate standard deviation bounds for checking amount of overlap

    bounds = [guess(:,1) - guess(:,3) * proxThresh, ...
        guess(:,1), guess(:,1) + guess(:,3) * proxThresh];

    % Loop through peak bounds, comparing current bound to that of next peak
    drop_inds =  [];

    for ind = 1:size(bounds,1)-1

        b_0 = bounds(ind,:);
        b_1 = bounds(ind + 1,:);

        % Check if bound of current peak extends into next peak
        if b_0(2) > b_1(1)
            % If so, get the index of the gaussian with the lowest height (to drop)
            drop_inds = [drop_inds (ind - 1 + find(guess(ind:ind+1,2) == ...
                min(guess(ind,2),guess(ind+1,2))))];
        end
    end
    % Drop any peaks guesses that overlap too much, based on threshold.
    guess(drop_inds,:) = [];
    
    % Readjust order by amplitude
    guess = sortrows(guess,2,'descend');
end

function peak_params = fit_peak_guess(guess, freqs, flat_spec, peak_type, guess_weight, std_limits, hOT)
%     Fits a group of peak guesses with a fit function.
%
%     Parameters
%     ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       peakType : {'gaussian', 'cauchy', 'best'}
%           Which types of peaks are being fitted.
%       guess_weight : 'none', 'weak', 'strong'
%           Parameter to weigh initial estimates during optimization.
%       std_limits: 1x2 array
%           Minimum and maximum standard deviations for distribution.
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       peak_params : mx3, where m =  No. of peaks.
%           Peak parameters post-optimization.

    
    if hOT % Use OptimToolbox for fmincon
        lb = [max([ones(size(guess,1),1).*freqs(1) guess(:,1)-guess(:,3)*2],[],2),zeros(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(1)];
        ub = [min([ones(size(guess,1),1).*freqs(end) guess(:,1)+guess(:,3)*2],[],2),inf(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(2)];
        options = optimset('Display', 'off', 'TolX', 1e-3, 'TolFun', 1e-5, ...
            'MaxFunEvals', 5000, 'MaxIter', 5000); % Tuned options       
        peak_params = fmincon(@error_model_constr,guess,[],[],[],[], ...
            lb,ub,[],options,freqs,flat_spec, peak_type);
    else % Use basic simplex approach, fminsearch, with guess_weight
        options = optimset('Display', 'off', 'TolX', 1e-5, 'TolFun', 1e-7, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);
        peak_params = fminsearch(@error_model,...
            guess, options, freqs, flat_spec, peak_type, guess, guess_weight);
    end
end

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

%% ===== ERROR FUNCTIONS =====
function err = error_expo_nk_function(params,xs,ys)
    ym = -log10(xs.^params(2)) + params(1);
    err = sum((ys - ym).^2);
end

function err = error_expo_function(params,xs,ys)
    ym = expo_function(xs,params);
    err = sum((ys - ym).^2);
end

function err = error_model(params, xVals, yVals, peak_type, guess, guess_weight)
    fitted_vals = 0;
    weak = 1E2;
    strong = 1E7;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    switch guess_weight
        case 'none'
            err = sum((yVals - fitted_vals).^2);
        case 'weak' % Add small weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 weak*sum((params(:,1)-guess(:,1)).^2) + ...
                 weak*sum((params(:,2)-guess(:,2)).^2);
        case 'strong' % Add large weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 strong*sum((params(:,1)-guess(:,1)).^2) + ...
                 strong*sum((params(:,2)-guess(:,2)).^2);
    end
end

function err = error_model_constr(params, xVals, yVals, peak_type)
    fitted_vals = 0;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    err = sum((yVals - fitted_vals).^2);
end

function err = err_fm_constr(params, xVals, yVals, aperiodic_mode, peak_type)
    switch (aperiodic_mode)
        case 'fixed'  % no knee
        npk = (length(params)-2)/3;
        fitted_vals = -log10(xVals.^params(2)) + params(1);
        case 'knee'
        npk = (length(params)-3)/3;
        fitted_vals = params(1) - log10(abs(params(2)) +xVals.^params(3));
    end
    for set = 1:npk
        switch peak_type 
            case 'gaussian' % gaussian only
                fitted_vals = fitted_vals + gaussian(xVals, params(3.*set), params(3*set+1), params(3*set+2));
            case 'cauchy' % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(3.*set), params(3*set+1), params(3*set+2));
        end
    end
    err = sum((yVals - fitted_vals).^2);
end

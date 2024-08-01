%% generate random parameters for fixed aperiodic PSDs
psd_os = struct('ap_pars',[],'peak_pars',[]);

for s = 0:4
    psd_os(s*1000+1:(s+1)*1000) = gen_rand_spec_par([0.5 2;-8.1 -1.5],[3 35],[0.1 1.5],[1 3],2,s,1000);
end

%% Generate spectra from parameters
specs = zeros(5000,200);
fs = 0.5:0.5:100;
for s = 1:5000
    specs(s,:) = gen_spectrum(fs, psd_os(s).ap_pars([end 1]), psd_os(s).peak_pars, 0, 'gaussian', 1);
end

%% Extract stats and peaks from ms_specparam
fgtmp = ms_specparam.Options.FOOOF.data;
rstmp = struct();
npeak_simulated_BIC = zeros(1000,5);
npeak_recovered_BIC = zeros(1000,5);
npeak_total_BIC = zeros(1000,5);

for sim = 1:length(fgtmp)
    
    [~, tmp] = min([fgtmp(sim).models.BIC]);
    peaks = fgtmp(sim).models(tmp).peak_params;
    rstmp(sim).aper_pars = fgtmp(sim).models(tmp).aperiodic_params;
    if sim > 1000
        rstmp(sim).exp_peaks = nan(size(psd_os(sim).peak_pars));
    end
    if ~any(peaks(1,:))
        continue
    end
    zn = [psd_os(sim).peak_pars(:,1)-2.*psd_os(sim).peak_pars(:,3) psd_os(sim).peak_pars(:,1)+2.*psd_os(sim).peak_pars(:,3)];
    if size(psd_os(sim).peak_pars,1)
        [tmp, order] = sort(psd_os(sim).peak_pars(:,2),'descend');
        pktmp = peaks;
        for ind = 1:length(order)
            pk = order(ind);
            expi = find((zn(pk,1) <= pktmp(:,1)) & (zn(pk,2) >= pktmp(:,1))); % filter by cf range only
            exp_peaks = pktmp(expi,:);
            if size(exp_peaks,1) > 1
                [tmp, rmv] = max(exp_peaks(:,2));
                exp_peaks = exp_peaks(rmv,:);
                pktmp(expi(rmv),:) = [];
            else
                pktmp(expi,:) = [];
            end
            if isempty(exp_peaks)
                exp_peaks = nan(1,3);
            end
            rstmp(sim).exp_peaks(pk,:) = exp_peaks;
        end
    end
    
    k = mod(sim,1000);
    if ~mod(sim,1000)
       k = 1000; 
    end
    npeak_simulated_BIC(k,ceil(sim./1000)) = ceil(sim./1000)-1;
    if sim > 1000
        npeak_recovered_BIC(k,ceil(sim./1000)) = length(rstmp(sim).exp_peaks(~isnan(rstmp(sim).exp_peaks(:,1)),1));
    end
    npeak_total_BIC(k,ceil(sim./1000)) = length(peaks(peaks(:,1)~=0,1));
    
end

err_ap_BIC = nan(5000,length(rstmp(1).aper_pars));
for sim = 1:length(psd_os)
    err_ap_BIC(sim,:) = abs(psd_os(sim).ap_pars(1:1:end)-rstmp(sim).aper_pars(end:-1:1)');
end
err_pk_BIC = nan(25000,3);
ctr = 1;
for sim = 1001:length(psd_os)
    err_tmp = abs(psd_os(sim).peak_pars-rstmp(sim).exp_peaks);
    err_pk_BIC(ctr:ctr+size(psd_os(sim).peak_pars,1)-1,:) = err_tmp;
    ctr = ctr+size(psd_os(sim).peak_pars,1);
end

err_pk_BIC(isnan(err_pk_BIC(:,1)),:) = [];

%% Extract stats and peaks from default specparam

fgtmp = default_specparam.Options.FOOOF.data;
rstmp2 = struct();
npeak_simulated = zeros(1000,5);
npeak_recovered = zeros(1000,5);
npeak_total = zeros(1000,5);

for sim = 1:length(fgtmp)
    peaks = fgtmp(sim).peak_params;
    rstmp2(sim).aper_pars = fgtmp(sim).aperiodic_params;
    if sim > 1000
        rstmp2(sim).exp_peaks = nan(size(psd_os(sim).peak_pars));
    end
    if ~any(peaks(1,:))
        continue
    end
    zn = [psd_os(sim).peak_pars(:,1)-2.*psd_os(sim).peak_pars(:,3) psd_os(sim).peak_pars(:,1)+2.*psd_os(sim).peak_pars(:,3)];
    if size(psd_os(sim).peak_pars,1)
        [tmp, order] = sort(psd_os(sim).peak_pars(:,2),'descend');
        pktmp = peaks;
        for ind = 1:length(order)
            pk = order(ind);
            expi = find((zn(pk,1) <= pktmp(:,1)) & (zn(pk,2) >= pktmp(:,1))); % filter by cf range only
            exp_peaks = pktmp(expi,:);
            if size(exp_peaks,1) > 1
                [tmp, rmv] = max(exp_peaks(:,2));
                exp_peaks = exp_peaks(rmv,:);
                pktmp(expi(rmv),:) = [];
            else
                pktmp(expi,:) = [];
            end
            if isempty(exp_peaks)
                exp_peaks = nan(1,3);
            end
            rstmp2(sim).exp_peaks(pk,:) = exp_peaks;
        end
    end
    
    k = mod(sim,1000);
    if ~mod(sim,1000)
       k = 1000; 
    end
    npeak_simulated(k,ceil(sim./1000)) = ceil(sim./1000)-1;
    if sim > 1000
        npeak_recovered(k,ceil(sim./1000)) = length(rstmp2(sim).exp_peaks(~isnan(rstmp2(sim).exp_peaks(:,1)),1));
    end
    npeak_total(k,ceil(sim./1000)) = length(peaks(peaks(:,1)~=0,1));
end

err_ap = nan(5000,length(rstmp2(1).aper_pars));
for sim = 1:length(psd_os)
    err_ap(sim,:) = abs(psd_os(sim).ap_pars(1:1:end)'-rstmp2(sim).aper_pars(end:-1:1));
end
err_pk = nan(25000,3);
ctr = 1;
for sim = 1001:length(psd_os)
    err_tmp = abs(psd_os(sim).peak_pars-rstmp2(sim).exp_peaks);
    err_pk(ctr:ctr+size(psd_os(sim).peak_pars,1)-1,:) = err_tmp;
    ctr = ctr+size(psd_os(sim).peak_pars,1);
end

err_pk(isnan(err_pk(:,1)),:) = [];

%% calculate sensitivity, PPV
sens = sum(npeak_recovered(:,2:5))./sum(npeak_simulated(:,2:5));
sens_bic = sum(npeak_recovered_BIC(:,2:5))./sum(npeak_simulated_BIC(:,2:5));

ppv = sum(npeak_recovered(:,2:5))./sum(npeak_total(:,2:5));
ppv_bic = sum(npeak_recovered_BIC(:,2:5))./sum(npeak_total_BIC(:,2:5));

%% Bootstrap CIs for sensitivity, PPV

sens_bs = zeros(1000,4);
sens_bic_bs = zeros(1000,4);
ppv_bs = zeros(1000,4);
ppv_bic_bs = zeros(1000,4);

for k = 1:1000
   idx = randi(1000,1000,1);
   sens_bs(k,:) = sum(npeak_recovered(idx,2:5))./sum(npeak_simulated(idx,2:5));
   sens_bic_bs(k,:) = sum(npeak_recovered_BIC(idx,2:5))./sum(npeak_simulated_BIC(idx,2:5));
   ppv_bs(k,:) = sum(npeak_recovered(idx,2:5))./sum(npeak_total(idx,2:5));
   ppv_bic_bs(k,:) = sum(npeak_recovered_BIC(idx,2:5))./sum(npeak_total_BIC(idx,2:5));
end

%% Peak detection rates by number of simulated peaks

npeak_mat = zeros(5,7);
npeak_bic_mat = zeros(5,7);
for s = 1:5
    for i = 1:1000
        npeak_mat(s,npeak_total(i,s)+1) =...
            npeak_mat(s,npeak_total(i,s)+1)+1;
        npeak_bic_mat(s,npeak_total_BIC(i,s)+1) =...
            npeak_bic_mat(s,npeak_total_BIC(i,s)+1)+1;
    end
end

vals = zeros(5000,3);
a = ms_specparam;
for s = 1:5000
    vals(s,1) = ceil(s/1000);
    [~, vals(s,2)] = min([a.Options.FOOOF.data(s).models.BIC]);
    vals(s,3) = length([a.Options.FOOOF.data(s).models.BIC]);
end

npeak_mat = zeros(5,7);
for s = 1:5000
    npeak_mat(vals(s,1),vals(s,2)) = npeak_mat(vals(s,1),vals(s,2))+1;
end

figure, hold on
swatch = [86 34 98;91 61 28;162 55 27;198 97 43;236 189 78;108 156 167]./255;
xs = [0 0 0 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4];
ys = [0 1 2 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6];
sz = 1200*[npeak_mat(1,:)./sum(npeak_mat(1,:)) npeak_mat(2,:)./sum(npeak_mat(2,:)) npeak_mat(3,:)./sum(npeak_mat(3,:)) npeak_mat(4,:)./sum(npeak_mat(4,:)) npeak_mat(5,:)./sum(npeak_mat(5,:))];
xs = xs(logical(sz));
ys = ys(logical(sz));
sz = sz(logical(sz));
for p = 1:length(xs)
    if ys(p) == 0
        scatter(xs(p),ys(p),sz(p),[52 72 148]./255,'Filled')
    else
        scatter(xs(p),ys(p),sz(p),[52 72 148]./255,'Filled')
    end
end
yticks(0:6)
xticks(0:4)
xlim([-0.35 4.5])
ylim([-0.35 6.1])
xlabel('# Simulated peaks','Fontsize',14)
ylabel('# Fit peaks','Fontsize',14)

%% Prepare error distributions for violin plots
err_ap_exp_van = log10(abs(err_ap(:,1)));
err_ap_off_van = log10(abs(err_ap(:,end)));
err_ap_cf_van = log10(abs(err_pk(:,1)));
err_ap_am_van = log10(abs(err_pk(:,2)));
err_ap_sd_van = log10(abs(err_pk(:,3)));

err_ap_exp_bic = log10(abs(err_ap_BIC(:,1)));
err_ap_off_bic = log10(abs(err_ap_BIC(:,end)));
err_ap_cf_bic = log10(abs(err_pk_BIC(:,1)));
err_ap_am_bic = log10(abs(err_pk_BIC(:,2)));
err_ap_sd_bic = log10(abs(err_pk_BIC(:,3)));

%% Code for comparing output errors from ms-specparam and default specparam
c1 = [1 50 150]./255;
c2 = [239 161 7]./255;
figure('Position',[400 400 700 400]),hold on
fill([-20 -20 -21],[-20 -21 -20],c1,'FaceAlpha',0.4)
fill([-20 -20 -21],[-20 -21 -20],c2,'FaceAlpha',0.4)

spc = 0.02;
sz = 5;

box_plot(err_ap_exp_van,c1,1,sz,spc,1);
box_plot(err_ap_exp_bic,c2,1,sz,spc,0);
box_plot(err_ap_off_van,c1,2,sz,spc,1);
box_plot(err_ap_off_bic,c2,2,sz,spc,0);
box_plot(err_pk_cf_van,c1,3,sz,spc,1);
box_plot(err_pk_cf_bic,c2,3,sz,spc,0);
box_plot(err_pk_am_van,c1,4,sz,spc,1);
box_plot(err_pk_am_bic,c2,4,sz,spc,0);
box_plot(err_pk_sd_van,c1,5,sz,spc,1);
box_plot(err_pk_sd_bic,c2,5,sz,spc,0);

xlim([0 6])
yticks([-4 -3 -2 -1 0 1])
yticklabels({'10^{-4}','10^{-3}','10^{-2}','0.1','1','10'})
ylim([-3.2 0.8])
ylabel('|error|')
xticks([1 2 3 4 5])
row1 = {'Exponent','Offset','Center frequency','Amplitude','Standard deviation'};
row2 = {'(Hz^{-1})','(a.u.)','(Hz)','(a.u.)','(Hz)'};
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels(tickLabels); 
text(0.1,0.2,'Noise = 0.15')
legend('default specparam','ms-specparam','Location','Northwest')

%% plot sensitivity intervals
flipxy = 1;
c3 = [21 102 5]./255;
figure('Position',[400 400 600 300]), hold on
fill([-20 -20 -21],[-20 -21 -20],c1,'FaceAlpha',0.4)
fill([-20 -20 -21],[-20 -21 -20],c3,'FaceAlpha',0.4)
if ~flipxy
    for k = 1:4
        draw_95eci(sens(k),sens_bs(:,k),-k+0.1,flipxy,c1)
        draw_95eci(sens_bic(k),sens_bic_bs(:,k),-k-0.1,flipxy,c3)
    end
else
    for k = 1:4
        draw_95eci(sens(k),sens_bs(:,k),((-flipxy)*-k)-0.1,flipxy,c1)
        draw_95eci(sens_bic(k),sens_bic_bs(:,k),((-flipxy)*-k)+0.1,flipxy,c3)
    end
end
if ~flipxy
    text(0.51,-3.7,'Noise = 0.15')
    xlim([0.5 1.05])
    xlabel('Sensitivity')
    yticks(-4:-1)
    yticklabels({'4 peaks','3 peaks','2 peaks','1 peak'})
    legend('default specparam','ms-specparam','Location','Southwest')
else
    text(3.7,0.6,'Noise = 0.15')
    ylim([0.5 1.05])
    ylabel('Sensitivity')
    xticks(1:4)
    xticklabels({'1 peaks','2 peaks','3 peaks','4 peak'})
    legend('default specparam','ms-specparam','Location','Southeast')
end


%% plot PPV intervals
flipxy = 1;
c3 = [21 102 5]./255;
figure('Position',[400 400 600 300]), hold on
fill([-20 -20 -21],[-20 -21 -20],c1,'FaceAlpha',0.4)
fill([-20 -20 -21],[-20 -21 -20],c3,'FaceAlpha',0.4)
if ~flipxy
    for k = 1:4
        draw_95eci(ppv(k),ppv_bs(:,k),-k+0.1,flipxy,c1)
        draw_95eci(ppv_bic(k),ppv_bic_bs(:,k),-k-0.1,flipxy,c3)
    end
else
    for k = 1:4
        draw_95eci(ppv(k),ppv_bs(:,k),((-flipxy)*-k)-0.1,flipxy,c1)
        draw_95eci(ppv_bic(k),ppv_bic_bs(:,k),((-flipxy)*-k)+0.1,flipxy,c3)
    end
end

if ~flipxy
    text(0.265,-3.7,'Noise = 0.15')
    xlim([0.25 1.05])
    xlabel('Positive predictive value')
    yticks(-4:-1)
    yticklabels({'4 peaks','3 peaks','2 peaks','1 peak'})
    legend('default specparam','ms-specparam','Location','Southwest')
else
    text(3.7,0.4,'Noise = 0.15')
    ylim([0.25 1.05])
    ylabel('Positive predictive value')
    xticks(1:4)
    xticklabels({'1 peaks','2 peaks','3 peaks','4 peak'})
    legend('default specparam','ms-specparam','Location','Southeast')
end


%% Plot recovery rates by number of simulated peaks, for default specparam and ms-specparam, as matrix
data = npeak_bic_mat; % swap with npeak_mat
c3 = [21 102 5]./255;
c2a = [249 217 156]./255;
c1a = [153 173 213]./255;
cx = c3;
figure('Position',[400 400 500 430]), hold on
xs = [0 1 2 3 4 5 6];
ys = [0 1 2 3 4];
sz = [data(1,:)./sum(data(1,:)) data(2,:)./sum(data(2,:)) data(3,:)./sum(data(3,:)) data(4,:)./sum(data(4,:)) data(5,:)./sum(data(5,:))];
reshape(sz,[7 5]);
rto = ((6.85./4.85));
imagesc(ys,xs,data')
for x = 1:7
   for y = 1:5
      text(ys(y),xs(x),num2str(round(data(y,x)./1000,3).*100),'HorizontalAlignment','center','VerticalAlignment','middle') 
   end
end
yticks(0:6)
xticks(0:4)
xlim([-0.5 4.5])
ylim([-0.5 6.5])
ylabel('Number of detected peaks')
xlabel('Number of simulated peaks')
colormap([linspace(1,cx(1),1000)' linspace(1,cx(2),1000)' linspace(1,cx(3),1000)']);
caxis([0 1000])
c = colorbar('Ticks',[0,250,500,750,1000],...
         'TickLabels',{'0','25','50','75','100'});
c.Label.String = 'Density (%)';
set(gca,'FontSize',14)

%% Functions to help code run
function draw_95eci(est,bootstrap,xloc,flipxy,colour)
    
    prcs = [prctile(bootstrap,2.5) prctile(bootstrap,97.5)];
    if ~flipxy
        scatter(est,xloc,30,colour,'filled')
        plot(prcs,[xloc xloc], 'Color', colour)
        text(prcs(2)+0.02,xloc,num2str(round(est,2)),'Color',colour)
    else
        scatter(xloc,est,30,colour,'filled')
        plot([xloc xloc],prcs, 'Color', colour)
        text(xloc,prcs(2)+0.02,num2str(round(est,2)),'Color',colour,'HorizontalAlignment','center')
    end
    
end

function box_plot(X, cl, offset,w,spc,flipx)

quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
Xs          = sort(X);
whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
if flipx
   spc = -spc; 
end

plot([offset offset]-(-1+2.*flipx).*w.*0.025+spc, whiskers,'Color',cl,'LineWidth',0.2)
fill([offset offset offset-(-1+2.*flipx).*w.*0.05 offset-(-1+2.*flipx).*w.*0.05]+spc,quartiles([1 2 2 1]),[1 1 1],'EdgeColor',cl)
plot([offset offset-(-1+2.*flipx).*w.*0.05]+spc,ones(1,2).*quartiles(3),'Color',cl,'LineWidth',0.2)
scatter(offset-(-1+2.*flipx).*w.*0.025+spc,mean(X),w.*3,cl,'Filled')

end

%% model functions
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

function noise = gen_noise(xs,sd)
%       Generate noise, without a 'knee'.
%
%       Parameters
%       ----------
%       xs : mxn matrix
%           Input spectrum matrix. Used for sizing only
%       sd : double
%           standard deviation of noise (used for SNR calculations)
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential (no-knee) function.

        noise = randn(size(xs)).*sd;

end

function spec = gen_spectrum(freqs, ap_params, peak_params, sd, peak_type, nS)
%       Generate spectra (in log10-power) with known parameters and zero-
%       mean gaussian-distributed noise.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       ap_params : 1x2(3) array (offset[, knee], exp )
%           Parameters that describe aperiodic component.
%       peak_params : kx3 array (mu, hgt, sigma [gamma])
%           Parameters that define k peaks (centre frequency,
%           height, and standard deviation [or gamma]).
%       sd : double
%           Standard deviation of noise (used for SNR calculations).
%       peak_type : 'gaussian' or 'cauchy'
%           Type of peak modelled in spectra.
%       nS : integer double
%           Number of simulations to be generated using the parameters
%           provided.
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential (no-knee) function.

    % Generate model aperiodic
    if length(ap_params) == 2
        spec = expo_nk_function(freqs,ap_params);
    elseif length(ap_params) == 3
        spec = expo_function(freqs,ap_params);
    else
        error('Invalid aperiodic parameter set.')
    end
    
    % determine peak model
    switch peak_type
        case 'gaussian'
            peak_function = @gaussian;
        case 'cauchy'
            peak_function = @cauchy;
        otherwise
            error('Invalid peak type.')
    end
    
    % Add peaks (if any exist)
    for pk = 1:size(peak_params,1)
        spec = spec + peak_function(freqs, peak_params(pk,1), peak_params(pk,2), peak_params(pk,3));
    end
    
    % Copy spectrum nS times
    spec = repmat(spec,[nS 1]);
    
    % Add zero-mean gaussian-distributed noise
    noise = gen_noise(spec,sd);
    
    % Combine models with noise
    spec = spec + noise;

end

function [outStruct] = gen_rand_spec_par(ap_ranges, pk_cf_range,...
    pk_amp_range, pk_sd_range, min_dist, n_peaks, n_sim)

% ap_ranges -           2x2 [exp_min exp_max; off_min off_max]
% pk_cf_range  -        1x2 [pk_cf_min pk_cf_max]
% pk_amp_range -        1x2 [pk_amp_min pk_amp_max]
% pk_sd_range -         1x2 [pk_sd_min pk_sd_max]
% min_dist -            Minimum distance from other peaks, in sd of other peaks
% max_peaks -           Maximum number of unique peaks in a simulation
% n_sim -               Number of unique simulations

    outStruct = struct();

    for s = 1:n_sim
        % construct aperiodic parameters
        ap_init = round(ap_ranges(:,1)+diff(ap_ranges')'.*rand(size(ap_ranges,1),1),2)';
        pk_cf = zeros(n_peaks,1); pk_amp = zeros(n_peaks,1); pk_sd = zeros(n_peaks,1); 
        for pk = 1:n_peaks
            % Ensure peaks don't overlap in frequency space.
            if pk == 1
                pk_cf(pk) = round(pk_cf_range(1)+diff(pk_cf_range).*rand,1);
                pk_amp(pk) = round(pk_amp_range(1)+diff(pk_amp_range).*rand,2);
                pk_sd(pk) = round(pk_sd_range(1)+diff(pk_sd_range).*rand,1);
            else
                pk_amp(pk) = round(pk_amp_range(1)+diff(pk_amp_range).*rand,2);
                pk_sd(pk) = round(pk_sd_range(1)+diff(pk_sd_range).*rand,1);
                pk_cf(pk) = round(pk_cf_range(1)+diff(pk_cf_range).*rand,1);
                zn = [pk_cf(1:pk-1)-min_dist.*pk_sd(1:pk-1) pk_cf(1:pk-1)+min_dist.*pk_sd(1:pk-1)];
                ctr = 1;
                overlaps_f = (zn(:,1) < pk_cf(pk) & pk_cf(pk) < zn(:,2));
                while any(overlaps_f)
                    ctr = ctr+1;
                    if ctr > 100
                        error('Randomization process failed: peaks overlap') 
                    end
                    pk_cf(pk) = pk_cf_range(1)+diff(pk_cf_range).*rand;
                    overlaps_f = (zn(:,1) < pk_cf(pk) & pk_cf(pk) < zn(:,2));
                end                   
            end
        end
        
        outStruct(s).ap_pars = ap_init;
        outStruct(s).peak_pars = [pk_cf pk_amp pk_sd];
        
    end
end

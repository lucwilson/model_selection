function opt = build_opt(varargin)
% build_opt: Helper function to build options struct for ms-specparam.
% Inputs: 
%          - optional: field-value pairs (eg. 'max_peaks',6)
%          - optional: partially completed opt struct
%          - List of available fields: freq_range, peak_width_limits,
%          max_peaks, min_peak_height, aperiodic_mode, peak_threshold,
%          cf_bound, border_threshold, return_spectrum,
%          proximity_threshold, optim_obj, peak_type, guess_weight, hOT
%          N.B. Order does not matter, as long as field-value pairs are
%          connected (i.e, struct not in between).
% Ouputs:
%          - Completed options struct, with presets unless specified above.

    % Pre-populate options with presets
    opt.freq_range          = [1 40];
    opt.peak_width_limits   = [0.5 12];
    opt.max_peaks           = 6;
    opt.min_peak_height     = 1 / 10;
    opt.aperiodic_mode      = 'fixed';
    opt.peak_threshold      = 2;
    opt.cf_bound            = 1.5;
    opt.border_threshold    = 1;
    opt.return_spectrum     = 0;
    opt.power_line          = '-5';
    opt.proximity_threshold = 0.75;
    opt.optim_obj           = 'negloglike';
    opt.peak_type           = 'gaussian';
    opt.guess_weight        = 'none';
    opt.hOT                 = any(any(contains(struct2cell(ver), 'Optimization Toolbox')));

    % Check if user has provided any modifiers
    nIn = length(varargin);
    if nIn % First pass: if user provided partially complete options
        structs = strcmp(cellfun(@class, varargin, 'UniformOutput', false),'struct');
        if any(~mod(find(structs),2)) % Ensure struct not between field-values
            error('Unexpected input formatting.')
        end
        if sum(structs) > 1 % Too many structs provided (expects 1 max)
            error('Unexpected multitude of struct arrays (>1).')
        elseif sum(structs) == 1 % opt provided, fill blanks
            opt_usr = varargin{structs};
            f1 = fieldnames(opt_usr);
            f2 = fieldnames(opt);
            [~, odds, fills] = setxor(f1,f2);
            if ~isempty(odds) % Remove unknown fields, mistyped or non-existent
                disp(['build_opt: Removing ' num2str(length(odds)) ' unexpected fieldnames from user.'])
                opt_usr = rmfield(opt_usr,f1(odds));
            end
            if ~isempty(fills) % Add missing fields to user-provided options.
                disp(['build_opt: Adding ' num2str(length(fills)) ' expected fieldnames from template.'])
                for f = 1:length(fills)
                   opt_usr.(f2{fills(f)}) = opt.(f2{fills(f)}); 
                end
            end
            opt = opt_usr;
            varargin(structs) = [];
        end
    end
    nIn = length(varargin);
    if nIn % Second pass: if user provided field-value pairs
        if mod(nIn,2) % If a field/value is not accompanied by value/field
            error('Missing field-value pairing.')
        else
            f2 = fieldnames(opt);
            for k = 1:2:nIn-1
                field = varargin{k};
                value = varargin{k+1};
                if isempty(intersect(field,f2)) % Skip if field does not exist
                    warning(['Skipping unknown field: ' field '.'])
                    continue
                elseif ~strcmp(class(value),class(opt.(field))) % Skip if wrong format for field
                    warning(['Skipping field: ' field ' (invalid value format).'])
                    continue
                else
                    opt.(field) = value;
                end
            end
        end
    end
end

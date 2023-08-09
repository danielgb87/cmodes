function  [psd] = dFC_utils_caps_psd(cfc, TR)


% This function takes computes the power-spectral density of time-courses
% of CAP-to-frame correlations fro each subject. It uses the same structure
% as in the MATLAB tutorial on PSD computation using the FFT (see ref[1].)
% Also compute group stats.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   CFC  : a subject x K-caps cell with the CFC timecourse in them.
%   TR: repetition time of data.
% OUTPUTS
%   psd  : a structure with all the spectral characteristics of teh
%   signals including the frequencies; psd; and stats.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%[psd] = caps_analysis_psd(CFC, TR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 14-11-17

% Detect amount of templates, subjects, and frames per subject. Also detect the
% amount of clusters  (k) in the analysis.

% function to compute windowed PSD using the FFT as in the matlab tutorial.
[Nsubs, k] = size(cfc);
Fs = 1/(TR);

for sub = 1:Nsubs
    for c = 1:k
        x = cfc{sub,c};
        Nobs = length(x);

        xdft = fft(x);
        xdft = xdft(1:Nobs/2+1);
        psd.psd{sub,c} = (1/(Fs*Nobs)) * abs(xdft).^2;
        psd.psd{sub,c}(2:end-1) = 2*psd.psd{sub,c}(2:end-1);
    end
    
    psd.freqs{sub,1} = 0:Fs/length(x):Fs/2;
    clear x
end

% compute statsbut first, if there are subjects with different amount of
% observations (timepoints), downsample them to the subject with the lowest
% amount of observations.
for sub = 1:Nsubs
    len_f(sub) = length(psd.freqs{sub});
end
[nf, sub_min] = min(len_f);

for c = 1:k
    psd_temp{c} = zeros(Nsubs,length(psd.freqs{sub_min}));
    for sub = 1:Nsubs
        if length(psd.freqs{sub}) ~= length(psd.freqs{sub_min})
                psd_temp{c}(sub,:) = interp1(psd.freqs{sub},psd.psd{sub,c},psd.freqs{sub_min});
        else
            psd_temp{c}(sub,:) = psd.psd{sub,c};
        end
    end
    for f = 1:length(psd.freqs{sub_min})
        psd_temp_mean{c}(f) = mean(psd_temp{c}(:,f));
        psd_temp_std{c}(f) = std(psd_temp{c}(:,f));
    end
end

% put psd_temp downsampled in the format of the original psd.
for sub = 1:Nsubs
    for c = 1:k
        psd.psd_ds{sub,c} = spm_vec(psd_temp{c}(sub,:));
    end
end
psd.min_freq_vec = psd.freqs{sub_min};
psd.psd_ds_mean = psd_temp_mean;
psd.psd_ds_std = psd_temp_std;

    
  
        

% end function.
end


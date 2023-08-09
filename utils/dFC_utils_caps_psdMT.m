function  [psd_MT] = dFC_utils_caps_psdMT(cfc, TR, NW)


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
        NFFT = length(x);
        [psd_MT.psd{sub,c},psd_MT.freqs{sub,c}] = pmtm(x,NW,NFFT,Fs);
        
    end
    
    clear x
end

% compute statsbut first, if there are subjects with different amount of
% observations (timepoints), downsample them to the subject with the lowest
% amount of observations.
for sub = 1:Nsubs
    for c = 1:k
        len_f(sub,c) = length(psd_MT.freqs{sub,c});
    end
end
[nf, sub_min] = min(min(len_f));

for c = 1:k
    psd_temp{c} = zeros(Nsubs,length(psd_MT.freqs{sub_min,c}));
    for sub = 1:Nsubs
        if length(psd_MT.freqs{sub,c}) ~= length(psd_MT.freqs{sub_min,c})
                psd_temp{c}(sub,:) = interp1(psd_MT.freqs{sub,c},psd_MT.psd{sub,c},psd_MT.freqs{sub_min,c});
        else
            psd_temp{c}(sub,:) = psd_MT.psd{sub,c};
        end
    end
    for f = 1:length(psd_MT.freqs{sub_min})
        psd_temp_mean{c}(f) = mean(psd_temp{c}(:,f));
        psd_temp_std{c}(f) = std(psd_temp{c}(:,f));
    end
end

% put psd_temp downsampled in the format of the original psd.
for sub = 1:Nsubs
    for c = 1:k
        psd_MT.psd_ds{sub,c} = spm_vec(psd_temp{c}(sub,:));
    end
end
psd_MT.min_freq_vec = psd_MT.freqs{sub_min};
psd_MT.psd_ds_mean = psd_temp_mean;
psd_MT.psd_ds_std = psd_temp_std;

    
  
        

% end function.
end


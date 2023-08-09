function  [X_filt] = dFC_utils_caps_NBfilter(X,f_low,f_up, FS)


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
[Nsubs, k] = size(X);

 d1 = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',f_low,'CutoffFrequency2',f_up,'SampleRate',FS);
for sub = 1:Nsubs
    for c = 1:k
        X_filt{sub,c} = filtfilt(d1,double(X{sub,c}(:)));
        X_filt{sub,c} = zscore(X_filt{sub,c}(:));
    end
end 
        

% end function.
end


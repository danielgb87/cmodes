function  [X_filt] = caps_analysis_NBFilter(X,f_low,f_up, FS)


% This function tfilters signals in a band using a custom dseigned filter
% for fMRI signals (order 20, double pass).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 23-03-2020

% Detect amount of templates, subjects, and frames per subject. Also detect the
% amount of clusters  (k) in the analysis.

% function to compute windowed PSD using the FFT as in the matlab tutorial.

 d1 = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',f_low,'CutoffFrequency2',f_up,'SampleRate',FS);

        X_filt{sub,c} = filtfilt(d1,double(X{sub,c}(:)));
        X_filt{sub,c} = zscore(X_filt{sub,c}(:));
    end
end 
        

% end function.
end


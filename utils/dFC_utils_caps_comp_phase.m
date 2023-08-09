function [phase] = dFC_utils_caps_comp_phase(X)

% This function computes theb instantaneous phase of the input signals.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   peaks_sub:  cell (Nsubjects x Ncaps) with the signals to be analysed.

% OUTPUTS
%   phase: structure with the phase angle; amplitude and analytical signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
% [phase] = caps_analysis_comp_phase(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 17-11-17

[Nsubs, k] = size(X);

for sub = 1:Nsubs
    for c = 1:k
            phase.asignal{sub,c}=hilbert(spm_vec(X{sub,c}));
            phase.phase{sub,c}=angle(hilbert(spm_vec(X{sub,c})));
            phase.amplitude{sub,c}=abs(hilbert(spm_vec(X{sub,c})));
                        
    end

end

%end function
end
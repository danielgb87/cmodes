function [dec_out, data_dec] = postproc_deconvolve(data_norm, params, temporal_mask)

% This function deconvolves data with a reference hemodynamic response function (HRF) as described in ref[1].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%        data_norm    : amtrix of time x nodes.
%        params       : 
        %         params.TR -  repetition time.
        %         params.T  -  temporal grid to subsample (params.T=3, 3x finer grid).
        %         params.T0   -  position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then params.T0=fix(para.T/2)
        %         params.dt   -  params.dt  = params.TR/params.T;  fine scale time resolution.
        %         params.TD_DD = 2; % time and dispersion derivative
        %         params.AR_lag = 1; % AR(1) noise autocorrelation.
        %         params.thr = 1; % (mean+) params.thr*standard deviation threshold to detect event.
        %         params.len = 24; % length of HRF, here 24 seconds
        %         params.lag  = fix(3/params.dt):fix(9/params.dt); % 3 to 9 seconds

%        temporal_mask : temporal_mask = []; % without mask, it means temporal_mask = ones(nobs,1); i.e. all time points included. nobs: number of observation = size(data,1). if want to exclude the first 1~5 time points, let temporal_mask(1:5)=0;

% OUTPUTS
%   data_out   : structure for with the deconvolution information and the voxel-wise HRF.
%   data_dec   : cell struct with the data deconvolved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%[dec_out, data_dec] = postproc_deconvolve(data_norm, params,temporal_mask);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REFS. 
% [1] Guo-Rong Wu, Wei Liao, Sebastiano Stramaglia, Ju-Rong Ding,
% Huafu Chen, Daniele Marinazzo*. "A blind deconvolution approach to
% recover effective connectivity brain networks from resting state fMRI
% data." Medical Image Analysis 17:365-374

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 14-11-17


[dec_out.beta_hrf dec_out.bf dec_out.event_bold] = wgr_rshrf_estimation_canonhrf2dd_par2(data_norm,params,temporal_mask);

hrfa = dec_out.bf*dec_out.beta_hrf(1:size(dec_out.bf,2),:); %HRF
 
[nobs nvar] = size(data_norm); PARA = zeros(3,nvar);
for i=1:nvar 
hrf1 = hrfa(:,i); 
dec_out.hrf = hrfa;
[dec_out.PARA(:,i)] = wgr_get_parameters(hrf1,params.TR/params.T);% estimate HRF parameter 	
end

% now deconvolve data timecourses
T = round(params.len/params.TR);
for voxel_id=1:nvar
        hrf=hrfa(:,voxel_id);
        H=fft([hrf; zeros(nobs-length(hrf),1)]);
        M=fft(data_norm(:,voxel_id));
        data_dec(:,voxel_id) = ifft(conj(H).*M./(H.*conj(H)+2));
end













end
function   [motion_info] = dFC_utils_check_motion(path, suffix, radius, divide, suite, scrubFD_thr)


% This function uses a modified version of the function
% "bramila_framewiseDisplacement.m" bramilla toolbox (ref [1]), see
% motion_FD_comp. https://version.aalto.fi/gitlab/BML/bramila
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   inputs
%    path  = path were the motion traces .txt are '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/07_motion_correction'
%    suffix: the suffix of files to select: ex:'*GqNegKNeg*_*pre_cno_bs*_*.txt'
%    radius           = 5;      % mouse brain radius (5mm by default)
%    divide           = 1;      % divide by 10 translation motion traces.
%    suite            = 'afni'; % suite used to obtain motion parameters ('fsl', 'afni', or 'spm')
%    scrubFD_thr        = 0.05;   % Framewise displacement censoring criteria (in mm).

% OUTPUTS
%   motion_info: structure with information about motion for the selected
%   dataset including framewise displacement and DVARS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REFS.
% [1] MATLAB code (BRAMILA pipeline v2.0, available at https://git.becs.aalto.fi/bml/bramila/).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1, 25/02/2020

cd(path)
list = dir(suffix);
Nsubs  = length(list);

for sub = 1:Nsubs
    sub_list{sub,1} = list(sub).name;
end
Nsubs     = length(sub_list);
radius    = radius;
suite     = suite;
FD_thr    = scrubFD_thr;

%% 1. Assuming that the motion traces are in the 'nuisance_regresion' folder of each subject's preprocessing, extract the motion traces.

for sub = 1:Nsubs
    motion_raw{sub,1} = load(sub_list{sub});
    
    % Divide translation traces by 10
    if divide == 1
        motion_raw{sub,1}(:,4:6) = motion_raw{sub,1}(:,4:6)/10;
    end
    
    % Compute framewise displacement; differential motion traces,
    % root-mean-square deviation and scaled motion traces  (in mm and
    % degrees).
    
    [FD{sub,1}, dMT{sub,1}, RMS{sub,1}, TS_dof{sub,1}] = dFC_utils_motion_FD_comp(motion_raw{sub,1}, radius, suite);
    
end
% Partial summary
motion_info.raw_motion = motion_raw;
motion_info.FD = FD;
motion_info.dMT = dMT;
motion_info.RMS = RMS;
motion_info.TS_dof = TS_dof;



%% 2. Compute some statistics; count above FD threshold volumes; and put motion warnings on flagged subjects.

motion_info.FD_warnings = zeros(Nsubs,1);
for sub = 1:Nsubs
    motion_info.FD_count(sub,1) = sum(FD{sub,1} > FD_thr);
    motion_info.FD_flag{sub,1} = zeros(length(FD{sub,1}),1);
    for t = 1:length(FD{sub,1})
        if FD{sub,1}(t) > FD_thr
            motion_info.FD_flag{sub,1}(t,1) = 1;
        end
    end
    
    motion_info.FD_flag_proportion(sub,1) = sum(motion_info.FD_flag{sub,1})/length(motion_info.FD_flag{sub,1});
    % flag subjects
    if  motion_info.FD_flag_proportion(sub) > 0.5
        motion_info.FD_warnings(sub) = 1;
    end
    
end

motion_info.FD_thr = FD_thr;





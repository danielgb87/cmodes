function    [FD, dMT, RMS, TS_dof] = dFC_utils_motion_FD_comp(x, radius, suite)

% This function is a modified version of the function
% "bramila_framewiseDisplacement.m" bramilla toolbox (ref [1]), see
% motion_FD_comp. https://version.aalto.fi/gitlab/BML/bramila

% The function takes the motion traces x(six dimensions, 3 rotations, 3
% translations depending on the suite), and a radius for the brain (5mm for
% mice) and computes Framewise displacement; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%  x       : motion parameter traces (time x 6 Degrees of freedom)
%  radius  : estimated radius of the brain (in mm).
%  suite   : suite used for motion estimation ('afni', 'fsl', or 'spm') in
%  order nto see if rotations are first or translations are
    
% OUTPUTS
%  FD    : Framewise displacement timecourse.
%  dMT   : abvsolute derivatives of raw motion
%  RMS   : root-mean square of each motion trace
%  TS_dof: converted motion traces with 3 first rotations (degrees), then
%  translations in mm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%[FD, dMT, RMS, TS_dof] = motion_FD_comp(x, radius, suite);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REFS. 
% [1] MATLAB code (BRAMILA pipeline v2.0, available at https://git.becs.aalto.fi/bml/bramila/). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 14-11-17

    ts = double(x);
    r = radius;
    prepro_suite = suite;
    
    % Assure that the motion parameter file has 6 columns.
    if(size(ts,2)~=6)
        error(['The motion traces must have 6 motion parameters in 6 columns; the size of the input given is ' num2str(size(ts))])
    end
    
    % Define preprocessing suite... if not define, report error
            
    if(strcmp(prepro_suite,'fsl'))
        % convert radians into motion in mm
		% in FSL the first 3 columns are rotations, and they come in
		% radians.
        temp=ts(:,1:3);
        temp=r*temp;
        ts(:,1:3)=temp;
        
        TS_dof(:,4:6) = ts(:,4:6);
        TS_dof(:,1:3) = ts(:,1:3)*180/pi;  % turn to degrees.
                
    elseif(strcmp(prepro_suite,'afni'))
        % convert degrees into radians and then motion in mm
 		% in AFNI the first 3 columns are rotation, and come in degrees, so
        % for the FD calculation, they must go into radians and then
        % re-scaled to degrees for TS_dof saving.
        temp=ts(:,1:3);
        temp=r*temp*pi/180;
        ts(:,1:3)=temp;   
        
        TS_dof(:,4:6) = ts(:,4:6);
        TS_dof(:,1:3) = ts(:,1:3)*180/pi;  % turn to degrees.
        
    elseif (strcmp(prepro_suite,'spm'))
        % SPM way
        % SPM rotations come also in radians, but rotations are the last 3 columns;
        temp=ts(:,4:6);
        temp2=ts(:,1:3);
        temp=r*temp;
        ts(:,1:3)=temp;
        ts(:,4:6) = temp2;
        
        TS_dof(:,4:6) = ts(:,4:6);
        TS_dof(:,1:3) = ts(:,1:3)*180/pi;  % turn to degrees.
        
    else
        error('Please define the preprocessing suite as inputs.pproc.motion_suite. It can be either afni, fsl, or spm.')
    end
    
    dts=diff(ts);
    dts=[
        zeros(1,size(dts,2)); 
        dts
        ];  % first element is a zero, as per Power et al 2014
    
    dMT = abs(dts);
    FD=sum(abs(dts),2);
    RMS=sqrt(mean(ts.^2,1));    % root mean square for each column
    
end

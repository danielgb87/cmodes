function   data = dFC_pproc_normalize_data(data)

% This function normalizes data to standard deviation units. For each
% voxel's time-course X(t), the function delivers the temporal vector of
% z-scores (X-MEAN(X)) ./ STD(X).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   data   : cell structure with a cell for every subject and data in
%   matrix form (time x nodes(voxels))
% OUTPUTS
%   data   : cell structure with a matrix (time x nodes(voxels)) for each subject with data normalized temporally.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%data = postproc_normalize(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 25-01-2020

if iscell(data) == 0
        error('error - Data must be in cell format')
end

[Nvols Nvox] = size(data{1});
for sub = 1:length(data)
    for i = 1:Nvox
        data{sub,1}(:,i) = zscore(data{sub}(:,i));
    end
end

end %function
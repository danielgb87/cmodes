function [result, tmpFile] = chkzip(fileName)
% Function from:
% https://it.mathworks.com/matlabcentral/fileexchange/36193-chkzip--check-whether-a-file-is-a-compressed-archive
% CHKZIP  check whether file is a compressed archive
%
%  chkzip checks whether the passed file has one of the following
%    extensions: gz, zip, tar, tgz, tar.gz.  If it does, chkzip
%    attempts to unpack the archive into a temporary location,
%    which is platform-dependent.
%
% Usage:
%   [result, tmpFile] = chkzip(fileName)
%
%     fileName: file to be checked
%       result: number of files unpacked from compressed file,
%                 or 0 if file extension was not recognized
%      tmpFile: full path to first extracted file if
%                 unpacking was performed, otherwise returns
%                 original filename
%
% See also:  GZIP, ZIP, TAR

% v0.1 (Apr 2012) by Andrew Davis:  addavis@gmail.com

% Check input filename
assert(exist(fileName, 'file') == 2, 'file not found: %s', fileName);

% Create temporary directory and unpack
tmpDir = tempname;
mkdir(tmpDir);
if strcmp(fileName(end-3:end), '.tar') || strcmp(fileName(end-3:end), '.tgz') || strcmp(fileName(end-6:end), '.tar.gz'),
   tmpFiles = untar(fileName, tmpDir);
elseif strcmp(fileName(end-3:end), '.zip'),
   tmpFiles = unzip(fileName, tmpDir);
elseif strcmp(fileName(end-2:end), '.gz'),
   tmpFiles = gunzip(fileName, tmpDir);
else,
   tmpFile = fileName;
   rmdir(tmpDir);
   result = 0;
   return;
end;

tmpFile = tmpFiles{1};
result = length(tmpFiles);
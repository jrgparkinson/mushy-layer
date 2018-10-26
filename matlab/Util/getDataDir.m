function [ data_dir ] = getDataDir( folder )
%GETDATADIR Summary of this function goes here
%   Detailed explanation goes here

if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of 
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
end

user_name = strcat(user_name);

if strcmp(user_name, 'parkinsonjl')  
    base_dir = '/home/parkinsonjl/mnt/sharedStorage/';
else
   base_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/'; 
end

%base_dir = '/media/parkinsonjl/DATA/';
if nargin < 1   
    folder = 'optimalStates-highRes-new/';
end

data_dir = fullfile(base_dir, folder); %[base_dir, folder];


%data_dir = ['', ...
%    'optimalStates-highRes-new/'];

end


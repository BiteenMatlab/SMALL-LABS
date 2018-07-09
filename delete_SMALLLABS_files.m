function delete_SMALLLABS_files(directory_name,del_avgsub,del_guesses,del_offframes,del_fits,del_viewfits)
%% SMALLLABS_file_deletion
% This function deletes the various .mat data files created by SMALL-LABS
% for storage management.
%
% directory_name is the name of the directory to search for the files
%
%%% The following optional parameters are Booleans determining whether or
%%% not to delete the various kinds of files. Set to true to delete, and to
%%% false to not delete. Default is true for all.
%
% del_avgsub delete avgsub movies
% del_guesses delete the guess .mat files
% del_offframes delete the Mol_off_frames .mat files
% del_fits delete the fits .mat files
% del_viewfits delete the ViewFits movies

%% Boolean defaults
if ~exist('del_avgsub','var')
    del_avgsub=true;
end
if ~exist('del_guesses','var')
    del_guesses=true;
end
if ~exist('del_offframes','var')
    del_offframes=true;
end
if ~exist('del_fits','var')
    del_fits=true;
end
if ~exist('del_viewfits','var')
    del_viewfits=true;
end

%% Find the files

%prepping the directory_name
if directory_name(end)=='\'
    directory_name(end)=[];
end

if del_avgsub
    tfiles=dir([directory_name,filesep,'*_avgsub.mat']);
    %concatenate the file structures
    if ~exist('files','var')
        files=tfiles;
    else
        files=[files;tfiles];
    end
    
    %for legacy SMALL-LABS with avgsub.tif
    tfiles=dir([directory_name,filesep,'*_avgsub.tif']);
    %concatenate the file structures
    if ~exist('files','var')
        files=tfiles;
    else
        files=[files;tfiles];
    end
end

if del_guesses
    tfiles=dir([directory_name,filesep,'*_guesses.mat']);
    %concatenate the file structures
    if ~exist('files','var')
        files=tfiles;
    else
        files=[files;tfiles];
    end
end

if del_offframes
    tfiles=dir([directory_name,filesep,'*_Mol_off_frames.mat']);
    %concatenate the file structures
    if ~exist('files','var')
        files=tfiles;
    else
        files=[files;tfiles];
    end
end

if del_fits
    tfiles=dir([directory_name,filesep,'*_fits.mat']);
    %concatenate the file structures
    if ~exist('files','var')
        files=tfiles;
    else
        files=[files;tfiles];
    end
end

if del_viewfits
    tfiles=dir([directory_name,filesep,'*_ViewFits.avi']);
    %concatenate the file structures
    if ~exist('files','var')
        files=tfiles;
    else
        files=[files;tfiles];
    end
end


%% Delete the files

for ii=(1:length(files))
    delete([files(ii).folder,filesep,files(ii).name])
end


end


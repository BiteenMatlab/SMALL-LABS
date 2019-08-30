function Movie2mat(file_or_directory_name)
%% Movie2mat
% Converts movie(s) to a .mat file containing the movie as a variable
% called mov
%
% Currently fully supported formats: .tif
%
% However, Movie2mat will attempt to the Bio-Formats library to open the
% movie, and if that fails will then try to use Matlab's built in
% VideoReader utility to read in any not specified format. Additionally,
% this program is intended to be modified to be supported for whatever
% format you need to use, simply add an elseif statement in the main
% section (where I've added the comment "add future elseif statements
% here")
%
%%%% Inputs %%%%
% file_or_directory_name is a string which is either the name of a file or
% a directory. If a filename, then Movie2mat converts that movie to a .mat
% file If a directory name, then Movie2mat opens up uigetfile to select the
% list of movies to be converted
%
%%%% Outputs %%%%
% no outputs, however Movie2mat saves a .mat file with the mov variable
%
%%%% Dependencies %%%%
% TIFFStack
% bfmatlab
%
%     Copyright (C) 2018  Benjamin P Isaacoff
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%% Setup
tic
% check if it's a directory or a file
if exist(file_or_directory_name)==7
    %Select movies with uigetfile. If you make an error in specifying the
    %directory, it opens in the current directory.
    disp('Select the movie(s)')
    try
        [datalist,dataloc,findex]=uigetfile([file_or_directory_name,filesep,'*.*'],'multiselect','on');
    catch
        curdir=pwd;
        [datalist,dataloc,findex]=uigetfile([curdir,filesep,'*.*'],'multiselect','on');
    end
    if findex==0
        error('No movies selected')
    end
    %convert to a list of directories and filenames
    if ~iscell(datalist); datalist={datalist}; end
    for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
    [dlocs,dnames,exts]=cellfun(@fileparts,datalist,'uniformoutput',false);
elseif exist(file_or_directory_name)==2
    %get the directory and filename, and format into the cell list as above
    [dname,fname,ext] = fileparts(file_or_directory_name);
    dlocs{1}=dname;
    dnames{1}=fname;
    exts{1}=ext;
    clear dname fname ext
else
    error('Please input either a directory name or a filename.')
end

% try to add bio-formats matlab directory to path
try
    addpath(genpath('bfmatlab'))
end

%% Read the movies and save the .mat file

for ii=1:numel(dlocs)
    %the full filename
    filename=[dlocs{ii},filesep,dnames{ii},exts{ii}];
    
    %check if it's a .mat file
    if strcmp(exts{ii},'.mat')
        isv73=getmatver([dlocs{ii},filesep,dnames{ii},exts{ii}]);
        if ~isv73
            error(['Please ensure that movie .mat files are version 7.3 and ',...
                'the movie data is saved in the ''mov'' variable.'])
        end
    else
        
        %setup the matfile
        mov='tempvariable';
        save([dlocs{ii},filesep,dnames{ii}],'mov','-v7.3');
        matio=matfile([dlocs{ii},filesep,dnames{ii},'.mat']);
        matio.Properties.Writable=true;
        
        %check how big the file is and compare to the available memory so
        %as to keep this moving quickly and avoid going into the swap the
        %size of the current movie
        try
        fmem=dir(filename);
        fmem=fmem.bytes;
        %the currently available memory
        [~,curmem]=memory;
        curmem=curmem.PhysicalMemory.Available;
        %determine a reasonable number of chunks to import by comparing the
        %filesize to 90% of the available memory
        numchunks=floor(fmem/(0.9*curmem))+1;
        catch
        numchunks=1;
        end
        %% go through and import the different filetypes
        %if it's a .tif stack use TIFFStack
        if strcmp(exts{ii},'.tif') || strcmp(exts{ii},'.tiff')
            %creat the TIFFStack object
            tfstk=TIFFStack(filename);
            movsz=size(tfstk);%the size of the movie
            
            %create the .mat file and the variable mov
            mov=zeros(movsz);%create the correct sized array
            %determine the image type by querying a single pixel of the
            %movie
            imtype=class(tfstk(1));
            %convert mov to the appropriate type
            mov=eval([imtype,'(mov)']);
            %overwrite the temp variable with the empty mov
            matio.mov=mov;
            
            chunklen=ceil(movsz(3)/numchunks);%the number of frames for each chunk
            for jj=1:numchunks
                if jj~=numchunks
                    matio.mov(:,:,((jj-1)*chunklen+1):(jj*chunklen))=tfstk(:,:,((jj-1)*chunklen+1):(jj*chunklen));
                else
                    matio.mov(:,:,((jj-1)*chunklen+1):end)=tfstk(:,:,((jj-1)*chunklen+1):end);
                end
            end
            
            % add future elseif statements here for instance: elseif
            % strcmp(ext{ii},'.h5')
        else
            try
                clear mov
                %try using the bio-formats reader
                mov=bfopen(filename);
                matio.mov=cat(3,mov{1,1}{:,1});
            catch
                try
                    clear mov
                    %try using Matlab's built in Matlab's built in
                    %VideoReader utility
                    vid = VideoReader(filename);    
                    %Read all video frames.
                    jj=1;%frame index
                    while hasFrame(vid)
                        curframe= readFrame(vid);
                        if numel(size(curframe))==3
                            curframe=mean(curframe,3);
                        end
                        mov(:,:,jj) = curframe;
                        jj=jj+1;
                    end
                    %save it
                    matio.mov=mov;
                    
                catch
                    error('Unable to read in the file. Please use a supported format, or update Movie2mat.mat for this format')
                end
            end
        end
    end
end

% try to remove bio-formats matlab directory to path
try
    rmpath(genpath('bfmatlab'))
end

%check if it's a version 7.3 matfile.
    function isv73 = getmatver(fname)
        x = evalc(['type(''', fname, ''')']);
        isv73 = strcmp(x(2:20), 'MATLAB 7.3 MAT-file');
    end

end

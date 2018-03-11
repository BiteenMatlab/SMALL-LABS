function SMALLLABS_main(file_or_directory_name,dfrlmsz,avgwin,moloffwin,varargin)
%% SMALLLABS_main
%
%%%% If you obtained this code from anywhere other than the Biteen Lab
%%%% GitHub please visit https://github.com/BiteenMatlab/SMALL-LABS to
%%%% obtain the most up-to-date version of this code.
%
% SMALLLABS_main is the wrapper function for the SMALL-LABS algorithm which
% does accurate background substraction of single molecule imaging movies
% so that the molecules can be fit and their intensity accurately measured.
%
% Note: that by setting bgsub=false SMALLLABS_main will do fitting without doing
% background subtraction. Doing this obviates a lot of parameters, it
% doesn't matter what they're set to. See the User Guide for more details.
%
%
%%%% Inputs %%%%
%%% required
%   file_or_directory_name   is the name of the directory where the movies
%   will be selected OR the filename of a single movie OR a list of movies
%   to be fit. See the User Guide for more details and/or Movie2mat for
%   details about the file formats currently supported.
%
%   dfrlmsz   is the size of a diffraction limited spot in pixels. It's the
%   nominal diameter, NOT the FWHM or something similar. Must be an
%   integer! For an expected diffraction limited standard deviation, std,
%   using the full width at 20% max, dfrlmsz = std*(2*sqrt(2*log(5)))
%
%   avgwin   is the length of the temporal window (in frames) to be used
%   for the average subtraction.
%
%   moloffwin   is the length of the temporal window (in frames) to be
%   checked to determine in which frames that molecule was off and are thus
%   safe to subtract. Needs to be an even integer.
%
%%% optional (varargin)
% See the list of optional parameters in parameters section below, details
% can be found in the User Guide. They are called with a name value pair as
% is standard for Matlab functions. e.g., 'bpthrsh',83 would set the
% parameter bpthrsh=83
%
%
%%%% Outputs %%%%
% This function will output a series of movies & .mat files depending on
% which functions are called. See each function or the User Guide for
% details
%
%%%% Dependencies %%%%
% AVGSUB_moves
% TIFFStack 
% Guessing
% bpass
% Mol_off_frames
% MLEwG
% gaussFit
% Subtract_then_fit
% Track_3D2
% Tracking
% Track_filter
% hungarian
% ViewFits
% ViewFitsTracking
% gpufit
% Movie2mat
%
%%%%
% Written by Benjamin P Isaacoff at the University of Michigan last update
% 3/10/18 BPI & SAL
%
%     Copyright (C) 2018  Benjamin P Isaacoff
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% Parameter Defaults
% You are of course welcome to change the default values, but I would
% strongly urge you to instead set them as inputs to the function using a
% name value pair (e.g., 'bpthrsh',83). Please see the User Guide for
% details about the meaning and function of these parameters. The default
% parameter values are:

%%% Actions %%%
% Do background subtraction
params.bgsub = true;
% Make the average subtracted movie
params.makeAvgsub = true;
% Do guessing
params.guessing = true;
% Check the guesses
params.checkGuesses = false;
% Make the off frames list
params.makeOffFrames = true;
% Do fitting
params.fitting = true;
% Do tracking
params.tracking = true;
% Make the ViewFits movie
params.makeViewFits = true;

%%% AVGSUB parameters %%%
% Subtract the temporal average? otherwise use median
params.do_avg = true;
% AVGSUB offset
params.offset = 1000;

%%% Guessing parameters %%%
% the threshold percentile of the banpassed movie
params.bpthrsh = 95;
% how many pixels to ignore around the edge of the frame
params.egdesz = dfrlmsz;
% compare brightnesses in each frame? otherwise use entire movie
params.pctile_frame = false;
% use a mask for guessing? Enter it's filename or set to true to look
% for .mat file w/ '_PhaseMask' appended to movie name
params.mask_fname=[];
% make a movie (.avi) of the guesses to check parameters
params.make_guessmovie = false;

%%% Fitting parameters %%%
% do MLE fitting? If not least squares fitting will be used
params.MLE_fit = false;
% Goodfit parameters. See Subtract_mol_off_frames for the details
params.stdtol = 1.5;
params.maxerr = 0.10; %if you change this please also change the if statement after the next loop
% subtract the mean of off frames? If not use median
params.do_avgsub = true;
% Which Gaussian function to fit to when using LSQ fitting? 1. symmetric,
% 2. fixed angle asymmetric, 3. free angle asymmetric
params.which_gaussian = 1;
% the angle to fit to for a fixed angle Gaussian
params.fit_ang= 0;
% use a GPU if one is a available?
params.usegpu= true;

%%% Tracking parameters %%%
% default is [0.01,200,0.5,3,3,1,0]
% save the separate tracks .mat file?
params.savetracks = false;
% minimum merit
params.trackparams(1)=0.01;
% Integration time (ms)
params.trackparams(2)=200;
% gamma
params.trackparams(3)=0.5;
% maximum step size
params.trackparams(4)=3;
% minimum track length
params.trackparams(5)=3;
% speed estimation window halfsize
params.trackparams(6)=1;
% time delay between consecutive frames (ms)
params.trackparams(7)=0;

%%% ViewFits parameters %%%
% use the original movie? if not, use the avgsub movie
params.orig_movie = true;
% diameter of the circles showing the fits
params.circ_D = dfrlmsz;
% linewidth of the circles
params.linewidth = 1;
% write a .avi movie showing the fits. If not, goes to debug mode
params.write_mov=true;
% autoscale frame by frame?
params.autoscale_on = false;
% use the tracking viewfits instead
params.trackingVF = false;

%% Evaluating the inputs

paramsnames=fieldnames(params);

% if any parameters are included as inputs, change the parameter mentioned
if nargin>4
    for ii=1:2:nargin-5
        whichField = strcmp(paramsnames,varargin{ii});       
        try
            eval(['params.' paramsnames{whichField} ' = varargin{ii+1};'])
        catch
            error([varargin{ii}, '  is not an input parameter. Check the spelling.'])
        end
    end
end

%changing the default error if doing MLE fitting
if params.MLE_fit && params.maxerr==3
    params.maxerr=0.1;
end

% verify that a GPU is available to use if usegpu==true
if params.usegpu
    params.usegpu=parallel.gpu.GPUDevice.isAvailable;
end
%% Checking the required inputs
% Erroring out if dfrlmsz isn't an integer, because it's a "strong"
% parameter. You should really input what you mean.
if dfrlmsz~=round(dfrlmsz);error('dfrlmsz must be an integer');end

%Rounding moloffwin and reseting it to proper values. Not
%erroring because it's a "weak" parameter. Only needed if doing bgsub
if params.bgsub && params.makeOffFrames
    if moloffwin~=(ceil(moloffwin/2)*2)
        moloffwin=(ceil(moloffwin/2)*2);
        warning(['moloffwin must be an even integer. It has been reset to avgwin = ',num2str(moloffwin)])
    end
end

%% Select the movies

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
    %if it's a .txt file then assume it's a list of filenames
    if strcmp(ext,'.txt')
        %open the file for reading
        fid=fopen(file_or_directory_name,'r');
        %initializing loop variables
        linetxt='foobar';
        ii=0;
        while all(linetxt~=-1)
            ii=ii+1;
            linetxt=fgetl(fid);
            if linetxt~=-1
                [dname,fname,ext]=fileparts(linetxt);
                dlocs{ii}=dname;
                dnames{ii}=fname;
                exts{ii}=ext;
            end
        end
        fclose(fid);%close the file
        %remove any leading or trailing blank spaces
        dlocs=strtrim(dlocs);
        exts=strtrim(exts);
    else
        dlocs{1}=dname;
        dnames{1}=fname;
        exts{1}=ext;
    end
    clear dname fname ext
else
    error('Please input either a directory name or a filename.')
end

%% Loop through all the movies

% turn off the warning if a variable isn't found in a .mat file
warning('off','MATLAB:load:variableNotFound');

% try to add gpufit to path
try
    addpath(genpath('gpufit'))
end

h2=waitbar(0);
set(findall(h2,'type','text'),'Interpreter','none');
waitbar(0,h2,'Starting SMALLLABS');
wholeshabang=tic;
for ii=1:numel(dlocs)
    clear goodframe
    try; waitbar((7*ii-6)/numel(dlocs)/7,h2); end
    
    %% Convert the movies to .mat
    % Movies must be version 7.3 mat files with the movie data saved in the
    % 'mov' variable. This function converts several standard scientific movie
    % data types into this format.
    Movie2mat([dlocs{ii},filesep,dnames{ii},exts{ii}])
    
    load([dlocs{ii},filesep,dnames{ii},'.mat'],'mov');
    mov=single(mov);
    movsz=size(mov);
    try
        load([dlocs{ii},filesep,dnames{ii},'.mat'],'goodframe');
    catch
    end
    if ~exist('goodframe','var')
        goodframe=true(movsz(3),1);
    end
    %% The Average Subtraction
    % Only if doing bgsub.     
    % Because this step can be somewhat slow, first there is a check to
    % determine if the avgsub movie has already been created. If so, then
    % there's a further check to determine whether or not any of the
    % parameters have changed. If they're identical, this step is skipped.
    % To bypass this process, simply delete the avgsub movie from the
    % directory.
    
    % AVGSUB_movs will save an average subtracted movie .mat file, called
    % moviename_avgsub.mat 
    if params.makeAvgsub && params.bgsub
        if exist([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'file')~=0
            avgsubparams=load([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'subwidth','offset','do_avg','goodframe');
            if avgsubparams.subwidth==avgwin && avgsubparams.offset==params.offset && avgsubparams.do_avg==params.do_avg && sum(avgsubparams.goodframe)==sum(goodframe)
                disp(sprintf(['Avgsub parameters unchanged from previous calculation.\r Skipping avgsub on ',dnames{ii}],1))
            else
                try; waitbar((7*ii-5)/numel(dlocs)/7,h2,{['Running AVGSUB on ',dnames{ii}],'Overall Progress'}); end
                bgsub_mov=AVGSUB_movs([dlocs{ii},filesep,dnames{ii}],mov,goodframe,...
                    params.do_avg,avgwin,params.offset);
            end
        else
            %do the avgsub
            try; waitbar((7*ii-5)/numel(dlocs)/7,h2,{['Running AVGSUB on ',dnames{ii}],'Overall Progress'}); end
            bgsub_mov=AVGSUB_movs([dlocs{ii},filesep,dnames{ii}],mov,goodframe,...
                params.do_avg,avgwin,params.offset);
        end        
    end
    
    %% Make Guesses
    % loop through each movie and make the guesses, which will be saved in a
    % .mat file called moviename_guesses.mat
    if params.guessing
        try; waitbar((7*ii-4)/numel(dlocs)/7,h2,{['Making guesses for ',dnames{ii}],'Overall Progress'}); end
        if params.bgsub
            % try loading in the bgsub movie
            try
                bgsub_mov=load([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'mov');
                bgsub_mov=bgsub_mov.mov;
            catch
                warning('No avgsub file was found, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
                keyboard
                params.bgsub=0;
            end
            if params.bgsub
                Guessing([dlocs{ii},filesep,dnames{ii},'_avgsub'],bgsub_mov,movsz,goodframe,dfrlmsz,...
                    params.bpthrsh,params.egdesz,params.pctile_frame,params.checkGuesses,...
                    params.mask_fname,params.make_guessmovie);
            else
                Guessing([dlocs{ii},filesep,dnames{ii}],mov,movsz,goodframe,dfrlmsz,...
                    params.bpthrsh,params.egdesz,params.pctile_frame,params.checkGuesses,...
                    params.mask_fname,params.make_guessmovie);
            end
        else
            Guessing([dlocs{ii},filesep,dnames{ii}],mov,movsz,goodframe,dfrlmsz,...
                params.bpthrsh,params.egdesz,params.pctile_frame,params.checkGuesses,...
                params.mask_fname,params.make_guessmovie);
        end
        clear bgsub_mov
    end   
    %% Make the off frames list
    % Only if doing bgsub. 
    % Loop through all of the movies and using the guesses .mat file, will write
    % the off frames list to a .mat file, called guessesname_Mol_off_frames.mat
    if params.makeOffFrames && params.bgsub
        try; waitbar((7*ii-3)/numel(dlocs)/7,h2,{['Making off-frames lists for ',dnames{ii}],'Overall Progress'}); end
        %try loading in the guesses
        try
            load([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],'guesses','dfrlmsz')
        catch
            warning('No avgsub_guesses file was found, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
            keyboard
            params.bgsub=0;
        end
        if params.bgsub
            Mol_off_frames([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],...
                guesses,goodframe,movsz,dfrlmsz,moloffwin);
        end
        clear guesses
    end
    %% Subtract and fit
    % If not doing bgsub then a string ('nobgsub') is sent to
    % Subtract_then_fit which will proceed accordingly.     
    % Loop through the guesses and fit the subtracted images using the off
    % frames list. If doing bgsub outputs a .mat file, called
    % moviename_AccBGSUB_fits.mat, otherwise if not doing bgsub, it's
    % called moviename_fits.mat.
    if params.fitting
        try; waitbar((7*ii-2)/numel(dlocs)/7,h2,{['Fitting ',dnames{ii}],'Overall Progress'}); end
        if params.bgsub
            %try loading in the mol_off_frames
            try
                load([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses_Mol_off_frames.mat'],'off_frames','dfrlmsz','moloffwin')
            catch
                warning('No Mol_off_frames file was found, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
                keyboard
                params.bgsub=0;
            end            
            if params.bgsub
                % load in the guesses
                load([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],'guesses');
                % fit it
                Subtract_then_fit([dlocs{ii},filesep,dnames{ii}],mov,movsz,...
                    off_frames,moloffwin,guesses,dfrlmsz,params.MLE_fit,params.stdtol,...
                    params.maxerr,params.do_avgsub,params.which_gaussian,params.fit_ang,params.usegpu);
            else
                % load in the guesses
                load([dlocs{ii},filesep,dnames{ii},'_guesses.mat'],'guesses','dfrlmsz')
                % fit it
                Subtract_then_fit([dlocs{ii},filesep,dnames{ii}],mov,movsz,...
                    'nobgsub','nobgsub',guesses,dfrlmsz,params.MLE_fit,params.stdtol,...
                    params.maxerr,params.do_avgsub,params.which_gaussian,params.fit_ang,params.usegpu);
            end
        else
            % load in the guesses
            load([dlocs{ii},filesep,dnames{ii},'_guesses.mat'],'guesses','dfrlmsz')
            % fit it
            Subtract_then_fit([dlocs{ii},filesep,dnames{ii}],mov,movsz,...
                'nobgsub','nobgsub',guesses,dfrlmsz,params.MLE_fit,params.stdtol,...
                params.maxerr,params.do_avgsub,params.which_gaussian,params.fit_ang,params.usegpu);
        end
        clear guesses off_frames
    end
    
    %% Tracking
    % Does tracking and appends the tracks to the fits .mat file. Also will
    % append a logical vector called trk_filt which indicates if the fit
    % passed was successfully tracked and wasn't the first or last frame in
    % a track.
    if params.tracking
        try; waitbar((7*ii-1)/numel(dlocs)/7,h2,{['Tracking ',dnames{ii}],'Overall Progress'}); end
        if params.bgsub
            % try loading in the fits
            try
                load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits')
            catch
                warning('Fitting file was without AccBGSUB, but bgsub was true.\r\t\t Either run new instance from beginning or continue without bgsub.',1)
                keyboard
                params.bgsub=0;
                load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits')
            end
            % track it
            if params.bgsub
                Track_filter([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],fits,...
                    1,params.trackparams,params.savetracks);
            else
                Track_filter([dlocs{ii},filesep,dnames{ii},'_fits.mat'],fits,...
                    1,params.trackparams,params.savetracks);
            end
        else            
            load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits')
            Track_filter([dlocs{ii},filesep,dnames{ii},'_fits.mat'],fits,...
                1,params.trackparams,params.savetracks);
        end
        clear fits
    end
    
    %% Make the ViewFits movie
    % Make a ViewFits movie, or just go into debug mode, to look at the
    % results. Outpits an avi file called moviename_ViewFits.avi
    if params.makeViewFits
        try; waitbar((7*ii)/numel(dlocs)/7,h2,{['Making Viewfits ',dnames{ii}],'Overall Progress'}); end
        if params.bgsub
            %try loading in the fits & tracking results
            load([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],'fits','trk_filt','tracks');
            
            %set trk_filt to empty if it doesn't exist            
            if ~exist('trk_filt','var')
                trk_filt=[];
            end
                        
            if params.orig_movie
                if ~params.trackingVF
                    ViewFits([dlocs{ii},filesep,dnames{ii},'.mat'],...
                        mov,trk_filt,movsz,goodframe,fits,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                else
                    if ~exist('bgsub_mov','var')
                        try; bgsubmov=load([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],'mov'); end
                        bgsub_mov=single(bgsubmov.mov);
                        clear bgsubmov
                    end
                    ViewFitsTracking([dlocs{ii},filesep,dnames{ii},'.mat'],...
                        mov,movsz,tracks,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                end
            else
                if ~params.trackingVF
                    ViewFits([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],...
                        bgsub_mov,trk_filt,movsz,goodframe,fits,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                else
                    ViewFitsTracking([dlocs{ii},filesep,dnames{ii},'_avgsub.mat'],...
                        bgsub_mov,movsz,tracks,...
                        params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
                end
            end
        else
            %try loading in the fits & tracking results
            load([dlocs{ii},filesep,dnames{ii},'_fits.mat'],'fits','trk_filt','tracks');
            
            %set trk_filt to empty if it doesn't exist            
            if ~exist('trk_filt','var')
                trk_filt=[];
            end
            
            if ~params.trackingVF
                ViewFits([dlocs{ii},filesep,dnames{ii},'.mat'],...
                    mov,trk_filt,movsz,goodframe,fits,...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            else
                ViewFitsTracking([dlocs{ii},filesep,dnames{ii},'.mat'],...
                    mov,movsz,tracks,...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            end
        end
        clear fits trk_filt tracks
    end
end
try
    delete(h2)
end
tictoc=toc(wholeshabang);
% disp(num2str(tictoc))

%turn warning back on
warning('on','MATLAB:load:variableNotFound')

% try to remove gpufit path
try
    rmpath(genpath('gpufit'))
end
end
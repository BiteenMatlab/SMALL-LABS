function SMALLLABS_main(directoryname,dfrlmsz,avgwin,moloffwin,varargin)
%% SMALLLABS_main
%
%%%% If you obtained this code from anywhere other than the Biteen Lab
%%%% GitHub please visit https://github.com/BiteenMatlab/SMALL-LABS to
%%%% obtain the most up-to-date version of this code.
%
% SMALLLABS_main is the wrapper function to do real background substraction of
% single molecule imaging movies so that the molecules can be fit and their
% intensity accurately measured using the SMALL-LABS algorithm.
%
% Note: that by setting bgsub=0 SMALLLABS_main will do fitting without doing
% background subtraction. Doing this obviates a lot of parameters, it
% doesn't matter what they're set to. See the User Guide for more details.
%
%
%%%% Inputs %%%%
%%% required
%   directoryname   is the name of the directory where the movies will be
%   selected, if there is an error finding the directory the program will
%   open uigetfile in the current working directory
%
%   dfrlmsz   is the size of a diffraction limited spot in pixels. It's the
%   nominal diameter, NOT the FWHM or something similar. Must be an integer!
%   For an expected diffraction limited standard deviation, std, using the
%   full width at 20% max, dfrlmsz = std*(2*sqrt(2*log(5)))
%
%   avgwin   is the width of the temporal window (in frames) to be used for
%   the average subtraction. Needs to be an odd integer
%
%   moloffwin   is the width of the temporal window (in frames) to be
%   checked to determine in which frames that molecule was off and are thus
%   safe to subtract. Needs to be an even integer
%
%%% optional (varargin)
% See the list of optional parameters in parameters section below, details
% can be found in the user guide. They are called with a name value pair as
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
% AVGSUB_tiffs 
% saveastiff 
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
%
%%%%
% Written by Benjamin P Isaacoff at the University of Michigan last update
% 1/11/17 BPI
%
%  Copyright 2016 Benjamin P Isaacoff
%
% Licensed under the Apache License, Version 2.0 (the "License"); you may
% not use this file except in compliance with the License. You may obtain a
% copy of the License at
%
%   http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
% WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
% License for the specific language governing permissions and limitations
% under the License.

%% Parameter Defaults
% You are of course welcome to change the default values, but I would
% strongly urge you to instead set them as inputs to the function using a
% name value pair (e.g., 'bpthrsh',83). Please see the user guide for
% details about the meaning and function of these parameters. The default
% parameter values are:

%%% Actions %%%
% Do background subtraction
params.bgsub = 1;
% Do guessing
params.makeGuesses = 1;
% Check guesses
params.check_guesses = 0;
% Make the off frames list
params.makeOffFrames = 1;
% Do fitting
params.fitting = 1;
% Do tracking
params.tracking = 1;
% Make the ViewFits movie
params.makeViewFits = 1;

%%% AVGSUB parameters %%%
% subtract the temporal average? otherwise use median
params.do_avg = 1;
% running average (or median) boolean
params.runningavg = 1;
% AVGSUB offset
params.offset = 1000;

%%% Guessing parameters %%%
% the threshold percentile on the banpassed movie
params.bpthrsh = 95;
% how many pixels to ignore around the edge of the frame
params.egdesz = dfrlmsz;
% compare brightnesses in each frame? otherwise use entire movie
params.pctile_frame = 0;
% use a mask for guessing? What is it's filename? Set to 1 to look for .mat
% file w/ '_PhaseMask' appended to movie name
params.mask_fname=[];
% make a movie (.avi) of the guesses to check parameters
params.make_guessmovie = 0;

%%% Fitting parameters %%%
% do MLE fitting? If not least squares fitting will be used
params.MLE_fit = 0;
% Goodfit parameters. See Subtract_mol_off_frames for the details
params.stdtol = 1.5;
params.maxerr = 2; %if you change this please also change the if statement after the next loop
% subtract the mean of off frames? If not use median
params.do_avgsub = 1;
% Which Gaussian function to fit to when using LSQ fitting? 1. symmetric,
% 2. fixed angle asymmetric, 3. free angle asymmetric
params.which_gaussian = 1;

%%% Tracking parameters %%%
% default is [0,0.01,200,1,3,3,1,0]
% save the tracks .mat file?
params.savetracks = 0;
% minimum merit
params.trackparams(1)=0.01;
% Integration time (ms)
params.trackparams(2)=200;
% gamma
params.trackparams(3)=1;
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
params.orig_movie = 1;
% diameter of the circles showing the fits
params.circ_D = dfrlmsz;
% linewidth of the circles
params.linewidth = 1;
% write a .avi movie showing the fits. If not, goes to debug mode
params.write_mov=1;
% autoscale frame by frame?
params.autoscale_on = 0;
% use the tracking viewfits instead
params.trackingVF = 0;

%% Evaluating the inputs
%This section nearly verbatim from DJR

paramsnames=fieldnames(params);

% if any parameters are included as inputs, change the parameter mentioned
if nargin>4
    for ii=1:2:nargin-5
        whichField = strcmp(paramsnames,varargin{ii});
%         if all(~whichField)
%             warning('Check spelling. Parameter change may have not occurred.')
%         end
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

%% Checking the required inputs
% Erroring out if dfrlmsz isn't an integer, because it's a "strong"
% parameter. You should really input what you mean.
if dfrlmsz~=round(dfrlmsz);error('dfrlmsz must be an integer');end

%Rounding avgwin & moloffwin and reseting them to proper values. Not
%erroring because they are "weak" parameters. Only needed if doing bgsub
if params.bgsub
    if avgwin~=(ceil(avgwin/2)*2-1)
        avgwin=ceil(avgwin/2)*2-1;
        warning(['avgwin must be an odd integer. It has been reset to avgwin = ',num2str(avgwin)])
    end
    if moloffwin~=(ceil(moloffwin/2)*2)
        moloffwin=(ceil(moloffwin/2)*2);
        warning(['moloffwin must be an even integer. It has been reset to avgwin = ',num2str(moloffwin)])
    end
end

%% Select the movies
%This section nearly verbatim from DJR

%Select movies with uigetfile. If you make an error in specifying the
%directory, it opens in the current directory.
disp('Select the movie(s)')
try
    [datalist,dataloc,findex]=uigetfile([directoryname filesep '*.tif*'],'multiselect','on');
catch
    curdir=pwd;
    [datalist,dataloc,findex]=uigetfile([curdir filesep '*.tif*'],'multiselect','on');
end
if findex==0
    fprintf('no data selected\n')
    return
end
if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);

%% The Average Subtraction
% Only if doing bgsub.
% check which of the selected movies already have the avgsub tif stack
% made, if not then make it. If you want to remake the avgsub movie, then
% delete it or movie it from the directory

% AVGSUB_tiffs will save an average subtracted .tif stack, called
% moviename_avgsub.tif and write a .txt file with the parameters, called
% moviename_avgsub_info.txt
if params.makeGuesses && params.bgsub
    bFiles=dir([dataloc,'*_avgsub.tif']);%list all the avgsub tif stacks
    for ii=1:numel(dlocs)
        %compare the chosen files with the avgsub list and choose ones
        %without an avgsub movie
        if ~any(ismember({bFiles.name},[dnames{ii},'_avgsub.tif']))&&...
                ~any(ismember({bFiles.name},dnames{ii}))
            %do the avgsub
            AVGSUB_tiffs([dlocs{ii},filesep,dnames{ii},'.tif'],...
                params.do_avg,params.runningavg,avgwin,params.offset);
        end
    end
end

%% Make Guesses
% Only if doing bgsub.
% loop through each movie and make the guesses, which will be saved in a
% .mat file called moviename_guesses.mat
if params.makeGuesses
    for ii=1:numel(dlocs)
        if params.bgsub
            Guessing([dlocs{ii},filesep,dnames{ii},'_avgsub.tif'],dfrlmsz,...
                params.bpthrsh,params.egdesz,params.pctile_frame,params.check_guesses,...
                params.mask_fname,params.make_guessmovie);
        else
            Guessing([dlocs{ii},filesep,dnames{ii},'.tif'],dfrlmsz,...
                params.bpthrsh,params.egdesz,params.pctile_frame,params.check_guesses,...
                params.mask_fname,params.make_guessmovie);
        end
    end
end

%% Make the off frames list
% loop through all of the movies and using the guesses .mat file will write
% the off frames list to a .mat file, called guessesname_Mol_off_frames.mat
if params.makeOffFrames && params.bgsub
    for ii=1:numel(dlocs)
        Mol_off_frames([dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],dfrlmsz,moloffwin);
    end
end

%% Subtract and fit
% If not doing bgsub then a string ('nobgsub') is sent to Subtract_then_fit
% which will proceed accordingly.
% loop through the movies and fit the subtracted images using the off
% frames list. If doing bgsub outputs a .mat file, called
% moviename_AccBGSUB_fits.mat, otherwise if not doing bgsub, it's called
% moviename_fits.mat.
if params.fitting
    for ii=1:numel(dlocs)
        if params.bgsub
            Subtract_then_fit([dlocs{ii},filesep,dnames{ii},'.tif'],...
                [dlocs{ii},filesep,dnames{ii},'_avgsub_guesses_Mol_off_frames.mat'],...
                [dlocs{ii},filesep,dnames{ii},'_avgsub_guesses.mat'],...
                params.MLE_fit,params.egdesz,params.stdtol,params.maxerr,...
                params.do_avgsub,params.which_gaussian);
        else
            Subtract_then_fit([dlocs{ii},filesep,dnames{ii},'.tif'],'nobgsub',...
                [dlocs{ii},filesep,dnames{ii},'_guesses.mat'],...
                params.MLE_fit,params.egdesz,params.stdtol,params.maxerr,...
                params.do_avgsub,params.which_gaussian);
        end
    end
end

%% Tracking
% loop through each movie and does tracking and appends the tracks to the
% fits .mat file. Also will append a logical vector called trk_filt which
% indicates if the fit passed was successfully tracked and wasn't the first
% or last frame in a track
if params.tracking
    for ii=1:numel(dlocs)
        if params.bgsub
            Track_filter([dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],1,params.trackparams,params.savetracks);
        else
            Track_filter([dlocs{ii},filesep,dnames{ii},'_fits.mat'],1,params.trackparams,params.savetracks);
        end
    end
end

%% Make the ViewFits movie
% loop through each movie and make a ViewFits movie, or just go into debug
% mode, to look at the results. Outpits an avi file called
% moviename_ViewFits.avi
if params.makeViewFits
    for ii=1:numel(dlocs)
        if params.bgsub
            if params.orig_movie
                VF_fname=[dlocs{ii},filesep,dnames{ii},'.tif'];
            else
                VF_fname=[dlocs{ii},filesep,dnames{ii},'_avgsub.tif'];
            end
            if ~params.trackingVF
                ViewFits(VF_fname,[dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            else
                ViewFitsTracking(VF_fname,[dlocs{ii},filesep,dnames{ii},'_AccBGSUB_fits.mat'],...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            end
        else
            if ~params.trackingVF
                ViewFits([dlocs{ii},filesep,dnames{ii},'.tif'],...
                    [dlocs{ii},filesep,dnames{ii},'_fits.mat'],...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            else
                ViewFitsTracking([dlocs{ii},filesep,dnames{ii},'.tif'],...
                    [dlocs{ii},filesep,dnames{ii},'_fits.mat'],...
                    params.circ_D,params.write_mov,params.autoscale_on,params.linewidth)
            end
        end
    end
end

end































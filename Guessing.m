function guesses=Guessing(mov_fname,dfrlmsz,bpthrsh,egdesz,pctile_frame,debugmode,mask_fname,make_guessmovie)
%% Guessing
% make a list of guesses for sinlge molecules. Using a bandpass filter to
% filter pixel noise first, then uses bwpropfilt to find blobs of the
% correct size
%
%%%% Inputs %%%%
% mov_fname is the full filename of the tif stack movie to be analyzed
%
% dfrlmsz is the  size of a diffraction limited spot in pixels. It's the
% nominal diameter, NOT the FWHM or something similar. Integer please!
%
% bpthrsh is the the percentile of brightnesses of the bandpassed image
% below which those pixels will be ignored. default value is 90
%
% edgesz is the number of pixels on the edge of the image that will be
% ignored. Default is egdesz = dfrlmsz
%
% pctile_frame is a boolean determining whether bpthrsh will be applied
% frame by frame, or to the entire movie. Using the entire movie (setting
% to 0) is more sensitive to low frequency noise and background changes,
% but is a more robust guessing method. Using each frame tends to produce a
% constant number of guesses per frame, regardless of their absolute
% brightness.
%
% debugmode is a boolean to determine if you want to go through and look at
% the guesses. Default is 0
%
% mask_fname is the filename of a mask to use for guessing. If no mask is
% being used just leave it empty. If mask_fname is set 1, then the program
% will look for a file in the same directory as the movie with '_PhaseMask'
% appened to the name of the movie. The mask is a .mat file which has a
% logical array (or at least where nonzero entries will be converted to 1s)
% called PhaseMask that is the same size as a frame in the current movie.
%
%%%% Outputs %%%%
% guesses is an array with columns 1. frame #, 2. row #, 3. column # of the
% guesses
%
% The program currently also writes a .mat file with guesses and all of the
% user parameters saved.
%
%
%%%% Dependencies %%%%
% bpass
% TIFFStack
%
%  Copyright 2016 Benjamin P Isaacoff
%
% Licensed under the Apache License, Version 2.0 (the "License"); you
% may not use this file except in compliance with the License. You may
% obtain a copy of the License at
%
%   http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
% implied. See the License for the specific language governing
% permissions and limitations under the License.
%
% default values
if nargin<3;bpthrsh=90;end
if nargin<4;egdesz=dfrlmsz;end
if nargin<5;pctile_frame=1;end
if nargin<6;debugmode=0;end

%did you not set dfrlmsz to an integer?
if dfrlmsz~=round(dfrlmsz);error('dfrlmsz must be an integer');end

%pad size for the bandpass function
pdsz=50;

tic;%for measuring the time to run the entire program
% last updated 8/12/16 BPI
%% Peak Guessing
tfstk=TIFFStack(mov_fname);
movsz=size(tfstk);%the size of the movie
[pathstr,fname] = fileparts(mov_fname);

%intializing the guess indices cell array
guesses=zeros(1,3);

%making the phasemask logical map
if ~isempty(mask_fname)
    if mask_fname
        %the strrep is to get rid of the avgsub, note that this shouldn't
        %do anything if bgsub=0
        load([pathstr,filesep,strrep(fname,'_avgsub',[]),'_PhaseMask'],'PhaseMask')
    else
        load(mask_fname,'PhaseMask')
    end
    PhaseMask(PhaseMask~=0)=1;
    PhaseMask=logical(PhaseMask);
else
    PhaseMask=true(movsz([1,2]));
end

%using the percentiles on the entire movie
if ~pctile_frame
    %initializing the bandpassed movie
    bimgmov=zeros(movsz);
    %looping through and making the bandpassed movie
    for ll=1:movsz(3)
        %padding the current frame to avoid the Fourier ringing associated
        %with the edges of the image
        curfrm=padarray(double(tfstk(:,:,ll)),[pdsz,pdsz],'symmetric');
        %bandpass parameters
        LP=1;%lnoise, should always be 1
        HP=round(dfrlmsz*1.5);%lobject, set by diffraction limit
        T=0;%threshold, now always zero
        lzero=egdesz;%how many pixels around the edge should be ignored, optional
        %bandpass it
        bimg=bpass(curfrm,LP,HP,T,lzero+pdsz);
        %removed the padded pixels around the edge
        bimgmov(:,:,ll)=bimg((pdsz+1):(movsz(1)+pdsz),(pdsz+1):(movsz(2)+pdsz));
    end
    
    %convert it to a logical movie by thresholding with the bpthrsh
    %percentile of the brightnesses for nonzero pixels
    bimgmov=logical(bimgmov.*(bimgmov>prctile(bimgmov(bimgmov>0 & ...
        repmat(PhaseMask,[1,1,movsz(3)])),bpthrsh)).*repmat(PhaseMask,[1,1,movsz(3)]));
end

if make_guessmovie
    v = VideoWriter([pathstr,filesep,fname,'_Guesses.avi'],'Uncompressed AVI');
    open(v);
    
    disp(['Making guesses movie for ',fname]);
end

h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,['Making guesses for ',fname]);

for ll=1:movsz(3)
    try;waitbar(ll/movsz(3),h1);end
    %using the percentile on each frame
    if pctile_frame
        %padding the current frame to avoid the Fourier ringing associated
        %with the edges of the image
        curfrm=double(tfstk(:,:,ll));
        curfrmbp=padarray(curfrm,[pdsz,pdsz],'symmetric');
        
        %bandpass parameters
        LP=1;%lnoise, should always be 1
        HP=round(dfrlmsz*1.5);%lobject, set by diffraction limit
        T=0;%threshold, now always zero
        lzero=egdesz;%how many pixels around the edge should be ignored, optional
        %bandpass it
        bimg=bpass(curfrmbp,LP,HP,T,lzero+pdsz);
        %pull out the actual data
        bimg=bimg((pdsz+1):(movsz(1)+pdsz),(pdsz+1):(movsz(2)+pdsz));
        
        %threshold with the bpthrsh percentile of the brightnesses for nonzero
        %pixels, then turn it into a logical array
        logim=logical(bimg.*(bimg>prctile(bimg(bimg>0 & PhaseMask),bpthrsh)).*PhaseMask);
    else
        logim=bimgmov(:,:,ll);
    end
    
    %search for shapes with an EquivDiameter of floor(dfrlmsz/2) to 2*dfrlmsz
    bw2=bwpropfilt(logim,'EquivDiameter',[floor(dfrlmsz/2),2*dfrlmsz]);
    rgps=regionprops(bw2,'centroid');% find the centroids of those shapes
    centroids = cat(1, rgps.Centroid);%just rearraging the array
    %filling the array for this frame
    if ~isempty(centroids)
        guesses=cat(1,guesses,[repmat(ll,size(centroids(:,2))),round(centroids(:,2)),round(centroids(:,1))]);
    end
    
    if debugmode || make_guessmovie %plot the guesses, for checking parameters
        if ~pctile_frame
            curfrm=double(tfstk(:,:,ll));
        end
        imshow(curfrm,prctile(curfrm(curfrm>0),[.1,99.8]))
        if ~isempty(centroids)
            %viscircles is reversed
            vcs=viscircles([centroids(:,1),centroids(:,2)],repmat(dfrlmsz,[length(centroids(:,2)),1]));
            set(vcs.Children,'LineWidth',1)
        end
        if debugmode
            title([fname,'   frame ',num2str(ll)],'Interpreter','none')
            keyboard
        elseif make_guessmovie
            frame = getframe;
            writeVideo(v,frame);
        end
    end
end
guesses=guesses(2:end,:);%get rid of first row of zeros

if make_guessmovie
    close(v)
    keyboard
end

tictoc=toc;%the time to run the entire program

[pathstr,name,~] = fileparts(mov_fname);
save([pathstr,filesep,name,'_guesses.mat'],'guesses','dfrlmsz','egdesz','pctile_frame','bpthrsh',...
    'movsz','tictoc','mask_fname');

try
    close(h1)%closing the waitbar
end

end











function off_frames=Mol_off_frames(guessfname,dfrlmsz,moloffwin)
%% Mol_off_frames
% Identifies frames in which two localized molecules are within a 2xdfrlmsz
% sized box of each other using guess results from an average subtracted
% movie
%
%%%% Inputs %%%%
% guessfname is the name of the .mat file with the guesses
%
% dfrlmsz is the  size of a diffraction limited spot in pixels. It's the
% nominal diameter, NOT the FWHM or something similar. Integer please!
%
% moloffwin is the number of frames around the current frame to use for the BGSUB. For
% example if moloffwin=50 then you would subtract 25 frames before and after the
% current frame. Even number please!
%
%%%% Outputs %%%%
% off_frames is a list of off frames around each guess
%
% The program currently also writes a .mat file with off list and all of the
% parameters saved.

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

% Erroring out if dfrlmsz isn't an integer, because it's a "strong"
% parameter. You should really input what you mean.
if dfrlmsz~=round(dfrlmsz);error('dfrlmsz must be an integer');end
%rounding moloffwin
if moloffwin~=(ceil(moloffwin/2)*2)
    moloffwin=(ceil(moloffwin/2)*2);
    warning(['moloffwin must be an even integer. It has been reset to avgwin = ',num2str(moloffwin)])
end

% last updated 8/12/16 BPI

% NOTE: if you're inspecting this code, you'll notice that I'm making use
% of ismembc instead of ismember. ismembc is an undocumented helper
% function that is much faster than ismember in this case. Basically
% ismember eventually calls ismembc, but wastes a lot of time before
% getting there that we've cut out.

tic;%for measuring the time to run the entire program
%% Constructing off frames lists

[pathstr,fname,~] = fileparts(guessfname);
%load in the guesses & the movie size
load(guessfname,'guesses','movsz');

%cell array of off frames vectors, for each localization (each row of
%guesses) a list of frames to include
off_frames=cell(size(guesses,1),1);

h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,['Creating off frames list for ',fname]);
for ii=1:movsz(3)
    try; waitbar(ii/movsz(3),h1); end
    %number of molecules in the current frame
    frmrows=find(guesses(:,1)==ii);
    nummol=length(frmrows);
    
    %skip it if there aren't any guesses in the current frame
    if nummol~=0
        %determine the frame list of frames to check for the current frame
        if ii<=(moloffwin/2)%the first group of frames
            frmlst=ii+[-(ii-1):-1,1:(moloffwin/2)];
        elseif ii>=(movsz(3)-moloffwin/2)%the last group of frames
            frmlst=movsz(3)+(-moloffwin:0);
        else %all the frames in the middle
            frmlst=ii+[(-moloffwin/2):-1,1:(moloffwin/2)];
        end
        
        %find the rows of the guesses that are in the current frame list,
        %then pull out the row & column #'s
        mols2frms=find(ismembc(guesses(:,1),frmlst));
        allr=guesses(mols2frms,2);
        allc=guesses(mols2frms,3);        
        
        for jj=frmrows(1):frmrows(end)
            %current molecule's position
            molr=guesses(jj,2);
            molc=guesses(jj,3);
            dists=(abs(allr-molr)<dfrlmsz) & (abs(allc-molc)<dfrlmsz);
            
            %save the lists of frames  in which the current localization is off.
            off_frames{jj}=frmlst(~ismembc(frmlst,guesses((mols2frms(dists,1)),1)));
            
            % NOTE : when I first wrote this I thought that there needed to
            % be a unique as shown below. In reviewing this though, I don't
            % think it's needed and there is a big speed boost from leaving
            % it out. Please let me know if you find that this is needed.
            %off_frames{jj}=frmlst(~ismembc(frmlst,guesses(unique(mols2frms(dists,1)),1)));
            
            if length(off_frames{jj})<(0.05*moloffwin)
                warning(['Guess ',num2str(jj),' only has ',...
                    num2str(length(off_frames{jj})),' off frames. Consider increasing moloffwin'])
            end
        end        
    end
end

try;close(h1);end

tictoc=toc;%the time to run the entire program

%save the data
save([pathstr,filesep,fname,'_Mol_off_frames.mat'],'off_frames','dfrlmsz','moloffwin','movsz','tictoc')

end
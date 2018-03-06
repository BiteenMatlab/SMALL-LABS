function off_frames=Mol_off_frames(guessfname,guesses,goodframe,movsz,dfrlmsz,moloffwin)
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

%     Copyright (C) 2017  Benjamin P Isaacoff
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

% Erroring out if dfrlmsz isn't an integer, because it's a "strong"
% parameter. You should really input what you mean.
if dfrlmsz~=round(dfrlmsz);error('dfrlmsz must be an integer');end

%setting moloffwin to movie length - 1, which is the maximum meaningful
%length for moloffwin
if moloffwin>=movsz(3)
    moloffwin=movsz(3)-1;
    %round to nearest even integer
    if moloffwin~=(floor(moloffwin/2)*2)
        moloffwin=(floor(moloffwin/2)*2);
    end
    warning(['moloffwin was >=',num2str(movsz(3)),' (the number of frames in the movie). ',...
        'It has been reset to ',num2str(moloffwin)])
end
%rounding moloffwin
if moloffwin~=(floor(moloffwin/2)*2)
    moloffwin=(floor(moloffwin/2)*2);
    warning(['moloffwin must be an even integer. It has been reset to avgwin = ',num2str(moloffwin)])
end

% last updated 3/6/16 BPI & SAL

% NOTE: if you're inspecting this code, you'll notice that I'm making use
% of ismembc instead of ismember. ismembc is an undocumented helper
% function that is much faster than ismember in this case. Basically
% ismember eventually calls ismembc, but wastes a lot of time before
% getting there that we've cut out.

tic;%for measuring the time to run the entire program
%% Constructing off frames lists

[pathstr,fname,~] = fileparts(guessfname);

disp([char(datetime),'   Making off-frames lists for ',fname])

%load in the guesses & the movie size
%cell array of off frames vectors, for each localization (each row of
%guesses) a list of frames to include
off_frames=cell(size(guesses,1),1);

for ii=1:movsz(3)
    
    %number of molecules in the current frame
    frmrows=find(guesses(:,1)==ii);
    nummol=length(frmrows);
    
    %skip it if there aren't any guesses in the current frame
    if nummol~=0
        if ii<=(moloffwin/2)%the first group of frames
            frmlst=1:(moloffwin+1);
            frmlst=frmlst(frmlst~=ii);
        elseif ii>=(movsz(3)-moloffwin/2)%the last group of frames
            frmlst=movsz(3)+((-moloffwin):0);
            frmlst=frmlst(frmlst~=ii);
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
            
            %save the lists of frames  in which the current localization is
            %off and it is safe to subtract (using goodframes)
            off_frames{jj}=frmlst(~ismembc(frmlst,guesses((mols2frms(dists,1)),1)) & goodframe(frmlst)');
            
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


tictoc=toc;%the time to run the entire program

%save the data
save([pathstr,filesep,fname,'_Mol_off_frames.mat'],'off_frames','dfrlmsz','moloffwin','movsz','tictoc')

end
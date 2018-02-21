function  bgsub_mov=AVGSUB_movs(filename,mov,goodframe,do_avg,subwidth,offset)
%% AVGSUB_tiffs
% updated BPI 1/30/17 This function does average (or median) subtraction
% on a batch of movies. Currently does the entire movie, I might add
% functionality to do a portion of the frames only.

% filename is the name of the movie to be subtracted, either .mat or tiff
% stack

% do_avg is a boolean. Set to 1 to use an average, or set to 0 to use a
% median. Note median is a lot slower than mean.

% runningavg is a boolean. Set to 1 to do a running average or set to 0 to
% do a static window background subtraction

% subwidth is the width of temporal window that will be averaged. Needs to
% be an odd integer

% offset is the intensity offset to deal with negative pixels. Default is
% 1000

% This functions doesn't output anything, instead saves the AVGSUB movie as
% a .mat file at the original file location with the original file name,
% but with _avgsub appended to the name. Also outputs a short text file
% with the parameters used


%%%% Dependencies %%%%

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

%default offset
tic;%for measuring the time to run the entire program

%% Setup
[pathstr,fname] = fileparts(filename);

    disp(['Running AVGSUB_movs for ',fname])

%change the bad frames as determined by the goodframe list to NaNs
mov(:,:,~goodframe)=NaN;

%% Do the avgsub
if do_avg
    bgsub_mov=mov-movmean(mov,subwidth,3,'omitnan')+offset;
else
    bgsub_mov=mov-movmedian(mov,subwidth,3,'omitnan')+offset;
end

% bgsub_mov(:,:,~gfs)=NaN;

%% Check and save

if any(bgsub_mov(:)<0)
    warning(['There are negative pixels in ',fname,'_avgsub, consider increasing the offset'])
end

tictoc=toc;%the time to do the calculations

%%%%Save the movie%%%%
%rename the bgsub movie for saving
mov=int16(bgsub_mov);
save([pathstr,filesep,fname,'_avgsub.mat'],'mov','goodframe','do_avg','subwidth','offset','tictoc','-v7.3')

%this is the old tif saving code, commented out for now
% options.overwrite=true;
% %save using saveastiff
% saveastiff(uint16(bgsub_mov), [pathstr,filesep,fname,'_avgsub.tif'],options);

% %save a .txt file with the parameters
% fileID = fopen([pathstr,filesep,fname,'_avgsub_info.txt'],'w');
% if do_avg
%     fprintf(fileID,'Average of background was subtracted\n');
% else
%     fprintf(fileID,'Median of background was subtracted\n');
% end
% fprintf(fileID,['running average:\t',num2str(runningavg),'\n']);
% fprintf(fileID,['subwidth:\t',num2str(subwidth),'\n']);
% fprintf(fileID,['offset:\t',num2str(offset),'\n']);
% fprintf(fileID,['time to run:\t',num2str(tictoc)]);
% fclose(fileID);

end


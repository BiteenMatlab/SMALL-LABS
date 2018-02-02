function  AVGSUB_movs(filename,do_avg,runningavg,subwidth,offset)
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

%rounding subwidth up to the next odd number down, in case you didn't read the
%instructions
if subwidth~=(ceil(subwidth/2)*2-1)
    subwidth=ceil(subwidth/2)*2-1;
    warning(['subwidth must be an odd integer. It has been reset to subwidth = ',num2str(subwidth)])
end
%default offset
if nargin<4;offset=1000;end

tic;%for measuring the time to run the entire program

%% Setup

% check if a GPU is available
try
    usegpu=parallel.gpu.GPUDevice.isAvailable;
catch
    usegpu=false;
end

[pathstr,fname,ext] = fileparts(filename);

if ~strcmp(ext,'.mat')
    error('Please use a .mat with the variable ''mov''')
end

matio=matfile(filename,'Writable',false);
%get the movie size
movsz=whos(matio,'mov');
movsz=movsz.size;

%look for a goodframe list, otherwise set all frames as goodframes
try
    goodframe=matio.goodframe;
catch
    goodframe=true(movsz(3),1);
end
 
if usegpu
    disp(['Running AVGSUB_movs on a GPU for ',fname])
    mov=gpuArray(int16(matio.mov));
else
    disp(['Running AVGSUB_movs for ',fname])
    mov=int16(matio.mov);
end
%change the bad frames as determined by the goodframe list to NaNs
mov(:,:,~goodframe)=NaN;

%% Do the avgsub
if do_avg
    bgsub_mov=mov-int16(movmean(mov,subwidth,3,'omitnan'))+offset;
else
    bgsub_mov=mov-int16(movmedian(mov,subwidth,3,'omitnan'))+offset;
end

% bgsub_mov(:,:,~gfs)=NaN;

%% Check and save

if any(bgsub_mov(:)<0)
    warning(['There are negative pixels in ',fname,'_avgsub, consider increasing the offset'])
end

tictoc=toc;%the time to do the calculations

%%%%Save the movie%%%%
%rename the bgsub movie for saving
if usegpu
    mov=int16(gather(bgsub_mov));
else
    mov=int16(bgsub_mov);
end
save([pathstr,filesep,fname,'_avgsub.mat'],'mov','goodframe','runningavg','subwidth','offset','tictoc','-v7.3')

%this is the old tif saving code, commented out for now
% options.overwrite=true;
% %save using saveastiff
% saveastiff(uint16(bgsub_mov), [pathstr,filesep,fname,'_avgsub.tif'],options);

%close the waitbar
try; close(h1); end

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


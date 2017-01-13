function  AVGSUB_tiffs(filename,do_avg,runningavg,subwidth,offset)
%% AVGSUB_tiffs
% updated BPI 10/13/16
% This function does average (or median) subtraction on a batch of tiff stacks
% movies. Currently does the entire movie, I might add functionality to do
% a portion of the frames only.

% filename is the name of the tiffstack to be subtracted

% do_avg is a boolean. Set to 1 to use an average, or set to 0 to use a
% median. Note median is a lot slower than mean.

% runningavg is a boolean. Set to 1 to do a running average or set to 0 to
% do a static window background subtraction

% subwidth is the width of temporal window that will be averaged. Needs to
% be an odd integer

% offset is the intensity offset to deal with negative pixels. Default is
% 1000

% This functions doesn't output anything, instead saves the AVGSUB movie as
% a .tif stack at the original file location with the original file name,
% but with _avgsub appended to the name. Also outputs a short text file
% with the parameters used


%%%% Dependencies %%%%
% TIFFStack
% saveastiff

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

%rounding subwidth up to the next odd number down, in case you didn't read the
%instructions
if subwidth~=(ceil(subwidth/2)*2-1)
    subwidth=ceil(subwidth/2)*2-1;
    warning(['subwidth must be an odd integer. It has been reset to subwidth = ',num2str(subwidth)])
end

tic;%for measuring the time to run the entire program

%% Do the AVGSUB

[pathstr,fname,~] = fileparts(filename);

%default offset
if nargin<4;offset=1000;end

%create A `TIFFStack` object  which behaves like a read-only memory
%mapped TIFF file
tfstk=TIFFStack(filename);
movsz=size(tfstk);%the size of the movie

%create the waitbar
h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,['writing avg sub frames for ',fname]);

%initialize the bgsub movie
bgsub_mov=zeros(movsz);
%different loops depending on running or static window
if runningavg
    if mod(subwidth,2)==0;error('subwidth needs to be an odd integer');end
    v=tfstk(:,:,1:subwidth); %the first subwidth frames
    for jj=1:movsz(3)
        try;waitbar(jj/movsz(3),h1);end
        curr_v=tfstk(:,:,jj);%current frame
        if jj>(movsz(3)-floor(subwidth/2))%the last subwidth frames
            if jj==(movsz(3)-floor(subwidth/2))+1%only need to fill it once
                v=tfstk(:,:,(movsz(3)-subwidth+1):movsz(3));
            end
        elseif jj>floor(subwidth/2)%frames in the middle
            v=cat(3,v(:,:,2:end),tfstk(:,:,jj+floor(subwidth/2)));
        end
        if do_avg
            bgsub_mov(:,:,jj)=bsxfun(@plus,double(curr_v),-mean(v,3))+offset;
        else
            bgsub_mov(:,:,jj)=bsxfun(@plus,double(curr_v),-double(median(v,3)))+offset;
        end
    end
else
    for jj=1:floor(movsz(3)/subwidth)
        try;waitbar(jj/movsz(3),h1);end
        if jj~=floor(movsz(3)/subwidth)
            v=tfstk(:,:,1+subwidth*(jj-1):subwidth*jj);
            if do_avg
                bgsub_mov(:,:,1+subwidth*(jj-1):subwidth*jj)=bsxfun(@plus,double(v),-mean(v,3))+offset;
            else
                bgsub_mov(:,:,1+subwidth*(jj-1):subwidth*jj)=bsxfun(@plus,double(v),-double(median(v,3)))+offset;
            end
        else
            v=tfstk(:,:,1+subwidth*(jj-1):movsz(3));
            if do_avg
                bgsub_mov(:,:,1+subwidth*(jj-1):movsz(3))=bsxfun(@plus,double(v),-mean(v,3))+offset;
            else
                bgsub_mov(:,:,1+subwidth*(jj-1):movsz(3))=bsxfun(@plus,double(v),-double(median(v,3)))+offset;
            end
        end
    end
end

if any(bgsub_mov(:)<0)
    warning(['There are negative pixels in ',fname,'_avgsub, consider increasing the offset'])
end

%%%%Save the movie%%%%
options.overwrite=true;
%save using saveastiff
saveastiff(uint16(bgsub_mov), [pathstr,filesep,fname,'_avgsub.tif'],options);

%close the waitbar
try; close(h1); end

tictoc=toc;%the time to run the entire program

%save a .txt file with the parameters
fileID = fopen([pathstr,filesep,fname,'_avgsub_info.txt'],'w');
if do_avg
    fprintf(fileID,'Average of background was subtracted\n');
else
    fprintf(fileID,'Median of background was subtracted\n');
end
fprintf(fileID,['running average:\t',num2str(runningavg),'\n']);
fprintf(fileID,['subwidth:\t',num2str(subwidth),'\n']);
fprintf(fileID,['offset:\t',num2str(offset),'\n']);
fprintf(fileID,['time to run:\t',num2str(tictoc)]);
fclose(fileID);

end


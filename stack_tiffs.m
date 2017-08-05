function  stack_tiffs(directoryname,stack_name)
%% stack_tiffs
% written BPI 10/20/16
% is a quick function that will convert a folder of individual tiff images
% into a tiffstack.
%
% directoryname is the name of the directory holding the tiff images.
%
% stack_name is the name of the tiff stack that will be saved. If no input
% is given the default is directorynamestack.tif and it will be saved in
% the directory containing the directoryname. This is also the save
% directory if no directory is specified in stack_name.
%
% This program is very simple, so use it cautiously. It assumes that the ls
% command lists the files in the correct order. It takes all of the tiffs
% in a folder and saves them as a tiff stack.
%
%%%% Dependencies %%%%
% TIFFStack
% saveastiffs

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
%%
%move to directory with the imagess
curdir=pwd;
cd(directoryname);

%setting up the save filename
if nargin<2
    save_name=[directoryname,'stack.tif'];
else
    [save_directory,snam,~]=fileparts(stack_name);
    if isempty(save_directory)
        cd ..
        save_directory=pwd;
        cd(directoryname);
    end
    save_name=[save_directory,filesep,snam,'.tif'];
end

%listing all of the tif images
fnames=cellstr(ls('*.tif'));

%the first frame just to get the size
frame=TIFFStack(fnames{1});
%initializing the movie
mov=zeros([size(frame,1),size(frame,2),length(fnames)]);

%loop through each frame and save it to mov
for ii=1:length(fnames)
    frame=TIFFStack(fnames{ii});
    mov(:,:,ii)=frame(:,:);
end

%save it
options.overwrite=true;
saveastiff(uint16(mov),save_name,options);

%return to original directory
cd(curdir);
end
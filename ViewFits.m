function  ViewFits(movfname,fits_fname,circ_D,write_mov,autoscale_on,linewidth)
%ViewFits plots or writes a view fits movie using the tiff stack movie
%specified by mov_fname, and the fits from Subtract_mol_off_frames mat file
%specficied by fits_fname.
%
% circ_D is the diameter of the circle to show the fit, default is 7
%
% write_mov is a boolean determining whether the viewfits movie will be
% written to an avi. If set to 0 this function will be used in debug mode
%
% autoscale_on is a boolean determining if the movie grayscale will
% be set frame by frame. If set to 0 a handful of frames throughout the movie
% are used to set the grayscale
% 
% linewidth is the linewidth of the circles in the movie. Default is 1
%
%%%% Color Scheme %%%%
% green circles are good fits
% red circles are bad fits 
% magenta circles are fits which passed the tracking filter

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

if nargin<3;circ_D=7;end
if nargin<4;write_mov=0;end
if nargin<5;autoscale_on=0;end
if nargin<6;linewidth=1;end
%% setup

%load in fits
load(fits_fname);

matio=matfile(movfname,'Writable',false);
%get the movie size
movsz=whos(matio,'mov');
movsz=movsz.size;

[pathstr,name,ext] = fileparts(movfname);

%look for a goodframe list, otherwise set all frames as goodframes
try
    goodframe=matio.goodframe;
catch
    goodframe=true(movsz(3),1);
end

if write_mov
    v = VideoWriter([pathstr,filesep,name,'_ViewFits.avi'],'Uncompressed AVI');
    open(v);
    
    disp(['Making ViewFits for ',name]);
end

%load in the movie
mov=double(matio.mov);

%% Make the ViewFits movie

if ~autoscale_on
    frms4scl=100 ;%number of frames for the scaling
    %pull out some frames to find the percentiles of
    frmscl=mov(:,:,round(linspace(1,movsz(3),frms4scl)));
    % the intensity bounds for not autoscaling
    int_bounds=prctile(frmscl(frmscl>0),[.1,99.8]);
end

figure
set(gcf,'Position',[1281,1,1280,948])
for ii=1:movsz(3)    
    
    curfrm=mov(:,:,ii);
    
    if autoscale_on
        int_bounds=prctile(curfrm(curfrm>0),[.1,99.8]);
    end
    imshow(curfrm,int_bounds)
    
    %which molecules appear in this frame
    thisfrm_gf=[fits.row(fits.frame==ii & fits.goodfit),...
        fits.col(fits.frame==ii & fits.goodfit)];    
    thisfrm_bf=[fits.row(fits.frame==ii & ~fits.goodfit),...
        fits.col(fits.frame==ii & ~fits.goodfit)];
    
    if exist('trk_filt','var')
        thisfrm_trk=[fits.row(fits.frame==ii & trk_filt),...
        fits.col(fits.frame==ii & trk_filt)];
    else
        thisfrm_trk=[];
    end 
    
    %note that for viscircles row and column are switched...
    if ~isempty(thisfrm_gf)
        vcs=viscircles([thisfrm_gf(:,2),thisfrm_gf(:,1)],repmat(circ_D,[length(thisfrm_gf(:,1)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color','green')
    end
    if ~isempty(thisfrm_bf)
        vcs=viscircles([thisfrm_bf(:,2),thisfrm_bf(:,1)],repmat(circ_D,[length(thisfrm_bf(:,1)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color','red')
    end
    if ~isempty(thisfrm_trk)
        vcs=viscircles([thisfrm_trk(:,2),thisfrm_trk(:,1)],repmat(circ_D+4,[length(thisfrm_trk(:,1)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color','magenta')
    end
    
    %check for a bad frame
    if ~goodframe(ii)
      rectangle('Position',[2,2,size(curfrm,2)-2,size(curfrm,1)-2],'EdgeColor','red',...
          'LineWidth',3)  
    end
    
    if write_mov
        frame = getframe;
        writeVideo(v,frame);
    else
        title([name,' frame #',num2str(ii)],'Interpreter', 'none')
        keyboard
    end    
end
if write_mov
    close(v)
end
end


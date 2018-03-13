function  ViewFits(movfname,mov,trk_filt,movsz,goodframe,fits,circ_D,write_mov,...
    autoscale_on,linewidth)
%ViewFits plots or writes a view fits movie
% 
%%%% Inputs %%%%
%
% movfname is the filename of the movie, for output naming purposes
%
% mov is the movie data as a 3D array where the third dimension is the
% frame number.
% 
% trk_filt is the track filter from Track_filter.m
% 
% movsz is the output of size(mov)
%
% goodframe is an optional logical vector indicating which frames are to be
% ignored. The length of the goodframe vector should be the number of
% frames in mov. To ignore a frame set the corresponding element in
% goodframe to false.
%
% fits is the fit structure outputted from Subtract_then_fit.m
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
% linewidth is the linewidth of the circles in the movie.
%
%%%% Color Scheme %%%%
% green circles are good fits
% red circles are bad fits 
% magenta circles are fits which passed the tracking filter

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
%% setup

[pathstr,name] = fileparts(movfname);
global verbose
if verbose
disp([char(datetime),'   Making ViewFits movie for ',name])
end
%look for a goodframe list, otherwise set all frames as goodframes
if write_mov
    v = VideoWriter([pathstr,filesep,name,'_ViewFits.avi'],'Uncompressed AVI');
    open(v);
end


%% Make the ViewFits movie

if ~autoscale_on
    frms4scl=100 ;%number of frames for the scaling
    %pull out some frames to find the percentiles of
    frmscl=double(mov(:,:,round(linspace(1,movsz(3),frms4scl))));
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
    close(gcf)
end
end


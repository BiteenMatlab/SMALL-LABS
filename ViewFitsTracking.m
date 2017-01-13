function  ViewFitsTracking(movfname,fits_fname,circ_D,write_mov,autoscale_on,linewidth)
%ViewFitstTracking plots or writes a view fits movie using the tiff stack movie
%specified by mov_fname, and the fits from Subtract_mol_off_frames mat file
%specficied by fits_fname.
%
%Current implementation uses 5 colors to identify different tracks. The
%program then loops through those 5 colors (e.g., track 6 is the same color
%as track 1). Only successfully tracked molecules are shown and identified.
%If you want to change the colors you'll need to modify the code, which
%I would rate as a fairly easy task.
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
if nargin<3;circ_D=7;end
if nargin<4;write_mov=0;end
if nargin<5;autoscale_on=0;end
if nargin<6;linewidth=1;end
%%

%load in fits
load(fits_fname);

%create A `TIFFStack` object  which behaves like a read-only memory
%mapped TIFF file
tfstk=TIFFStack(movfname);
movsz=size(tfstk);%the size of the movie
[pathstr,name,~] = fileparts(movfname);

if write_mov
    v = VideoWriter([pathstr,filesep,name,'_ViewFitsTracking.avi'],'Uncompressed AVI');
    open(v);
    
    disp(['Making ViewFits for ',name]);
end

if ~autoscale_on
    frms4scl=100 ;%number of frames for the scaling
    %pull out some frames to find the percentiles of
    frmscl=double(tfstk(:,:,round(linspace(1,movsz(3),frms4scl))));
    % the intensity bounds for not autoscaling
    int_bounds=prctile(frmscl(frmscl>0),[.1,99.8]);
end

%colors to use for the tracks
colors=jet(5);

figure
set(gcf,'Position',[1281,1,1280,948])
for ii=1:movsz(3)
    
    curfrm=double(tfstk(:,:,ii));
    
    if autoscale_on
        int_bounds=prctile(curfrm(curfrm>0),[.1,99.8]);
    end
    imshow(curfrm,int_bounds)
    
    %which molecules appear in this frame and which track do they belong to
    thisfrm_c1=tracks(tracks(:,1)==ii & mod(tracks(:,4),5)==1,:);%goodfits
    thisfrm_c2=tracks(tracks(:,1)==ii & mod(tracks(:,4),5)==2,:);%goodfits
    thisfrm_c3=tracks(tracks(:,1)==ii & mod(tracks(:,4),5)==3,:);%goodfits
    thisfrm_c4=tracks(tracks(:,1)==ii & mod(tracks(:,4),5)==4,:);%goodfits
    thisfrm_c5=tracks(tracks(:,1)==ii & mod(tracks(:,4),5)==0,:);%goodfits
    
    %note that for plotting row and column are switched...
    if ~isempty(thisfrm_c1)
        vcs=viscircles([thisfrm_c1(:,3),thisfrm_c1(:,2)],repmat(circ_D,[length(thisfrm_c1(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color',colors(1,:))
    end
    if ~isempty(thisfrm_c2)
        vcs=viscircles([thisfrm_c2(:,3),thisfrm_c2(:,2)],repmat(circ_D,[length(thisfrm_c2(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color',colors(2,:))
    end
    if ~isempty(thisfrm_c3)
        vcs=viscircles([thisfrm_c3(:,3),thisfrm_c3(:,2)],repmat(circ_D,[length(thisfrm_c3(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color',colors(3,:))
    end
    if ~isempty(thisfrm_c4)
        vcs=viscircles([thisfrm_c4(:,3),thisfrm_c4(:,2)],repmat(circ_D,[length(thisfrm_c4(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color',colors(4,:))
    end
    if ~isempty(thisfrm_c5)
        vcs=viscircles([thisfrm_c5(:,3),thisfrm_c5(:,2)],repmat(circ_D,[length(thisfrm_c5(:,3)),1]));
        set(vcs.Children,'LineWidth',linewidth,'Color',colors(5,:))
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


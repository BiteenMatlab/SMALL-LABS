function tracks = Tracking(fits_fname,trackparams,savetracks)
%% Tracking
% written BPI 6/7/16
% This function is just a wrapper for Track_3D2 to interface with the fit
% fits used in the SMALL-LABS fitting protocol

%%%% Inputs %%%%
% fits_fname is the name of the fits .mat file, importantly containing an
% array called fits with the fit information with column 9 being the
% goodfit boolean

% trackparams are the tracking parameters, definitions and defaults below

%%%% Outputs %%%%
% tracks is an array with the track information.

%%%% Dependencies %%%%
% Track_3D2
% hungarian

%% Default tracking parameters
if nargin<2
    % TRACKING PARAMETERS
    % minimum merit
    trackparams(1)=0.01;
    % Integration time (ms)
    trackparams(2)=200;
    % gamma
    trackparams(3)=1;
    % maximum step size
    trackparams(4)=3;
    % minimum track length
    trackparams(5)=3;
    % speed estimation window halfsize
    trackparams(6)=1;
    % time delay between consecutive frames (ms)
    trackparams(7)=0;
end
alpha=-log(trackparams(1))/trackparams(4);
%% Tracking

load(fits_fname,'fits')

%filling the goodfit for the tracking function array
goodfits=zeros(sum(fits.goodfit),23);
goodfits(:,1)=fits.frame(fits.goodfit);%frame number
goodfits(:,9)=fits.row(fits.goodfit);%x position
goodfits(:,11)=fits.col(fits.goodfit);%y position
goodfits(:,14)=fits.sum(fits.goodfit);%intensity

trfile=Track_3D2(goodfits,trackparams(1),alpha,trackparams(3),trackparams(5),trackparams(6),...
    1,trackparams(7),trackparams(2));

if savetracks
    save([fits_fname,'_tracks'],'trfile','trackparams')
end

if ~isempty(trfile)
    %get rid of useless fits from the tracking program
    % tracks is made of 1: frame #, 2: x (px), 3: y (px), 4: track #
    tracks=trfile(:,[2,4,5,1]);
else
    warning('Not enough goodfits to contruct tracks')
    tracks=[];
end

end


function trk_filt= Track_filter(fits_fname,append_vec,trackparams,savetracks)
%% Track_filter
% written BPI 6/7/16
% Track_filter is a function to filter based on tracking. Currently just
% tracks and removes and the first and last frame from the track

% fits_fname is the name of the fits .mat file, importantly containing an
% array called fits with the fit information with column 9 being the
% goodfit boolean

% append_vec is a boolean determining whether trk_filt will be appended to
% the fit .mat file. Default is 0

% trackparams are the tracking parameters, definitions and defaults below

%%%% Outputs %%%%
% trk_filt is a logical vector which indicates whether or not the fit passed
% the filter
%
% tracks is array with the track information and columns
% tracks is made of 1: frame #, 2: x (px), 3: y (px), 4: track #

% the function also can append the fit .mat file with trk_filt

if nargin<2;append_vec=0;end
if nargin<4;savetracks=0;end

%% Default tracking parameters
if nargin<3
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

%% Track it
load(fits_fname,'fits')

%the logical vector of whether or not the fit passed the tracking filtering
trk_filt=false(size(fits.frame,1),1);

[~,fname] = fileparts(fits_fname);
disp(['Tracking ',fname]);

tracks = Tracking(fits_fname,trackparams,savetracks);
if ~isempty(tracks)
    save(fits_fname,'tracks','-append')
    
    % remove the first and last entries of each track
    for ii=1:max(tracks(:,4))
        %the molecules in the current track
        mols_inds=find(tracks(:,4)==ii);
        %remove the first and last entries
        tracks(mols_inds([1,end]),:)=[];
    end
    % trk_filt(ismember(fits(:,1:3),tracks(:,1:3),'rows'))=1;
    trk_filt(ismember([fits.frame,fits.row,fits.col],tracks(:,1:3),'rows'))=1;
else
    trk_filt=false(size(fits.frame,1),1);
end

if append_vec
    save(fits_fname,'trk_filt','-append')
end

end


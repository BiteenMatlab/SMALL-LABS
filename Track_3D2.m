function track_data=Track_3D2(fitfile1,min_merit,alpha,gamma,min_tr_length,speed_boxcar_halfsize,pxsize,....
    timedelay,itgtime)
% Written by David J. Rowland, The University of Michigan
%
% This code takes in the good fits file returned by the
% 'Astig3D_master_program.m' code and does single particle tracking in 3D
% by minimizing the distance between all possible pairs of fits.

% INPUTS:

% goodfitsfile_full_path: Full path of the .dat file that stores fitting
% parameters of good fits. See the 'Astig3D_master_program.m' code for the
% information stored in each column of the spreadsheet.

% track_file_path: Full path (including the file extension '.dat')
% specifying the location where the tracking file will be stored. It is
% recommended for the tracking file to take the form of [tifdirectory_name,
% '_sROI1_tROI1_tracking.dat'] for ease of batch data post-processing.

% sROI: Spatial ROI. Fits not belonging to this sROI won't be tracked.

% tROI: Temporal ROI. Fits not belonging to this tROI won't be tracked.

% min_merit: Minimal merit in pairing 2 fits. Any pair with a merit value
% smaller than this threshod won't be linked in any case.

% alpha: Penalty for displacements between 2 fits. The higher the alpha
% value, the lower the merit between 2 fits will be. Also, this value
% specifies the maximum allowable displacement (in px) of a molecule
% according to this relation: "alpha = -ln(min_merit)/max_step_size".

% gamma: Penalty for molecules disappearing between frames. Similar idea as
% 'alpha'.

% min_tr_length: Minimal track length (in # of frames, including dark
% frames). Tracks having fewer steps than this value won't be recorded.

% speed_boxcar_halfsize: Specify the half sliding box size that is used in
% calculating and smoothing the instantaneous speed, e.g.,
% "speed_boxcar_halfsize = 2" will smooth the instantaneous speed over 5
% frames;

% pxsize: Pixel size (nm/px).
% timedelay: Time delay between each frame, in ms.
% itgtime: Integration time of each frame, in ms.

% OUTPUTS:

% A tracking file is created at the designated location.

fitfile.data=fitfile1;

% fitfile.data=fitfile.data(abs(fitfile.data(:,22)-sROI)<0.01&...
%     abs(fitfile.data(:,23)-tROI)<0.1,:);
% Retain only data from the specified sROI and tROI.

% If the # of good fits is fewer than the minimum track length
if size(fitfile.data,1)>=min_tr_length
    enough_goodfits=1;
    
    % Total # of fits (of specified sROI and tROI) in the fit file
    num_spots=length(fitfile.data(:,1));
    
    % The frame number of the last frame
    max_fr=max(fitfile.data(:,1));
    
    % Create a (sparse) matrix for storing trajectories (each fit is
    % indexed into this variable according to which row it appears in the fit
    % file)
    track_row=spalloc(num_spots,max_fr,max([num_spots,max_fr]));
    
    sz_track=size(track_row);
    %% -----------------------------------------------------------------------%
    %  Construct Trajectories Based on Merit Functions
    %  -----------------------------------------------------------------------%
    
    % Trajectories can still be constructed (if merit is high enough)
    % even if a particle disappears temporarily, given that the number of
    % continuous frames during which this particle has disappeared does not
    % exceed the value of "max_dark_fr". Notice that here "log" is actually
    % "ln", and when max_dark_fr = 0 there is no dark frame tolerance.
    max_dark_fr=floor(-log(min_merit)/gamma);
    
    fr=fitfile.data(:, 1);
    
    %Frame change identifier. For example, if the frames containing
    % fits in specified ROIs are [1 3 3 3 4 6 6 9]', the frame change
    % identifier "fr_change_id" would be [1 1 0 0 1 1 0 1 1]'. Notice the
    % last "1" in the vector is a dummy purposely placed there such that in
    % the following loop the program knows how many row there are for the
    % last frame of current ROI.
    fr_change_id=vertcat(1,(fr(2:length(fr))-fr(1:length(fr)-1))~=0,1);
    
    
    % The row #s in the fit file that correspond to a change in frame
    % number. If the fr_change_id is [1 1 0 0 1 1 0 1]', "row_fr_change" would
    % be [1 2 5 6 8 9]', meaning rows 1, 2, 5, 6 and 8 in the fit file
    % correspond to spots appearing in a different frame than the previous
    % ones. The last element "9" is simply a result of the dummy # "1" that is
    % addressed above.
    row_fr_change=find(fr_change_id);
    
    delta_row_fr_change=row_fr_change(2:end)-row_fr_change(1:end-1);
    
    % This is used as an index to speed up the "max" command called
    % below. It does so by narrowing down the number of rows over which the
    % "max" command is performed.
    row_search_start=1;
    
    % Loop through rows in the fit file. In reality, it does not loop
    % row by row, but a group of rows (corresponding to spots from the
    % previous frame) to another group of rows (corresponding to spots from
    % the current frame), thus in essence, it loops frame by frame.
    for row_fr_change_id=2:length(row_fr_change)-1
        % List the row #s in the fit file of all spots from the previous
        % frame
        pre_fr_spots_row=row_fr_change(row_fr_change_id-1):...
            row_fr_change(row_fr_change_id-1)+...
            delta_row_fr_change(row_fr_change_id-1)-1;
        
        % List the row #s in the fit file of all spots from the current
        % frame
        curr_fr_spots_row=row_fr_change(row_fr_change_id):...
            row_fr_change(row_fr_change_id)+delta_row_fr_change(row_fr_change_id)-1;
        
        % The current frame number
        curr_fr_num=fitfile.data(curr_fr_spots_row(1),1);
        
        % This while loop is used to update the starting row number
        % from which the subsequent "max" is performed. For example, if
        % current frame number of 10 and "max_dark_fr" is 0, this loop
        % ensures that the subsequent "max" command does not look through
        % rows of data corresponding to frames earlier than frame#9, and
        % thus helps to reduce computational time.
        while curr_fr_num-fitfile.data(row_search_start)-1>max_dark_fr
            row_search_start=row_search_start+1;
        end
        
        % Frame #s, x-positions, y-positions, integrated intensities, cell
        % # indices and z-positions of spots from the current frame
        curr_fr_spots_data=fitfile.data(curr_fr_spots_row,[1 9 11 14 21 19]);
        
        % The maximum row number of each current track (i.e., the row
        % # of the most recent spot in each track or the "trajectory tail").
        max_row_curr_tr=nonzeros(max(track_row(row_search_start:...
            curr_fr_spots_row(end),1:curr_fr_num),[],2));
        
        % Find the corresponding frame number of the most recent spot
        % in each trajectory. Trajectories that have been broken for too long
        % (longer than "max_dark_fr") will be discarded and no longer be eligible
        % to compete for spots from the current frame. The "tailing spot" of
        % each remaining trajectory forms this variable "valid_max_row_tr".
        if ~isempty(max_row_curr_tr)
            valid_max_row_tr=((curr_fr_num-fitfile.data(max_row_curr_tr,1))...
                -1<=max_dark_fr).*max_row_curr_tr;
        else valid_max_row_tr=double.empty(0,1);
        end
        
        % Get rid of the zeros resulted from the previous logical statement.
        valid_max_row_tr=valid_max_row_tr(valid_max_row_tr~=0);
        
        % This includes all the candidate spots that will compete for
        % the spots from current frame. These candiate spots are either
        % "tailing spots" from trajecotries serveral frames ago, or all the spots
        % from the frame that immediately preceeds the current frame in the
        % fit file.
        pre_and_broken_row=unique(vertcat(pre_fr_spots_row',valid_max_row_tr));
        
        % Frame #s, x-positions, y-positions, integrated intensities,cell
        % number indices and z-positions of all candiate spots that will compete
        % for spots from the current frame.
        pre_and_broken_data=fitfile.data(pre_and_broken_row,[1,9,11,14,21,19]);
        
        % The matrix containing values of displacement between spots from
        % the current frame and the candidates spots from the previous frames Note
        % that trajectories are not constructed at this point, and thus the
        % displacement values are just "possible" values, i.e., all possible
        % combinations between the spots from different frames of interest
        % (current, previous, and broken frames) are considered. Combinations that
        % give smaller displacements will be favored later in deciding which pair
        % of spots are actually the same particle from different frames.
        disp_matrix=(((log(exp(pre_and_broken_data(:,2))...
            *exp(-curr_fr_spots_data(:,2)'))).^2)+...
            ((log(exp(pre_and_broken_data(:,3))...
            *exp(-curr_fr_spots_data(:,3)'))).^2)+...
            ((log(exp(pre_and_broken_data(:,6))...
            *exp(-curr_fr_spots_data(:,6)'))).^2)).^0.5;
        
        % Matrix of frame number (time) difference from all possible
        % combinations of spots between frames of interest. Again, "10^7"s are
        % inserted to prevent numbers from blowing up.
        fr_matrix=round(log(exp(-pre_and_broken_data(:,1)/10^7)...
            *exp(curr_fr_spots_data(:,1)'/10^7))*10^7-1);
        
        % This matrix will be used subsequently such that spots from
        % different cells (base on the phase mask of cells) will not form a
        % track.
        cellNo_matrix=round(log(exp(-pre_and_broken_data(:,5)/10^7)...
            *exp(curr_fr_spots_data(:,5)'/10^7))*10^7);
        cellNo_matrix(cellNo_matrix~=0)=1;
        cellNo_matrix=~cellNo_matrix;
        
        % Merit values for all possible combinations
        merit = exp(-alpha*disp_matrix).*exp(-gamma*fr_matrix).*cellNo_matrix;
        
        % Matrix of merit values from all possible combinations of
        % spots between frames of interest. All merit values below the
        % threshold "min_merit" are set to NaN.
        merit(merit<min_merit)=nan;
        
        % Maximization of merits. is the same thing as minimization of
        % penalty. the hungarian algorithm takes penalties instead of
        % merits, and the merit runs from 0 to 1, so:
        [assignment,~]=hungarian(1-merit);
        
        % The two-column matrix "row_spots_pair" contains the row #s
        % of spots from two different frames that most likely correspond to
        % the same particle, for example, the [2 12; 7 15] matrix means the
        % spot from row 2 in the earlier frame correspond to the spot from
        % row 12 in the current frame, and spots from row 7 and row 15 are
        % likely the same particle in different frames as well.
        [~,col_assignment,val_assignment]=find(assignment);
        row_spots_pair=[pre_and_broken_row(col_assignment)';...
            curr_fr_spots_row(val_assignment)]';
        
        % If pairing signals from different frames are detected
        % Check to see if spot pairs from different frames are found.
        if ~isempty(row_spots_pair)
            % Find the row#s and the indices of all tailing spots from
            % "track_row".
            % 'max_row_val': Actual values of non-zero elements (serves as row# indices
            % in 'track_row').
            % 'max_row_ind': Indices of locations of non-zero
            % elements (serves as column# indices in 'track_row').
            [max_row_val,max_row_ind]=...
                max(track_row(1:curr_fr_spots_row(end),1:curr_fr_num),[],2);
            
            % Create dummy indices from "row_spots_pair". For example, if
            % row_spots_pair = [2 12; 7 15], ia = [2; 7], which will directs the
            % trajectory "2-->12" to be appended to the 2nd row of "track_row", and the
            % trajectory "7-->15" to be appended to the 7th row of "track_row".
            % However, if the starting point of the latter trajectory (i.e., the spots
            % from row 7) is already a tailing spot of an earlier trajectory (e.g.,
            % "4-->6-->7"), the "7-->15" will be directly appended wherever that
            % earlier trajectory is in the "track_row" (here it is the 4th row) as
            % opposed to the default 7th row, and this is accomplished as coded below.
            max_row_val=full(max_row_val);
            ia=row_spots_pair(:,1);
            
            % Find whether the starting spots in "row_spots_pair" have already
            % been the tailing spots of earlier trajectories. If that is the case, the
            % indices "ia" have to be changed so that these spots won't be directly
            % appended to the default rows in "track_row".
            [intersect_val,intersect_ind,~]=intersect(max_row_val,row_spots_pair(:,1));
            
            % Change the "ia" indices of spots that are already tailing spots of
            % earlier trajectories so that they'll be appended to corresponding rows
            % where these trajectories already exist.
            [~,intersect_ind2,~]=intersect(ia,intersect_val);
            ia(intersect_ind2)=intersect_ind;
            
            % Row and column subscripts specifying where each of the fit in the linked
            % pair will be stored in 'track_row'.
            first_sub=[ia,max_row_ind(ia)];
            second_sub=[ia,max_row_ind(ia)+1];
            
            % Extract all row numbers that currently present in 'track_row', which will
            % be combined with data from the newly paired fits ('first_sub' and
            % 'second_sub') to re-generate an updated version of 'track_row'. I do this
            % instead of directly change the numbers in 'track_row' because it is
            % highly inefficient to directly change the content of a sparse matrix by
            % indexing as it alters the nonzero pattern of the sparse matrix.
            nonzero_value_track_row=nonzeros(track_row);
            nonzero_linearIndex_track_row=find(track_row);
            
            [nonzero_sub_r_track_row,nonzero_sub_c_track_row]=...
                ind2sub(sz_track,nonzero_linearIndex_track_row);
            
            update_r=vertcat(nonzero_sub_r_track_row,first_sub(:,1),second_sub(:,1));
            update_c=vertcat(nonzero_sub_c_track_row,first_sub(:,2),second_sub(:,2));
            update_value=vertcat(nonzero_value_track_row,row_spots_pair(:,1),...
                row_spots_pair(:,2));
            
            % Avoid updating the same element in 'track_row' twice (which
            % happens when a pair is appended to an existing track) by getting rid of
            % repeated rows of the update matrix ([row, column, value(row#)]).
            update_mat=[update_r,update_c,update_value];
            update_mat=unique(update_mat,'rows');
            update_r=update_mat(:,1);
            update_c=update_mat(:,2);
            update_value=update_mat(:,3);
            
            % row_fr_change_id % FOR DEBUGGING
            
            % Generate a new 'track_row' with information updated.
            track_row=sparse(update_r,update_c,update_value,...
                num_spots,max_fr,max([num_spots,max_fr]));
            % full_tr=full(track_row); % FOR DEBUGGING
        end
    end
else
    enough_goodfits=0;
end % If the fit file has enough fits from the specified s/t ROIs.

%% ------------------------------------------------------------------------
%  Truncate Short Trajectories and Generate Trajetory Data
%  ------------------------------------------------------------------------
% keep_track_row = track_row; % The original (untruncated) 'track_row'

% If there're enough good fits to construct 1 track
if enough_goodfits==1
    % First Step in getting rid of all tracks that don't meet the
    % required length. This step truncates rows within "track_row" that do not
    % contain any spots (i.e., all-zero row).
    track_row=track_row(track_row(:,1)~=0,:);
    
    % The number of frames (including dark frames) each track lasts.
    fr_last=(fitfile.data(nonzeros(max(track_row,[],2)),1)...
        -(fitfile.data(nonzeros(track_row(:,1)),1)))+1;
    
    % Get rid of tracks that are too short.
    track_row=track_row(fr_last>=min_tr_length,:);
    
    % The # of data points each trajectory, dark frames don't count.
    num_data_pts=sum((track_row)'~=0)';
    
    % Cumulative track length (or, more precisely, cumulative number of
    % data points). It is used as row indices to identify the end of a
    % trajectory in the "track_data" below.
    cum_tr_length=cumsum(num_data_pts);
    
    % Create a vector goes from 1 to the total number of trajectories
    % (e.g., [1; 2; 3; ... 50]). This will be used subsequently to assign a
    % track number in each trajectory.
    tr_num_ind=(1:length(cum_tr_length))';
    
    % Linearize "track_row"
    track_row=reshape(track_row',numel(track_row),1);
    track_row=track_row(track_row~=0);
    
    % Convert 'track_row' from sparse to a full vector
    track_row=full(track_row);
    
    % Create a matrix for storing trajecotry data. The 16 columns wil
    % be used to store (from left to right): (1) track # (2) frame # (3)
    % time(s) (4) x (px) (5) y (px) (6) displacements(px) (7) displacements(nm)
    % (8) x-error (px) (9) sROI# (10) instantaneous speed (11) y-error (px)
    % (12) tROI# (13) cell# (14) amplitude (15) z-position (px) (16) z-error
    % (px)
    track_data=zeros(size(track_row,1),16);
    
    track_data(cum_tr_length,1)=tr_num_ind;                         % track number
    cum_tr_length=vertcat(0,cum_tr_length);
    
    % Assign a track number for each trajectory
    for ii=1:length(cum_tr_length)-1
        track_data((cum_tr_length(ii)+1):cum_tr_length(ii+1),1)=...
            track_data(cum_tr_length(ii+1));
    end

    track_data(:,2)=fitfile.data(track_row,1);                      % Frame #
    track_data(:,3)=(track_data(:,2)-1).*(itgtime+timedelay)*0.001; % Time(s)
    track_data(:,4)=fitfile.data(track_row,9);                      % x-position(px)
    track_data(:,5)=fitfile.data(track_row,11);                     % y-position(px)
    track_data(:,15)=fitfile.data(track_row,19);                    % z-position(px)
    track_data(2:end,6)=real(sqrt(sum(track_data(2:end,[4,5,15])-...
        track_data(1:end-1,[4,5,15]),2)));                          % Displacement(px)
    track_data(cum_tr_length(1:end-1)+1,6)=nan;                     % ?
    track_data(:,7)=track_data(:,6)*pxsize;                         % Displacement(nm)
    track_data(:,8)=fitfile.data(track_row,10);                     % x-error (px)
    track_data(:,9)=fitfile.data(track_row,22);                     % sROI#
    track_data(:,10)=nan;                                           % Speed (TBD)
    track_data(:,11)=fitfile.data(track_row,17);                    % y-error
    track_data(:,12)=fitfile.data(track_row,23);                    % tROI#
    track_data(:,13)=fitfile.data(track_row,21);                    % Cell#
    track_data(:,14)=fitfile.data(track_row,3);                     % Amplitude
    track_data(:,16)=fitfile.data(track_row,20);                    % z-error(px)
    
    %% ------------------------------------------------------------------------
    %  Below: Calculate instantaneous speed
    %  ------------------------------------------------------------------------
    
    % Create a mask over which x- and y- displacements are to be
    % averaged respectively
    B=ones(2*speed_boxcar_halfsize+1,1);
    
    dx=track_data(2:end,4)-track_data(1:end-1,4);
    dy=track_data(2:end,5)-track_data(1:end-1,5);
    dz=track_data(2:end,15)-track_data(1:end-1,15);
    
    % Convolve x- y- and z-displacements using the mask B (i.e., sum
    % displacements over the range of the boxcar mask). Notice that dx_sum and
    % dy_sum will have a size that equals to the (size of "track_data" -
    % 2*speed_boxcar_halfsize - 1) because we don't use zero padding for the
    % convolution here and an additional row is omitted when taking the
    % differences dx and dy.
    dx_sum=conv2(dx,B,'valid');
    dy_sum=conv2(dy,B,'valid');
    dz_sum=conv2(dz,B,'valid');
    
    % Calculate the time difference over the boxcar range. For example,
    % with "speed_boxcar_halfsize" = 2 (mask size = 5 frames) and a time series
    % like [1 2 3 4 2 3]', the dt_boxcar would give "2"(3-1), corresponding to
    % the boxcar time difference centered on the 4th row of the time series
    % vector.
    dt_boxcar=track_data(2*speed_boxcar_halfsize+2:end,3)-...
        track_data(1:end-2*speed_boxcar_halfsize-1,3);
    
    conv_disp=conv2(track_data(:,7),...
        B,'valid');
    conv_disp=conv_disp(2: end);
    conv_disp(~isnan(conv_disp))=1;
    
    % Just to get rid of the time differences in "dt_boxcar" that
    % "cross-talk" between different tracks, with the help of the
    % displacment(nm) column (i.e., the 7th column in "track_data") to identify
    % the beginning of each track; For example, with a time series like [1 2 3
    % 4 2 3]', obviously the 5th element "2" indicates the starting time of a
    % new track, and with "speed_boxcar_halfsize" = 2 (mask size = 5 frames),
    % the time difference of "2" assigned to the 4th element does not mean
    % anything since it is produced by subtracting a time from one track to a
    % time from another track, needs to be set to "NaN";
    try
        dt_boxcar=conv_disp.*dt_boxcar;
        
        % Instantaneous speeds for x-, y- and z- displacements
        sp_x = dx_sum./dt_boxcar(1:end);
        sp_y = dy_sum./dt_boxcar(1:end);
        sp_z = dz_sum./dt_boxcar(1:end);
        
        % Instantaenous speed for the particle in um/s.
        velo = sqrt(sp_x.^2 + sp_y.^2 + sp_z.^2) * pxsize * 0.001;
    catch
        velo=nan;
    end
    % Instantaenous speed
    track_data(speed_boxcar_halfsize+2:end-speed_boxcar_halfsize,10)=velo;
else % If the # of good fits is fewer than the minimum track length
    fprintf(['The number of good fits is not enough to construct ', ...
        'a trajectory of minimal track length.\n']);
    track_data=[];
end
end
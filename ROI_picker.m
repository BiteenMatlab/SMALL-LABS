function ROI_picker( directoryname,mov_prd,prct_frms,append_str)
%% ROI_picker
% updated BPI 8/12/16
% This function allows a user to pick an ROI in a series of movies and
% write a new tif stack movie of just that ROI.

%   directoryname   is the name of the directory where the movies will be
% selected, if there is an error finding the directory the program will
% open uigetfile in the current working directory
%
% mov_prd is the period of movies that you will be prompted to
% choose/check the ROI of, e.g., if you choose 7 movies and mov_prd=3 then
% you will be shown movies 1,4,7. If you are choosing only one movie this
% input doesn't matter
%
% prct_frms is the percentage of frames that will be used when averaging.
% If you want to use the entire movie simply set prct_frms to 100. Default
% is 10
%
% append_str is the string which will be appended to the input movie name.
% Default is '_ROI'
%
% the function will also write a .txt file called filename_ROI_info.txt
% with the ROI coordinates.
%
%%%% optional use to specify ROI %%%%
% if you have the coordinates of the ROI that you want, you can input those
% instead of going through the clicking process. The way to do this is to
% input the coordinates in the place of mov_prd. The coordinates are input
% as [ row # right side , col # top ; row # left side , col # bottom]


%%%% Dependencies %%%%
% TIFFStack
% saveastiff

if nargin<3; prct_frms=10;end
if nargin<4; append_str='_ROI';end
%% Select the movies
display('Select the movies.')
try
    [datalist,dataloc,findex]=uigetfile([directoryname filesep '*.tif*'],'multiselect','on');
catch
    curdir=pwd;
    [datalist,dataloc,findex]=uigetfile([curdir filesep '*.tif*'],'multiselect','on');
end
if findex==0
    fprintf('no data selected\n')
    return
end
if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);


%% Pick the ROIs

if isscalar(mov_prd)
    
    %the movies to show
    movs_show=1:mov_prd:numel(dlocs);
    
    ROIs=cell(1,length(movs_show));
    
    for ii=1:length(movs_show);
        filename=[dlocs{movs_show(ii)},filesep,dnames{movs_show(ii)},'.tif'];
        %create A `TIFFStack` object  which behaves like a read-only memory
        %mapped TIFF file
        tfstk=TIFFStack(filename);
        movsz=size(tfstk);%the size of the movie
        %the frames to use
        frms_show=round(linspace(1,movsz(3),prct_frms/100*movsz(3)));
        
        movshow=mean(tfstk(:,:,frms_show),3);
        
        disp(dnames{movs_show(ii)});
        
        firstthru=1;%for showing the previous ROI
        ROItxtin = 110; %110 is the ASCII code for the letter n
        while ROItxtin ~= 121
            close(figure(11))
            figure(11);
            imshow(movshow,prctile(movshow(:),[.1,99.8]))
            %this doesn't work with small frame sizes, need to fix this
%             htit=title(['The new ROI will be a box \newline First click on the upper right corner',...
%                 '\newline Then click on the lower left corner']);
%             axpos = get(gca,'pos');
%             set(htit,'units','normalized');
%             extent = get(htit,'extent');
%             set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-0.66*extent(4)])
            
            if ii>1 && firstthru
                ROIclicks=ROIs{ii-1};
                ROIpos=[ROIclicks(2,1),ROIclicks(1,2),ROIclicks(1,1)-ROIclicks(2,1),...
                    ROIclicks(2,2)-ROIclicks(1,2)];
                
                hold on
                rectangle('Position',ROIpos,'EdgeColor','r')
                hold off
                
                ROItxtin=input('Does this ROI selection look OK to you (y/n)?  ','s');
                if isempty(ROItxtin);ROItxtin=110;end
                firstthru=0;
            end
            if ROItxtin ~= 121
                % Click and choose fiduciaries. ginput returns
                % [col_1,row_1;col_2,row_2;...;col_n,row_n];
                ROIclicks = round(ginput(2));
                %check the boundaries
                if ROIclicks(1,1)>movsz(2);ROIclicks(1,1)=movsz(2);end
                if ROIclicks(1,2)<1;ROIclicks(1,2)=1;end
                if ROIclicks(2,2)>movsz(1);ROIclicks(2,2)=movsz(1);end
                if ROIclicks(2,1)<1;ROIclicks(2,1)=1;end
                
                %the position vector for the rectangle object [x,y,w,h], where x,y
                %is the lower left corner
                ROIpos=[ROIclicks(2,1),ROIclicks(1,2),ROIclicks(1,1)-ROIclicks(2,1),...
                    ROIclicks(2,2)-ROIclicks(1,2)];
                
                hold on
                rectangle('Position',ROIpos,'EdgeColor','r')
                hold off
                
                ROItxtin=input('Does this ROI selection look OK to you (y/n)?  ','s');
                if isempty(ROItxtin);ROItxtin=110;end
            end
            
            ROIs{ii}=ROIclicks;
        end
    end
else
    ROIs=mov_prd;
end

%% Pull out the selection and write the new movies

for ii=1:numel(dlocs);
    tic;
    
    filename=[dlocs{ii},filesep,dnames{ii},'.tif'];
    [pathstr,fname,~] = fileparts(filename);
    
    disp(['Pulling out and writing ',fname,append_str])
    
    %create A `TIFFStack` object  which behaves like a read-only memory
    %mapped TIFF file
    tfstk=TIFFStack(filename);
    
    if isscalar(mov_prd)
        %which ROI click corresponds to this movie?
        ROIclicks=ROIs{find(movs_show<=ii,1,'last')};
    else
        ROIclicks=ROIs
    end
    
    newmov=tfstk(ROIclicks(1,2):ROIclicks(2,2),ROIclicks(2,1):ROIclicks(1,1),:);
    
    %%%%Save the movie%%%%
    options.overwrite=true;
    %save using saveastiff
    saveastiff(newmov,[pathstr,filesep,fname,append_str,'.tif'],options);
    
    tictoc=toc;%the time to run the entire program
    
    %save a .txt file with the parameters
    fileID = fopen([pathstr,filesep,fname,'_ROI_info.txt'],'w');
    fprintf(fileID,['ROI corners:\t',mat2str(ROIclicks),'\n']);
    fprintf(fileID,['time to run:\t',num2str(tictoc)]);
    fclose(fileID);
end



end


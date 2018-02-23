
filename='T:\Lab Members\Ben Isaacoff\Data\1_17_18\movies2Bfit.txt';

dir1='T:\Lab Members\Ben Isaacoff\Data\1_17_18';
dir2='T:\Lab Members\Ben Isaacoff\Data\1_22_18';
dir3='T:\Lab Members\Ben Isaacoff\Data\1_26_18';
dir4='T:\Lab Members\Ben Isaacoff\Data\1_27_18';


%% List all the movies
%dir1
files=dir([dir1,filesep,'mov_*.mat']);
files(cellfun('length',{files.name})~=11)=[];
%dir2
tfiles=dir([dir2,filesep,'mov_*.mat']);
tfiles(cellfun('length',{tfiles.name})~=11)=[];
%concatenate the file structures
files=[files;tfiles];
%dir3
tfiles=dir([dir3,filesep,'mov_*.mat']);
tfiles(cellfun('length',{tfiles.name})~=11)=[];
%concatenate the file structures
files=[files;tfiles];
%dir4
tfiles=dir([dir4,filesep,'mov_*.mat']);
tfiles(cellfun('length',{tfiles.name})~=11)=[];
%concatenate the file structures
files=[files;tfiles];

% check for outliers in filesize
allbytes=cell2mat({files.bytes});
outind=find(abs(allbytes-median(allbytes))>(0.1*median(allbytes)));
if ~isempty(outind)
    warning(['Possible non-movie file included in movie list: ',...
        files(outind).folder,filesep,files(outind).name])
end

%% Write the files in a list to a .txt file

fid=fopen(filename,'w');
for ii=1:length(files)
    fprintf(fid,'%s \r\n',[files(ii).folder,filesep,files(ii).name]);    
end
fclose(fid);

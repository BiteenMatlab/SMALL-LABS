%% SimulateDataPoisson
% last update 10/5/16 BPI
% A simple script to simulate molecule data with various backgrounds. Current
% implementation uses Poissonian probability distributions (poissrnd) for
% shot noise, and Guassian distributions for molecule brightness, and track
% length. Molecule positions are also random using matlab's randi function.
%
% You can change all of the parameters of the simulated movie by changing
% the parameters in the Control Panel below. The script will then save the
% simulated movie as a tiff stack. It will also save a .mat file with the
% ground truths of the molecules.
%
% You can check how accurate your fitting was by using the included script
% Check_Fitting_Simulated.

%% Control Panel

%filename to save tiff stack movie  and .mat file with parameter to, no
%file extension please, leave blank if you don't want to save this, like
%filename='';
directoryname='C:\Users\isaacoff\Documents\SMALL-LABS\Test data and simulations\';
filename='SimData_wNPs';

movsz=[256,256,500];%frame size then frame numbers

%background distribution
bgmean=1040;%BG mean
%with a Poissonian noise profile the noise amplitude with this background
%is sqrt(bgmean)

%molecule brightness distribution
molmean=65;%molecule mean
molstd=100;%molecule brightness std
molmin=40;%minimum brightness
constant_brightness=0;%if one then use a constant brightness and not a distribution

%molecule turning on distribution per frame
onmean=1;
onstd=3;

%molecule track distribution
trkmean=3;
trkstd=7;

%diffraction limit std in pixels
difstd=2;

%%%% Complex Backgrounds %%%%
%Settings for some constant PL GNRs
make_GNRs=1;%Boolean for GNR addition
numGNRs=7;%number of GNRs to add
%GNR brightness distribution
GNRmean=350;%GNR mean
GNRstd=100;%GNR std

%Settings for a large centered Gaussian background
make_GaussBG=0;%Boolean for Gaussian BG addition
gaussbrght=10000;%brightness of Gaussian BG
gaussstd=200;%std of Gaussian BG

%% Make the background

%make the background frames
mov=zeros(movsz)+bgmean;

if make_GNRs
    %what is their positions
    GNRr=randi(round([difstd*5,movsz(1)-difstd*5]),numGNRs,1);
    GNRc=randi(round([difstd*5,movsz(2)-difstd*5]),numGNRs,1);
    %what is their brightnesses
    GNRbrght=round(abs(GNRstd.*randn(numGNRs,1)+GNRmean));
    for ii=1:numGNRs
        %make the GNR, using an ~5*difstd size square
        hlfpts=round(5*difstd/2);
        [mshx,mshy]=meshgrid((-hlfpts:hlfpts),(-hlfpts:hlfpts));
        molecule=round(GNRbrght(ii).*exp(-mshx.^2/(2*difstd^2)-mshy.^2/(2*difstd^2)));
        %add the GNR brightness to the movie at the specified location
        mov(GNRr(ii)+(-hlfpts:hlfpts),GNRc(ii)+(-hlfpts:hlfpts),:)=...
            bsxfun(@plus,mov(GNRr(ii)+(-hlfpts:hlfpts),GNRc(ii)+(-hlfpts:hlfpts),:),molecule);
    end
else
    GNRr=0;GNRc=0;GNRbrght=0;    
end

if make_GaussBG
    %make the spot
    [mshx,mshy]=meshgrid(1:movsz(1),1:movsz(2));
    molecule=round(gaussbrght.*exp(-(mshx-movsz(1)/2).^2/(2*gaussstd^2)-(mshy-movsz(2)/2).^2/(2*gaussstd^2)));
    mov=bsxfun(@plus,mov,molecule);
    
end

%% Add the molecules

% for the .mat file, columns:
% 1. frame number, 2. row (px),3. col (px), 4. width (px) ,5. offset,
% 6. amplitude, 7. integral (sum) of Gaussian
sim_mols=zeros(1,7);

%prepare the meshgrid for the molecule
hlfpts=round(5*difstd/2);
[mshx,mshy]=meshgrid((-hlfpts:hlfpts),(-hlfpts:hlfpts));

for ii=1:(movsz(3))
    %how many molecules in this frame
    molon=round(abs(onstd.*randn()+onmean));
    %what is their positions
    molr=randi(round([difstd*5,movsz(1)-difstd*5]),molon,1);
    molc=randi(round([difstd*5,movsz(2)-difstd*5]),molon,1);
    
    %what is their brightnesses
    brghts=round(molstd.*randn(molon,1)+molmean);    
    %dealing with any below the minimum brightness
    while any(brghts<molmin)
        brghts(brghts<molmin)=round(molstd.*randn(sum(brghts<molmin),1)+molmean);        
    end
    if constant_brightness
       brghts=repmat(molmean,size(brghts));
    end
    %what is their track lengths
    trklens=round(abs(trkstd.*randn(molon,1)+trkmean));
    
    %add them to the movie
    for jj=1:molon        
        %make sure none of the molecules are too close
        tooclose=1;%the loop flag
        while tooclose==1
            framer=sim_mols(sim_mols(:,1)==ii,2);framec=sim_mols(sim_mols(:,1)==ii,3);
            dists=(abs(framer-molr(jj))<=(difstd*5)) & (abs(framec-molc(jj))<=(difstd*5));
            
            if any(dists)
                %remake the position
                molr(jj)=randi(round([difstd*5,movsz(1)-difstd*5]));
                molc(jj)=randi(round([difstd*5,movsz(2)-difstd*5]));
            else
                tooclose=0;%leave the while loop
            end    
        end        
        
        %make the molecule, using an ~5*difstd size square        
        molecule=round(brghts(jj).*exp(-mshx.^2/(2*difstd^2)-mshy.^2/(2*difstd^2)));
        %integrate
        mol_int=sum(molecule(:));
        %change the track length if it would have the molecule on longer
        %than the movie
        if (trklens(jj)+ii)>movsz(3);trklens(jj)=movsz(3)-ii;end
        %add the molecule brightness to the movie at the specified location
        %and time points
        mov(molr(jj)+(-hlfpts:hlfpts),molc(jj)+(-hlfpts:hlfpts),ii+(0:trklens(jj)))=...
            bsxfun(@plus,mov(molr(jj)+(-hlfpts:hlfpts),molc(jj)+(-hlfpts:hlfpts),ii+(0:trklens(jj))),molecule);
        
        %add it to the sim_mols array
        sim_mols=cat(1,sim_mols,[(ii+(0:trklens(jj)))',repmat(molr(jj),trklens(jj)+1,1),...
            repmat(molc(jj),trklens(jj)+1,1),repmat(difstd,trklens(jj)+1,1),zeros(trklens(jj)+1,1),...
            repmat(brghts(jj),trklens(jj)+1,1),repmat(mol_int,trklens(jj)+1,1)]);
    end
end
sim_mols=sim_mols(2:end,:);%get rid of first row of zeros

mov=poissrnd(mov);


%% Save it
if ~isempty(filename)
    %save all the variables and the sim_mols array
    save([directoryname,filename,'_params'],'bgmean','molmean','molstd','onmean','onstd','trkmean',...
        'trkstd','difstd','make_GNRs','numGNRs','GNRmean','GNRstd','make_GaussBG','gaussbrght','gaussstd',...
        'GNRc','GNRr','GNRbrght','sim_mols');
    
    %%%%Save the movie%%%%
    options.overwrite=true;
    %save using saveastiff
    saveastiff(uint16(mov), [directoryname,filename,'.tif'],options);
end



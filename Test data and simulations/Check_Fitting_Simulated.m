%% Check_Fitting_Simulated
%%%control panel%%%
% The only section you should have to change anything
%Filename assuming that you're using the defualt filenames in SMALL-LABS
%The directory where the data is (with a \ at the end)
directory_name='C:\Users\isaacoff\Documents\SMALL-LABS\Test data and simulations\';

%the name of the movie
namebase='SimData_wGNRs';

%background subtraction?
bgsub=1;

%This lets you click through the frames to see the individual fits
viewfits_bool=0;

%this step can take a while, it is only necessary to do once for each
%dataset. So if you've already run this script with this made, then no need
%to do it again.
make_corgess=1;%make the corgess vector

%plot commands, these booleans let you plot different metrics of the
%fitting
cohort_stats=0;%cohort statistics on the intensity
indv_cmprs=1;%individual comparisons
%options for the individual comparisons
percent_error=1;%show percent errors for each fit
cohort_error=0;%show percent error based on mean of simulated population

%% Loading in the data

simfname=[directory_name,namebase,'_params.mat'];
load(simfname,'sim_mols');

if bgsub
    ges_fname=[directory_name,namebase,'_avgsub_guesses.mat'];
    fits_fname=[directory_name,namebase,'_AccBGSUB_fits.mat'];
else
    ges_fname=[directory_name,namebase,'_guesses.mat'];
    fits_fname=[directory_name,namebase,'_fits.mat'];
end

load(ges_fname,'guesses','dfrlmsz','corgess','sim2gess')
load(fits_fname,'fits','trk_filt')

movfname=[directory_name,namebase,'.tif'];

%create A `TIFFStack` object  which behaves like a read-only memory
%mapped TIFF file
tfstk=TIFFStack(movfname);
movsz=size(tfstk);%the size of the movie

%% Compare with simulated data
if make_corgess
    corgess=zeros(size(guesses,1),1);%correct guesses
    %indices of simulated molecules corresponding to correct guesses
    sim2gess=zeros(size(guesses,1),1);
    
    frms=unique(sim_mols(:,1));
    
    ctoff=dfrlmsz/2;%cutoff distance
    intsct=0;%intersection
    idx=0;%gess index
    for ii=1:size(guesses,1)
        % for ii=1:100 %debugging
        gess=guesses(guesses(:,1)==ii,:);
        grndtrth=sim_mols(sim_mols(:,1)==ii,:);%simulated ground truth
        
        crfges=gess;crfgrt=grndtrth;%temp for this frame
        
        intsctfrm=0;
        for jj=1:size(gess,1)
            idx=idx+1;
            for kk=1:size(grndtrth,1)
                if  abs(crfges(jj,2)-crfgrt(kk,2))<ctoff && abs(crfges(jj,3)-crfgrt(kk,3))<ctoff
                    corgess(idx)=1;%adding it to the correct guess list
                    %filling the sim2gess list, the (1) is for the rare case
                    %when two molecules overlap. It's arbitrary and will
                    %mess up the intenstiy comparisons, but it's so rare
                    %that it shouldn't matter. Note that the new
                    %simulations have no overlapping molecules
                    idxs=find(sim_mols(:,1)==ii & sim_mols(:,2)==crfgrt(kk,2) & ...
                        sim_mols(:,3)==crfgrt(kk,3));
                    sim2gess(idx)=idxs(1);
                    
                    crfgrt(kk,:)=NaN;%remove it from the list to avoid double counting
                    crfges(jj,:)=NaN;%remove it from the list to avoid double counting
                    intsctfrm=intsctfrm+1;
                end
            end
        end
        intsct=intsct+intsctfrm;
        
        if viewfits_bool
            curfrm=double(tfstk(:,:,ii));
            imshow(curfrm,prctile(curfrm(curfrm>0),[.1,99.8]))
            
            if ~isempty(grndtrth)
                vcs=viscircles([grndtrth(:,3),grndtrth(:,2)],repmat(dfrlmsz+2,[length(grndtrth(:,3)),1]));
                set(vcs.Children,'LineWidth',1,'Color','blue')
            end
            if ~isempty(gess)
                vcs=viscircles([gess(:,3),gess(:,2)],repmat(dfrlmsz,[length(gess(:,3)),1]));
                set(vcs.Children,'LineWidth',1,'Color','red')
            end
            title(['frame ',num2str(ii)])
            
            disp(['frame ',num2str(ii),' fn = ',num2str(size(grndtrth,1)-intsctfrm)])
            keyboard
        end
    end
    corgess=logical(corgess);%converting to a logical vector
    
    save(ges_fname,'corgess','sim2gess','-append')
    
    
end

intsct=size(find(corgess),1);

jaccard=intsct/(size(guesses,1)+size(sim_mols,1)-intsct);
fp=(size(guesses,1)-intsct)/size(guesses,1);
fn=(size(sim_mols,1)-intsct)/size(sim_mols,1);

disp(['Guess Jaccard Index = ',num2str(jaccard)])
disp(['Guess false positive rate = ',num2str(fp)])
disp(['Guess false negative rate = ',num2str(fn)])

% break %debugging

% end
%% Intensity comparisons

intsct=size(find(corgess & fits.goodfit),1);

jaccard=intsct/(size(find(fits.goodfit),1)+size(sim_mols,1)-intsct);
fp=(size(find(fits.goodfit),1)-intsct)/(size(find(fits.goodfit),1));
fn=(size(sim_mols,1)-intsct)/size(sim_mols,1);

disp(['Fit Jaccard Index = ',num2str(jaccard)])
disp(['Fit false positive rate = ',num2str(fp)])
disp(['Fit false negative rate = ',num2str(fn)])

% break %debugging

%cohort statistics
if cohort_stats
    % % histogram of simulated molecules
    %histogram binning
    simm=round(prctile(sim_mols(:,6),98));
    fitm=round(prctile(gdfts(:,6),98));
    
    figure
    histogram(sim_mols(:,6),(0:round(simm/100):simm))
    xlim([0,simm])
    xlabel('Gaussian Amplitude')
    ylabel('Occurrences')
    title('Simualted Molecules Gaussian Amplitudes')
    %histogram of goodfits
    figure
    histogram(gdfts(:,6),(0:round(fitm/100):fitm))
    xlim([0,fitm])
    xlabel('Gaussian Amplitude')
    ylabel('Occurrences')
    title('Good Fits Gaussian Amplitude')
    
    % histogram of simulated molecules
    %histogram binning
    simm=round(prctile(sim_mols(:,7),98));
    fitm=round(prctile(gdfts(:,8),98));
    
    figure
    histogram(sim_mols(:,7),(0:round(simm/100):simm))
    xlim([0,simm])
    xlabel('Integral')
    ylabel('Occurrences')
    title('Simualted Molecules Integrals')
    %histogram of goodfits
    figure
    histogram(gdfts(:,8),(0:round(fitm/100):fitm))
    xlim([0,fitm])
    xlabel('Integral')
    ylabel('Occurrences')
    title('Good Fits Integrals')
    
    disp(['Gaussian Amplitude ',char(10),' simulated mean = ',num2str(mean(sim_mols(:,6))),...
        '    fit mean = ',num2str(mean(gdfts(:,6))),char(10),' simulated median = ',...
        num2str(median(sim_mols(:,6))),'    fit median = ',num2str(median(gdfts(:,6))),...
        char(10),' simulated mode = ',num2str(mode(sim_mols(:,6))),'    fit mode = ',...
        num2str(mode(gdfts(:,6)))])
    
    disp(['Gaussian Integral ',char(10),' simulated mean = ',num2str(mean(sim_mols(:,7))),...
        '    fit mean = ',num2str(mean(gdfts(:,8))),char(10),' simulated median = ',...
        num2str(median(sim_mols(:,7))),'    fit median = ',num2str(median(gdfts(:,8))),...
        char(10),' simulated mode = ',num2str(mode(sim_mols(:,7))),'    fit mode = ',...
        num2str(mode(gdfts(:,8)))])
    
end



if indv_cmprs
    %Individiual molecule comparisons
    % 1. x (px), 2. y (px), 3. width (px) ,4. offset,
    % 5. amplitude, 6. integral (sum) of Gaussian
    
    sim_gdfts=sim_mols(sim2gess(corgess & fits.goodfit),2:7);
%     [2:6,8]
    cor_gdfts=[fits.row(corgess & fits.goodfit),fits.col(corgess & fits.goodfit),...
        fits.widthr(corgess & fits.goodfit),fits.offset(corgess & fits.goodfit),...
        fits.amp(corgess & fits.goodfit),fits.sum(corgess & fits.goodfit)];
    
    dif=cor_gdfts-sim_gdfts;
    
    %The Gaussian volume
    sim_vol=sim_gdfts(:,5).*(sim_gdfts(:,3).^2);
    fit_vol=cor_gdfts(:,5).*(cor_gdfts(:,3).^2);
    
    figure
    histdat=dif(:,1);
    histlim=prctile(histdat,[2,98]);
    histogram(histdat,linspace(histlim(1),histlim(end),100))
    xlabel('X difference (pixels)')
    ylabel('Occurrences')
    title('X position (fit - simualted)')
    
    figure
    histdat=dif(:,2);
    histlim=prctile(histdat,[2,98]);
    histogram(histdat,linspace(histlim(1),histlim(end),100))
    xlabel('Y difference (pixels)')
    ylabel('Occurrences')
    title('Y position (fit - simualted)')
    
    figure
    histdat=sqrt(dif(:,1).^2+dif(:,2).^2);
    histlim=prctile(histdat,[0,98]);
    histogram(histdat,linspace(0,histlim(end),100))
    xlabel('magnitude difference (pixels)')
    ylabel('Occurrences')
    title('magnitude position (fit - simualted)')
    
    if percent_error
        figure
        histdat=100*dif(:,3)./sim_gdfts(:,3);
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Width % error (pixels)')
        ylabel('Occurrences')
        title('Gaussian std % error (fit - simualted)')
        axis tight
        
        %         figure
        %         histdat=100*dif(:,4)./sim_gdfts(:,5);
        %         histlim=round(prctile(histdat,[2,98]));
        %         histogram(histdat,linspace(histlim(1),histlim(end),100))
        %         xlabel('Offset % error (intensity)')
        %         ylabel('Occurrences')
        %         title('Offset intensity % error (offset/amplitude) (fit - simualted)')
        %         axis tight
        
        figure
        histdat=100*dif(:,5)./sim_gdfts(:,5);
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Amplitude % error (intensity)')
        ylabel('Occurrences')
        title('Gaussian amplitude % error (fit - simualted)')
        axis tight
        
        figure
        histdat=100*dif(:,6)./sim_gdfts(:,6);
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Integral % error (intensity)')
        ylabel('Occurrences')
        title('Integral % error (fit - simualted)')
        axis tight
        
        figure
        histdat=100*(fit_vol-sim_vol)./sim_vol;
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Gaussian Volume % error (intensity)')
        ylabel('Occurrences')
        title('Gaussian Volume % error (fit - simualted)')
        axis tight
        
    elseif cohort_error
        figure
        histdat=100*dif(:,3)/mean(sim_gdfts(:,3));
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Width % error (pixels)')
        ylabel('Occurrences')
        title('Gaussian std cohort % error (fit - simualted)')
        axis tight
        
        %         figure
        %         histdat=100*dif(:,4)/mean(sim_gdfts(:,5));
        %         histlim=round(prctile(histdat,[2,98]));
        %         histogram(histdat,linspace(histlim(1),histlim(end),100))
        %         xlabel('Offset % error (intensity)')
        %         ylabel('Occurrences')
        %         title('Offset intensity cohort % error (offset/amplitude) (fit - simualted)')
        %         axis tight
        
        figure
        histdat=100*dif(:,5)/mean(sim_gdfts(:,5));
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Amplitude % error (intensity)')
        ylabel('Occurrences')
        title('Gaussian amplitude cohort % error (fit - simualted)')
        axis tight
        
        figure
        histdat=100*dif(:,6)/mean(sim_gdfts(:,6));
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Integral % error (intensity)')
        ylabel('Occurrences')
        title('Integral cohort % error (fit - simualted)')
        axis tight
        
        figure
        histdat=100*(fit_vol-sim_vol)./mean(sim_vol);
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Gaussian Volume % error (intensity)')
        ylabel('Occurrences')
        title('Gaussian Volume cohort % error (fit - simualted)')
        axis tight
        
    else
        figure
        histdat=dif(:,3);
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Width difference (pixels)')
        ylabel('Occurrences')
        title('Gaussian std difference (fit - simualted)')
        axis tight
        
        %         figure
        %         histdat=dif(:,4);
        %         histlim=round(prctile(histdat,[2,98]));
        %         histogram(histdat,linspace(histlim(1),histlim(end),100))
        %         xlabel('Offset difference (intensity)')
        %         ylabel('Occurrences')
        %         title('Offset intensity difference (fit - simualted)')
        %         axis tight
        
        figure
        histdat=dif(:,5);
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Amplitude difference (intensity)')
        ylabel('Occurrences')
        title('Gaussian amplitude difference (fit - simualted)')
        axis tight
        
        figure
        histdat=dif(:,6);
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Integral difference (intensity)')
        ylabel('Occurrences')
        title('Integral difference (fit - simualted)')
        axis tight
        
        figure
        histdat=fit_vol-sim_vol;
        histlim=round(prctile(histdat,[2,98]));
        histogram(histdat,linspace(histlim(1),histlim(end),100))
        xlabel('Gaussian Volume difference (intensity)')
        ylabel('Occurrences')
        title('Gaussian Volume difference (fit - simualted)')
        axis tight
    end
end






















function fits_mat=Legacy_fits_strcut2mat(directoryname,save_it)
%% Legacy_fits_strcut2mat 
% This is a function to change the fits structure output from
% Subtract_then_fit to a matrix in the same configuration as the fits
% matrix output in legacy versions of Subtract_then_fit. If you need this
% functionality simply ctrl+f replace fits with fits_mat in your analysis
% code.
%
% directoryname   is the name of the directory where the fits .mat files
% will be selected, if there is an error finding the directory the program
% will open uigetfile in the current working directory
%
%  save_it is a boolean asking whether or not you want to save the fits_mat
%  matrix to the fits .mat file
%
% The configuration of the legacy fits matrix is
% 1. frame number, 2. row pos (px),3. column pos (px), 4. width (px) ,5. offset,
% 6. amplitude, 7. error, 8. sum(:), 9. goodfit boolean

%default value of save it
if nargin<2;save_it=1;end

%% Select the fits .mat files
display('Select the fit .mat files.')
try
    [datalist,dataloc,findex]=uigetfile([directoryname filesep '*.mat*'],'multiselect','on');
catch
    curdir=pwd;
    [datalist,dataloc,findex]=uigetfile([curdir filesep '*.mat*'],'multiselect','on');
end
if findex==0
    fprintf('no data selected\n')
    return
end
if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);

%% Make it and save it

for ii=1:numel(dlocs)
     fits_fname=[dlocs{ii},filesep,dnames{ii}];
    load(fits_fname,'fits');
    
    fits_mat=[fits.frame,fits.row,fits.col,fits.widthc,fits.offset,fits.amp,...
        fits.err,fits.sum,fits.goodfit];
    
    if save_it
        save(fits_fname,'fits_mat','-append')%append the new matrix
    end
end


end


function Change_GoodFits(directoryname,new_maxerr,new_stdtol,new_dfrlmsz,new_citol)
%% Change_GoodFits 
% This is a function to change the goodfits vector in the fits array. If
% you don't want to change either maxerr or stdtol, then enter the
% corresponding input (new_maxerr or new_stdtol) as []. See either
% Subtract_then_fit or the User Guide for an explanation of what
% maxerr and stdtol are.
%
% directoryname   is the name of the directory where the fits .mat files
% will be selected, if there is an error finding the directory the program
% will open uigetfile in the current working directory
%
%  new_maxerr is the new max_err that you want to use. If you don't want to
%  change this then enter [] for new_maxerr
%
%  new_stdtol is the new stdtol that you want to use. If you don't want to
%  change this then enter [] for new_stdtol
%
%  new_dfrlmsz is the new dfrlmsz that you want to use. If you don't want to
%  change this then enter [] for new_dfrlmsz
%
%  citol is the modifier to tighten (<1) or loosen (>1) constraints on the
%  confidence interval limit. If you don't want to change this then enter 
%  [] for citol
%
%     Copyright (C) 2017  Benjamin P Isaacoff
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

%% Make the new good_fit vector

h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,'Changing GoodFits');

for ii=1:numel(dlocs);
    try; waitbar(ii/numel(dlocs),h1); end
    
    fits_fname=[dlocs{ii},filesep,dnames{ii}];
    load(fits_fname,'fits','stdtol','maxerr','dfrlmsz')%load in the necessary variables
    % change the parameters if they were sent
    if ~isempty(new_maxerr);maxerr=new_maxerr;end
    if ~isempty(new_stdtol);stdtol=new_stdtol;end
    if ~isempty(new_dfrlmsz);dfrlmsz=new_dfrlmsz;end
    if ~isempty(new_citol);citol=new_citol;end
    if isempty(new_citol);citol=1;end
    
    %the conversion between dfrlmsz and the STD of the Gaussian, reccomended
    %using the full width at 20% max given by (2*sqrt(2*log(5)))
    dfD2std=(2*sqrt(2*log(5)));
    %the guessed std
    gesss=dfrlmsz/dfD2std;
    
    %%% make the new goodfits vector %%%
    goodfits=(mean([fits.widthc,fits.widthr],2)<=(stdtol*gesss) & mean([fits.widthc,fits.widthr],2)>=(gesss/stdtol)) & ...
        all([fits.row,fits.col,fits.amp,fits.sum]>=0,2) & fits.amp<fits.sum & fits.err>maxerr & fits.rowCI<=citol*dfrlmsz ...
        & fits.colCI<=citol*dfrlmsz;
    
    %update the fits array
    fits.goodfit=goodfits;
    %save it
    save(fits_fname,'fits','stdtol','maxerr','dfrlmsz','-append')%append the new vector and the used parameters
end

%close the waitbar
try; close(h1); end
end


function [fitPars, conf95, guesses, outPut]=gaussFit(img, varargin)
%
% NAME:
%       gaussFit
% PURPOSE:
%       Fits a generalized gaussian function to 2d imaging data. This code
%       produces results in units of pixels for the center position and
%       widths.
% CATEGORY:
%       Image Processing
% CALLING SEQUENCE:
%       [fitPars, conf95] = gaussFit(img,findTheSpot);
% INPUTS:
%       img:            The two-dimensional array to be fit to a gaussian
%
%       varargin:       use paired inputs to set the property (input 1) to the
%           value (input 2) desired.
%
%       Properties:     Descriptions:
%
%       searchBool:    1 or 0. Default behavior is to fit an
%           ROI in the center of the image. If the spot is not near the
%           center or the image is very large, findTheSpot enables the code
%           to first roughly locate the spot and then use that location as
%           the ROI center.
%
%       showGuessing:   1 or 0. show output. default is 0.
%
%       widthGuess:     set the starting value for the width of the
%                       Gaussian in units of pixels.
%
%       nPixels         pixel width of ROI to be selected from img. default
%                       is 11. the value should be odd.
%
% OUTPUTS:
%       fitPars:        fitting coefficient vector, units are pixels.
%       fitCI:          95% confidence interval of fitting coefficients at
%                       end of fitting
% PROCEDURE:
%       1. Peak guessing and/or data ROI selection of local area inside img
%       2. Non-linear least squares minimization for 7 (or 6 or 5) -
%       parameter Gaussian function on the ROI selected.
%
% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 3/16.
% NOTES:
%       This code 'gaussFit.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.
%
%       For testing purposes, run this script:
%
%       img = exp(-x.^2/2/2^2-y.^2/2/3^2)+.02*randn(size(x));
%       p = gaussFit(img,'widthGuess',2);
opts = optimset('Display','off');
warning('off','MATLAB:singularMatrix');
% warning('off','all')
imSize = size(img);

%% default parameters
% peak guessing parameters
params.spotSizeLB = 1.2;
params.spotSizeUB = 10;
params.intThresh = 300;
params.lZero = 10;
params.hMax = 200;

% other parameters
params.searchBool = 1;
params.checkVals = 0;
params.widthGuess = 2;
params.frameNumber = 1;

% fitting window width; should be odd valued
params.nPixels = 11;

% 3 is a symmetric gaussian (5 parameters) fit
% 2 is a fixed angle asymmetric gaussian fit
% 1 is a 7 parameter asymmetric gaussian fit
params.ffSwitch = 3;

fNames=fieldnames(params);

% if any sim parameters are included as inputs, change the simulation
% parameters mentioned
if nargin>1
    for ii=1:2:nargin-2
        whichField = strcmp(fNames,varargin{ii});
        
        if all(~whichField)
            warning('Check spelling. Parameter change may have not occurred.')
        end
        
        eval(['params.' fNames{whichField} ' = varargin{ii+1};'])
    end
end

%% fitting functions
switch params.ffSwitch
    case 1
        % freely rotating bivariate gaussian function for least squares minimization
        % parameters: [xCenter, yCenter, angle, xSD, ySD, amplitude, offset]
        xR=@(x,y,xc,yc,th)(x-xc)*cos(th)-(y-yc)*sin(th);
        yR=@(x,y,xc,yc,th)(x-xc)*sin(th)+(y-yc)*cos(th);
        fFun=@(p,X) exp( -xR(X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(4)^2 + ...
            -yR( X(:,1), X(:,2), p(1), p(2), p(3)).^2/2/p(5)^2 ) *p(6) + p(7);
        pStart = [];
        lb = [];
        ub = [];
    case 2
        % fixed angle fit
        % parameters: [xCenter, yCenter, xSD, ySD, amplitude, offset]
        th = 0;
        xR = @(x,y,xc,yc)(x-xc)*cos(th)-(y-yc)*sin(th);
        yR = @(x,y,xc,yc)(x-xc)*sin(th)+(y-yc)*cos(th);
        fFun = @(p,X) exp( -xR(X(:,1), X(:,2), p(1), p(2)).^2/2/p(3)^2 + ...
            -yR( X(:,1), X(:,2), p(1), p(2)).^2/2/p(4)^2 ) *p(5) + p(6);
        pStart = [];
        lb = [-inf, -inf, 0, 0, -inf, -inf];
        ub = [inf, inf, inf, inf, inf, inf];
    case 3
        % symmetric gaussian
        fFun = @(p,X) exp( -((X(:,1)-p(1)).^2 + (X(:,2)-p(2)).^2)/2/p(3).^2) * p(4) + p(5);
        pStart = [];
        lb = [-params.nPixels, -params.nPixels, 0, 0, -inf];
        ub = [params.nPixels*2, params.nPixels, inf, inf, inf];
end

% starting guesses
pStart(1) = 0;
pStart(2) = 0;
pStart(3) = params.widthGuess;

%% rough localization of molecules
if params.searchBool
    % band pass
    bIm = bpassDJR(img, params.spotSizeLB, params.spotSizeUB, params.intThresh, params.lZero);
    
    % watershed
    extImg = imextendedmax(bIm,params.hMax);
    
    % failed watershed can result in all ones
    if all(extImg(:))
        extImg = extImg-1;
    end
    
    % shrink to a point. this is the estimated location of the spot
    sIm = bwmorph(extImg,'shrink',inf);
    
    % if shrinking the image produces rings, remove the rings
    cc = bwconncomp(sIm);
    if cc.NumObjects < sum(sIm(:))
        whichBad = cellfun(@numel,cc.PixelIdxList) > 1;
        sIm(cc.PixelIdxList{whichBad}) = 0;
    end
    
    if params.checkVals
        subplot(2,4,2)
        imshow(img,[]);
        
        subplot(2,4,3)
        imshow(bIm,[])
        
        subplot(2,4,6)
        imshow(extImg,[])
        
        subplot(2,4,7)
        imshow(sIm,[])
    end
    
    % the index of the one pixel is a good guess for the particle location
    [guesses(:,1),guesses(:,2)] = find(sIm);
else
    % otherwise, assume the spot is in near the center of the image
    guesses = round(size(img)/2);
end

% number of fits
nFits = size(guesses,1);

% output initialization
if nFits > 50
    warning(['way too many fits in frame number ' num2str(params.frameNumber)])
    fitPars = nan(1,5);
    conf95 = nan(1,5);
    guesses = nan(1,2);
    outPut.firstorderopt = [];
    outPut.iterations = [];
    outPut.funcCount = [];
    outPut.cgiterations = [];
    outPut.algorithm = [];
    outPut.stepsize = [];
    outPut.message = [];
    return
elseif nFits > 0
    fitPars = nan(nFits,numel(lb));
    conf95 = nan(nFits,numel(lb));
    
    outPut(nFits).firstorderopt = [];
    outPut(nFits).iterations = [];
    outPut(nFits).funcCount = [];
    outPut(nFits).cgiterations = [];
    outPut(nFits).algorithm = [];
    outPut(nFits).stepsize = [];
    outPut(nFits).message = [];
else
    fitPars = nan(1,5);
    conf95 = nan(1,5);
    guesses = nan(1,2);
    outPut.firstorderopt = [];
    outPut.iterations = [];
    outPut.funcCount = [];
    outPut.cgiterations = [];
    outPut.algorithm = [];
    outPut.stepsize = [];
    outPut.message = [];
end

% pad the img(s) with nans (removed at end).
padsize = params.nPixels(ones(1,2));
img = padarray(img,padsize,nan,'both');
guesses = guesses+params.nPixels;

%% fit the data
% fitting domain
[x,y] = ndgrid(1:params.nPixels, 1:params.nPixels);
X = cat(2,x(:),y(:)) - params.nPixels/2 - .5;

for ii = 1:nFits
    % find the selection domain
    [sDom1,sDom2] = ndgrid(guesses(ii,1)-(params.nPixels-1)/2:guesses(ii,1)+(params.nPixels-1)/2, ...
        guesses(ii,2)-(params.nPixels-1)/2:guesses(ii,2)+(params.nPixels-1)/2);
    inds = sub2ind(size(img),sDom1(:),sDom2(:));
    
    % select the data
    truImg = img(inds);
    
    % amplitude, offset starting values
    mVals = [nanmax(truImg(:)),nanmin(truImg(:))];
    pStart(4) = mVals(1)-mVals(2);
    pStart(5) = mVals(2);
    
    % fit the data
    [fitPars(ii,:), ~, residual, ~, ~, ~, jacobian] = ...
        lsqcurvefit(fFun,pStart,X(~isnan(truImg(:)),:),truImg(~isnan(truImg(:))),lb,ub,opts);
%     fitPars(ii,:) = lsqcurvefit(fFun,pStart,X(~isnan(truImg(:)),:),truImg(~isnan(truImg(:))),lb,ub,opts);
    % confidence intervals
    conf95(ii,:) = diff(nlparci(fitPars(ii,:), residual, 'jacobian', jacobian), 1, 2);
    
end
% disabled confidence intervals
% conf95 = fitPars;

%% shift origins back to lab frame
fitPars(:,[1,2]) = fitPars(:,1:2) + guesses - params.nPixels;
guesses = guesses - params.nPixels;

%% plot the output
if params.checkVals
    [x,y] = ndgrid(1:imSize(1), 1:imSize(2));
    X=cat(2,x(:),y(:));
    fImg = zeros(imSize);
    
    for ii = 1:nFits
        fImg = fImg + reshape(fFun(fitPars(ii,:),X),imSize);
    end
    fImg = fImg-nanmin(fImg(:));
    
    subplot(2,4,4)
    imshow(fImg,[])
    subplot(2,4,8)
    imshow(img,[])


set(gcf,'NextPlot','add');
axes
h = title(['Frame number ' num2str(params.frameNumber)]);
set(gca,'Visible','off');
set(h,'Visible','on');
end
warning('on','MATLAB:singularMatrix');
end
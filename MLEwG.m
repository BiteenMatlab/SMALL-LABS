% Maximum likelihood fit to a 2D Gaussian with a constant background
%
% From
% K. I. Mortensen, L. S. Churchman, J. A. Spudich, and H. Flyvbjerg, Nat.
% Methods 7, 377 (2010) doi:10.1038/nmeth.1447 
%
%
%   N* 1/(2*pi*s^2) * exp (-( (x-ux).^2+(y-uy).^2 ) / (2*s^2)) + b^2
%
% Input Parameters:
% Required:
% data -- the image of the isolated probe.
% params0 -- The user's initial guess of the parameters to be fit:
%           [ux, uy, s, b, N]
%           ux, uy and s should be specified in nanometers.
%           note: as the fit minimizes only localy it is important for this
%           inital guess to be fairly close to the true value.
% a -- pixel size in nanometers
%
% Optional:
% plot_on -- a binary parameter that determines whether the outcome is
% plotted.
% find_variance -- a binary parameter that determines whether the variance
% should be calculated. If find_variance is 0, then -1 is returned as the
% variance.
% noise_floor -- a constant offset added by the camera before the pixel is
% read out. This is necessary if an accurate spot amplitude is desired.
%
% Output parameter:
% paramsF -- the result of the fit.
%
% Optional:
% variance -- the variance of the localization based on eq.(5) of Mortensen
% et al.

function [paramsF, variance] = MLEwG (data, params0, a, plot_on, ...
                    find_variance, noise_floor)

if nargin < 6
    data = double (data);
else
    data = double (data-noise_floor);
end

if nargin < 4, plot_on = 0; end
if nargin < 5
    find_variance = 0;
    variance  =-1;
end

sizeyx = size (data);
[x, y] = meshgrid (0.5:sizeyx (2)-.5, .5:sizeyx (1)-.5);
x = x*a;
y = y*a;

% The funtion to be minimized is the negative of the log likelihood
% -1 * SOM Eq (49)
datafun = @(params)(sum (sum ((expected (x,y,params,a))))...
                        -sum (sum (data.*log (expected (x,y,params,a)))));
options = optimset ('MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-5,'Display','off');
% fminsearch performs the multivariable minimization
[paramsF,fval,exitflag,output]  = fminsearch (datafun, params0, options);

if find_variance
    b2 = paramsF (4)^2;
    N = paramsF (5);
    sa = paramsF (3);
    F = (@(t)log (t)./(1+(N*a^2*t/(2*pi*sa^2*b2)) ) );
    integral = quadgk (F, 0,1);
    variance = sa^2/N*(1+integral)^-1;
end

if plot_on
    figure (10)
    subplot (1,2,1)
    mesh (data)
    subplot (1,2,2)
    mesh (expected (x,y,paramsF,a))
    figure (11)
    imagesc (data)
    hold on
    plot (paramsF (1)/a+.5, paramsF (2)/a+.5, '*g')
end


end

function p = twoDGauss (x,y,ux,uy,s)
% 2D Gaussian. (SOM Eq. 3)
p = 1/(2*pi*s^2) * exp (-( (x-ux).^2+(y-uy).^2 ) / (2*s^2));

end

function E = expected (x,y,params,a)
% The expected counts per pixel. (SOM Eq. 12)
E = params (5)*a^2*twoDGauss (x,y,params (1),params (2),params (3))+params (4)^2;
end

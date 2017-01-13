function res = bpass(image_array,lnoise,lobject,threshold,lzero)
% 
% NAME:
%               bpass
% PURPOSE:
%               Implements a real-space bandpass filter that suppresses 
%               pixel noise and long-wavelength image variations while 
%               retaining information of a characteristic size.
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               res = bpass( image_array, lnoise, lobject )
% INPUTS:
%               image:  The two-dimensional array to be filtered.
%               lnoise: Characteristic lengthscale of noise in pixels.
%                       Additive noise averaged over this length should
%                       vanish. May assume any positive floating value.
%                       May be set to 0 or false, in which case only the
%                       highpass "background subtraction" operation is 
%                       performed.
%               lobject: (optional) Integer length in pixels somewhat 
%                       larger than a typical object. Can also be set to 
%                       0 or false, in which case only the lowpass 
%                       "blurring" operation defined by lnoise is done,
%                       without the background subtraction defined by
%                       lobject.  Defaults to false.
%               threshold: (optional) By default, after the convolution,
%                       any negative pixels are reset to 0.  Threshold
%                       changes the threshhold for setting pixels to
%                       0.  Positive values may be useful for removing
%                       stray noise or small particles.  Alternatively, can
%                       be set to -Inf so that no threshholding is
%                       performed at all.
%
% OUTPUTS:
%               res:    filtered image.
% PROCEDURE:
%               simple convolution yields spatial bandpass filtering.
% NOTES:
% Performs a bandpass by convolving with an appropriate kernel.  You can
% think of this as a two part process.  First, a lowpassed image is
% produced by convolving the original with a gaussian.  Next, a second
% lowpassed image is produced by convolving the original with a boxcar
% function. By subtracting the boxcar version from the gaussian version, we
% are using the boxcar version to perform a highpass.
% 
% original - lowpassed version of original => highpassed version of the
% original
% 
% Performing a lowpass and a highpass results in a bandpassed image.
% 
% Converts input to double.  Be advised that commands like 'image' display 
% double precision arrays differently from UINT8 arrays.

% MODIFICATION HISTORY:
%               Written by David G. Grier, The University of Chicago, 2/93.
%
%               Greatly revised version DGG 5/95.
%
%               Added /field keyword JCC 12/95.
% 
%               Memory optimizations and fixed normalization, DGG 8/99.
%               Converted to Matlab by D.Blair 4/2004-ish
%
%               Fixed some bugs with conv2 to make sure the edges are
%               removed D.B. 6/05
%
%               Removed inadvertent image shift ERD 6/05
% 
%               Added threshold to output.  Now sets all pixels with
%               negative values equal to zero.  Gets rid of ringing which
%               was destroying sub-pixel accuracy, unless window size in
%               cntrd was picked perfectly.  Now centrd gets sub-pixel
%               accuracy much more robustly ERD 8/24/05
%
%               Refactored for clarity and converted all convolutions to
%               use column vector kernels for speed.  Running on my 
%               macbook, the old version took ~1.3 seconds to do
%               bpass(image_array,1,19) on a 1024 x 1024 image; this
%               version takes roughly half that. JWM 6/07
%
%       This code 'bpass.pro' is copyright 1997, John C. Crocker and 
%       David G. Grier.  It should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.  

if nargin < 3, lobject = false; end
if nargin < 4, threshold = 0; end

normalize = @(x) x/sum(x);

image_array = double(image_array);

if lnoise == 0
  gaussian_kernel = 1;
else      
  gaussian_kernel = normalize(...
    exp(-((-ceil(5*lnoise):ceil(5*lnoise))/(2*lnoise)).^2));
end

if lobject  
  boxcar_kernel = normalize(...
      ones(1,length(-round(lobject):round(lobject))));
end
  
% JWM: Do a 2D convolution with the kernels in two steps each.  It is
% possible to do the convolution in only one step per kernel with 
%
  % gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
  % bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
% 
% but for some reason, this is slow.  The whole operation could be reduced
% to a single step using the associative and distributive properties of
% convolution:
%
  % filtered = conv2(image_array,...
  %   gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
  %   'same');
%
% But this is also comparatively slow (though inexplicably faster than the
% above).  It turns out that convolving with a column vector is faster than
% convolving with a row vector, so instead of transposing the kernel, the
% image is transposed twice.

gconv = conv2(image_array',gaussian_kernel','same');
gconv = conv2(gconv',gaussian_kernel','same');

if lobject
  bconv = conv2(image_array',boxcar_kernel','same');
  bconv = conv2(bconv',boxcar_kernel','same');

  filtered = gconv - bconv;
else
  filtered = gconv;
end

% commented out because why do it in the first place?
% % Zero out the values on the edges to signal that they're not useful.     
% lzero = max(lobject,ceil(5*lnoise));
% lzero=0;

filtered(1:(round(lzero)),:) = 0;
filtered((end - lzero + 1):end,:) = 0;
filtered(:,1:(round(lzero))) = 0;
filtered(:,(end - lzero + 1):end) = 0;

% JWM: I question the value of zeroing out negative pixels.  It's a
% nonlinear operation which could potentially mess up our expectations
% about statistics.  Is there data on 'Now centroid gets subpixel accuracy
% much more robustly'?  To choose which approach to take, uncomment one of
% the following two lines.
% ERD: The negative values shift the peak if the center of the cntrd mask
% is not centered on the particle.

% res = filtered;
filtered(filtered < threshold) = 0;
res = filtered;
function res = saveastiff(data, path, options)
% options.color
%   : true or FALSE
%   : If this is true, third dimension should be 3 and the data is saved as a color image.
% options.comp
%   : 'no', 'lzw', 'jpeg' or 'adobe'.
%     Compression type.
%       'no'    : Uncompressed(Default)
%       'lzw'   : lossless LZW
%       'jpeg'  : lossy JPEG
%       'adobe' : lossless Adobe-style
% options.message
%   : true or FALSE.
%     If this is false, all messages are skipped. 
% options.append
%   : true or FALSE
%     If path is exist, the data is appended to an existing file.
%     If path is not exist, this options is ignored.
% options.overwrite
%   : true or FALSE
%     Overwrite to an existing file.
% options.big 
%   : true or FALSE, 
%     Use 64 bit addressing and allows for files > 4GB
% 
% Defalut value of 'options' is
%     options.color     = false;
%     options.comp      = 'no';
%     options.message   = false;
%     options.append    = false;
%     options.overwrite = false;
%     options.big       = false;
% 
% res : Return value. It is 0 when the function is finished with no error.
%       If an error is occured in the function, it will have a positive
%       number (error code).
%
% Copyright (c) 2012, YoonOh Tak
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Gwangju Institute of Science and Technology (GIST), Republic of Korea nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

tStart = tic;
errcode = 0;
try
if isempty(data)
    errcode = 1; assert(false);
end
if isreal(data) == false
    errcode = 2; assert(false);
end
if nargin < 3
    options.color = false;
    options.comp = 'no';
    options.message = false;
    options.append = false;
    options.overwrite = false;
end
if ~isfield(options, 'message')
    options.message = false;
end
if ~isfield(options, 'append')
    options.append = false;
end
if ~isfield(options, 'comp')
    options.comp = 'no';
end
if ~isfield(options, 'color')
    options.color = false;
end
if ~isfield(options, 'overwrite')
    options.overwrite = false;
end
if (options.color == false && ndims(data) > 3) ...
    || (options.color == true && ndims(data) > 4)
    errcode = 3; assert(false);
end
if isfield(options, 'big') == 0 
    options.big = false; 
end

if ~options.color
    if ndims(data) >= 4, errcode = 3; assert(false); end;
    [height, width, depth] = size(data);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
else
    if ndims(data) >= 5, errcode = 3; assert(false); end;
    [height, width, rgb, depth] = size(data);
    if rgb ~= 3, errcode = 4; assert(false); end;
    if isa(data, 'uint8'), data = uint8(data); end;
    tagstruct.Photometric = Tiff.Photometric.RGB;
end

tagstruct.ImageLength = height;
tagstruct.ImageWidth = width;
tagstruct.SamplesPerPixel = (options.color+1)*(options.color+2)/2;
tagstruct.RowsPerStrip = height;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% Compresstion type : http://en.wikipedia.org/wiki/Tagged_Image_File_Format
switch lower(options.comp)
    case 'no'
        tagstruct.Compression = Tiff.Compression.None;
    case 'lzw'
        tagstruct.Compression = Tiff.Compression.LZW;
    case 'jpeg'
        tagstruct.Compression = Tiff.Compression.JPEG;
    case 'adobe'
        tagstruct.Compression = Tiff.Compression.AdobeDeflate;
    otherwise
        tagstruct.Compression = Tiff.Compression.None;
end

switch class(data)
    case {'uint8', 'uint16', 'uint32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'int16', 'int32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        if options.color == 3
            errcode = 5; assert(false);
        end
    case {'single', 'double', 'uint64', 'int64'}
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        errcode = 6; assert(false);
end

switch class(data)
    case {'uint8', 'int8'}
        tagstruct.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tagstruct.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tagstruct.BitsPerSample = 32;
    case {'single'}
        tagstruct.BitsPerSample = 32;
    case {'double', 'uint64', 'int64'}
        tagstruct.BitsPerSample = 64;
        data = double(data);
    otherwise
        errcode = 6; assert(false);
end

if exist(path, 'file') && options.append == false
    if ~options.overwrite
        errcode = 7; assert(false);
    end
end

path_parent = pwd;
[pathstr, fname, fext] = fileparts(path);
if ~isempty(pathstr)
    if ~exist(pathstr, 'dir')
        mkdir(pathstr);
    end
    cd(pathstr);
end

if options.append == false
    s=whos('data'); 
    if s.bytes > 2^32-1 || options.big 
        tfile = Tiff([fname, fext], 'w8');
    else 
        tfile = Tiff([fname, fext], 'w');
    end
    for d = 1:depth
        tfile.setTag(tagstruct);
        if ~options.color
            tfile.write(data(:, :, d));
        else
            tfile.write(data(:, :, :, d));
        end
        if d ~= depth
            tfile.writeDirectory();
        end
    end
else
    if ~exist([fname, fext], 'file')
        s=whos('data'); 
        if s.bytes > 2^32-1 || options.big 
            tfile = Tiff([fname, fext], 'w8');
        else 
            tfile = Tiff([fname, fext], 'w');
        end
    else
        tfile = Tiff([fname, fext], 'r+');
        while ~tfile.lastDirectory();
            tfile.nextDirectory();
        end
        tfile.writeDirectory();
    end
    
    for d = 1:depth
        tfile.setTag(tagstruct);
        if ~options.color
            tfile.write(data(:, :, d));
        else
            tfile.write(data(:, :, :, d));
        end
        if d ~= depth
            tfile.writeDirectory();
        end
    end
end
tfile.close();
if exist('path_parent', 'var'), cd(path_parent); end

tElapsed = toc(tStart);
if options.message
    display(sprintf('File saved successfully. Elapsed time : %.3f s.', tElapsed));
end

catch exception
    if exist('tfile', 'var'), tfile.close(); end
    if options.message
        switch errcode
            case 1
                error '''data'' is empty.';
            case 2
                error 'Data type of ''data'' should be real numbers';
            case 3
                error 'Data dimension is too large.';
            case 4
                error 'Third dimesion (color depth) should be 3.';
            case 5
                error 'RGB color image cannot have int8, int16 or int32 format.';
            case 6
                error 'Unsupported data type.';
            case 7
                error 'File already exists.';
            otherwise
                rethrow(exception);
        end
    end
    if exist('path_parent', 'var'), cd(path_parent); end
end
res = errcode;
end
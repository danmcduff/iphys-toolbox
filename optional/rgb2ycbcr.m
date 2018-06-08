% Copyright (C) 2011 Santiago Reyes Gonz??lez <mailalcuete@gmail.com>
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% -*- texinfo -*-
% @deftypefn {Function File} @var{B} =  rgb2ycbcr(@var{A})
% Perform colour convertion from RGB to YCbCr on a given image.
%
% The image @var{A} must be a NxMx3 image. The conversion
% The convertion changes the image from the RGB color model to YCbCr e.g.
% @example
% imOut = rgb2ycbcr(imIn);
% @end example
% Currently this function only works with @samp{uint8} and will always return
% an @samp{uint8} matrix.
% @seealso{cmunique}
% @end deftypefn

function im_out =  rgb2ycbcr(im)
  if (nargin ~= 1)
    print_usage;
  elseif (length(size(im)) ~= 3 || size(im,3) ~= 3)
    error("image must be NxMx3");
  end

  im            = im2double(im);
  im_out(:,:,1) = uint8(floor(77*im(:,:,1)    + 150*im(:,:,2) + 29*im(:,:,3)));
  im_out(:,:,2) = uint8(floor(((-44*im(:,:,1) - 87*im(:,:,2)  + 131*im(:,:,3))/256 + 0.5)*256));
  im_out(:,:,3) = uint8(floor(((131*im(:,:,1) - 110*im(:,:,2) - 21*im(:,:,3))/256 + 0.5)*256));

end
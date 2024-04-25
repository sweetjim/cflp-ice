function out = noisedither(img,ntype)
%   NOISEDITHER(INPICT,{TYPE})
%       Apply simple noise threshold dither to an I/RGB image
%
%   INPICT is a 2-D intensity image
%       if fed an RGB image, its luma channel will be extracted
%
%   TYPE is the noise type (default 'blue')
%       'white' for noise with a flat power spectrum
%       'blue' for noise with a power density roughly proportional to f
%
%   Output class is logical
%
% Webdocs: http://mimtdocs.rf.gd/manual/html/noisedither.html
% See also: dither, zfdither, orddither, arborddither, linedither

if ~exist('ntype','var')
	ntype = 'blue';
end

if size(img,3) == 3
	img = mono(img,'y');
end

img = imcast(img,'double');
s = size(img);

switch lower(ntype)
	case 'white'
		mask = rand(s(1:2));
	case 'blue'
		mask = rand(s(1:2));
		
		% this replaces imgaussfilt()
		sg = 2; R = 4;
		[xx yy] = meshgrid(-R:R,-R:R);
		fk = exp(-(xx.^2 + yy.^2)/(2*sg^2));
		fk = fk/sum(sum(fk));
		
		maskdiff = simnorm(mask-imfilterFB(mask,fk,'replicate'));
		mask = adapthisteqFB(maskdiff);
	otherwise
		error('NOISEDITHER: unknown noise type')
end

out = mask <= img;

end


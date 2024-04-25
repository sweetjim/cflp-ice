function [colorchans alphachans] = splitalpha(inpict)
%  [COLOR ALPHA] = SPLITALPHA(INPICT)
%    Split an image into its color and alpha channels.
%
%  INPICT is an I/IA/RGB/RGBA/RGBAAA image of any standard image class.
%    Multiframe images are supported.
%
%  Output class is inherited from the input.
%
%  Webdocs: http://mimtdocs.rf.gd/manual/html/splitalpha.html
%  See also: splitchans, chancount

[nc,~] = chancount(inpict);

alphachans = inpict(:,:,nc+1:end,:);
colorchans = inpict(:,:,1:nc,:);

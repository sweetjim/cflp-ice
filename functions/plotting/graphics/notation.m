function string = notation(value)
add = 0;
if value<0
   value    = abs(value);
   add      = 1;
end
exponent    = floor(log10(value));
number      = value*10^(-exponent);

string      = sprintf('%.1f$%s10^{%s}$',number,'\times',num2str(exponent));

if add
   string = strcat('-',string) ;
end
end


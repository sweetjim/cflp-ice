function input = nonempty(input)
%NONEMPTY returns the nonempty elements of the input array (e.g. empty
%cells, character-less strings, nan doubles)
if isa(input,'cell')
    input = input(~cellfun(@isempty,input));
elseif isa(input,'string')
    input = input(input~='');
elseif isa(input,'double')
    input = input(~isnan(input));
end

end


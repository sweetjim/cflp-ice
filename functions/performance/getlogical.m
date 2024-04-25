function value = getlogical(value)
%GETLOGICAL accepts logicals values and character and returns values

if isa(value,'char')||isa(value,'string')
    try
        value = eval(lower(value));
    catch
        error('Value is not of logical format')
    end
else
    value = logical(value);
end
end


function createFunctionSignatures(fcnhandle,argtable)
%CREATEFUNCTIONSIGNATURES 
s0      = struct('_schemaVersion','1.0.0');
argtable= []; 



    function j = encode(name,kind,type,choice)
        if numel(choice)>0
            type = struct(type,choice);
        end
        s = struct('name',name,'kind',kind,'type',type);
        j = jsonencode(s);
        %{"name":"input1",  "kind":"required", "type":["numeric"]}

    end
end
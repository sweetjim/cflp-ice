function str = symbol(label)
%SYMBOL is a collection of special character symbols

%% Greek letters
GLlower = {'alpha';'beta';'gamma';'delta';'epsilon';'zeta';'eta';'theta';'iota';'kappa';'lambda';'mu';'nu';'xi';'omicron';'pi';'rho';'sigma';'sigma2';'tau';'upsilon';'phi';'chi';'psi';'omega'};
GLloweridx = 945-1;
GLupper = cellfun(@capitalize,GLlower,'UniformOutput',false);
GLupperidx = 913-1;

%% Super+subscripts

if any(contains(label,{'^','_'}))
 str = sscript(label,label(contains(label,{'^','_'})));
 return
end

%%
switch label

    case {'permil','perthousand'}
        str = char(8240);
    case {'degC','degreecelsius','degreeCelsius'}
        str = [char(176) 'C'];
    case {'degClatex'}
        str = '$^\circ$C';
    case {'pm'}
        str = char(177);
    case GLlower
        char(GLloweridx+find(strcmp(GLlower,label)))
    case GLupper
        char(GLupperidx+find(strcmp(GLupper,label)))
end

%% Nested functions
    function str = sscript(label,type)
        operators = {'-','+','='};
        switch type
            case '^'
                typeval=7;
            case '_'
                typeval=8;
        end
        % Superscript
        str0 = extractAfter(label,type);
        if contains(label,operators)
            for i=1:numel(operators)
                switch operators{i}
                    case '-'
                        uc = sprintf('20%iB',typeval);
                    case '+'
                        uc = sprintf('20%iA',typeval);
                    case '='
                        uc = sprintf('20%iC',typeval);
                end
                if any(label==operators{i})
                    labels=split(extractAfter(label,type),operators);
                    tmpstr=[{''};cellfun(@(x) symbol([type x]),labels(2:end),'UniformOutput',false)];
                    str=join(tmpstr,char(hex2dec(uc)));
                    str=str{1};
                    idx=extractAfter(label,type)==operators{i};
                    str0(idx)=str(idx);
                    if i==1
                        str0(~idx)=str(~idx);
                    end
                else
                    continue
                end
            end
            str=str0;
            return
        end

        if contains(label,'.')
            labels = split(label,'.');
            labels{2} = [type labels{2}];
            str=[symbol(labels{1}) char(hex2dec('22c5')) symbol(labels{2})];
            return
        end

        val = str2double(extractAfter(label,type));
        if val>9
            val = num2str(val)';
            str=arrayfun(@(x) symbol([type x]),val)';
            return
        end
        lead = sprintf('20%i',typeval);
        if strcmp(type,'^')
            lead = '00b';
            switch val
                case 1
                    val = 9;
                case {[2],[3]}
                otherwise
                    lead = sprintf('20%i',typeval);
            end
        end
        hexval = hex2dec(sprintf('%s%i',lead,val));
        str = char(hexval);
    end
end
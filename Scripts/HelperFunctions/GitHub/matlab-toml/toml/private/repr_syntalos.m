% REPR TOML representation of a MATLAB object
%
%   REPR(obj) returns the TOML representation of `obj`.

function str = repr_syntalos(obj, parent, grandparent)

if ispc
    newline = sprintf('\r\n');
else
    newline = sprintf('\n');
end

switch class(obj)
    
    % strings
    case 'char'
        % A Syntalos error occured with a strangely formatted framerate. To
        % keep the 12 digit precision I converted it to a string first,
        % which should then however be displayed as a number. The format
        % version (e.g. '1'), should however be displayed as a string.
        if ~isnan(str2double(obj))
            obj_double = str2double(obj);
            if numel(num2str(abs(obj_double))) ~= 1
                str = lower(obj);
            else
                str = ['''', obj, ''''];
            end
        elseif isrow(obj) || isempty(obj)
            str = ['''', obj, ''''];
        else
            str = repr_syntalos(reshape(cellstr(obj), 1, []));
        end
        
        % Booleans
    case 'logical'
        reprs = {'false', 'true'};
        str = reprs{obj + 1};
        
        % numbers
    case 'double'
        if numel(obj) == 1
            str = lower(num2str(obj));
        else
            cel = arrayfun(@repr_syntalos, obj, 'uniformoutput', false);
            str = ['[', strjoin(cel, ', '), ']'];
        end
        
        % cell arrays
    case 'cell'
        if all(cellfun(@isstruct, obj))
            fmtter = @(a) sprintf('\n    [[%s.%s]]%s%s', grandparent, parent, newline, repr_syntalos(a,parent));
            cel_str = cellfun(fmtter, obj, 'uniformoutput', false);
            str = strjoin(cel_str);
        else
            cel_mod = cellfun(@repr_syntalos, obj, 'uniformoutput', false);
            str = ['[', strjoin(cel_mod, ', '), ']'];
        end
        
        % structures
    case 'struct'
        fn = fieldnames(obj);
        vals = struct2cell(obj);
        str = '';
        for indx = 1:numel(vals)
            new_parent = fn{indx};
            if isstruct(vals{indx})
                if nargin > 1
                    fmt_str = ['%1$s\n[', parent, '.%2$s]%4$s%3$s'];
                    new_parent = [parent, '.', fn{indx}];
                else
                    fmt_str = '%1$s\n[%2$s]%4$s%3$s';
                end
            elseif iscell(vals{indx}) && all(cellfun(@isstruct, vals{indx}))
                fmt_str = '%1$s%3$s%4$s';
            elseif exist("parent","var") && strcmp(parent, 'parts')
                fmt_str = '%1$s    %2$s = %3$s%4$s';
            else
                fmt_str = '%1$s%2$s = %3$s%4$s';
            end
            if exist("parent","var")
                str = sprintf(fmt_str, str, fn{indx}, repr_syntalos(vals{indx}, new_parent, parent), newline);
            else
                str = sprintf(fmt_str, str, fn{indx}, repr_syntalos(vals{indx}, new_parent), newline);
            end
        end
        
        % datetime objects
    case 'datetime'
        obj.Format = 'yyyy-MM-dd''T''HH:mm:ss';
        if ~isempty(obj.TimeZone)
            obj.Format = [obj.Format, 'XXX'];
        end
        str = char(obj);
        
        % unrecognized type
    otherwise
        error('toml:NonEncodableType', ...
            'Cannot encode type as TOML: %s', class(obj))
end
end
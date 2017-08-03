function values = nnConditions(params)
%% Make the values used for the Conditions file. 
%
%   values = nnConditions(params)
%
% Input
%   params - a struct with many fields
%
% Return
%   values - a cell array that can be passed to 
%
%        rtbWriteConditionsFile(conditionsFile,names,values);
%
% The fields should be the names of the parameters in the conditions file. 
%
% We use num2cell or mat2str to turn the numerical values (numbers or vectors)
% into cells containing strings. 
%
% If the input is a matrix, then we assume the conditions are the rows and the
% vector for this condition is across the column.
%
% We will use ndgrid to generate all of the permutations of conditions.
% 
% HB/BW SCIEN Team, 2017

% Get the field names
names = fieldnames(params);
nNames = length(names);

% Turn the parameters into cell arrays.  The number of rows determines the
% number of dimensions in this cell array.
for ii=1:nNames
    thisVal = params.(names{ii});
    if isnumeric(thisVal), thisVal = double(thisVal); end
    thisCell = cell(size(thisVal,1),1);    % Number of rows
    switch class(thisVal)
        case {'double'}
            for jj=1:size(thisVal,1)
                thisCell{jj} = mat2str(thisVal(jj,:));
            end
        case {'char'}
            for jj=1:size(thisVal,1)
                thisCell{jj} = thisVal(jj,:);
            end
        case {'cell'}
            for jj=1:size(thisVal,1)
                thisCell{jj} = thisVal{jj,:};
            end
        otherwise
            error('Unknown value conversion class %s\n',class(thisVal));
    end
    params.(names{ii}) = thisCell;
end

% Count how many conditions we expect to have
nConds = 1;
for jj=1:nNames
    nConds = nConds*size(params.(names{jj}),1);
end
values = cell(nConds,nNames);
fprintf('Initially: %d conditions by %d parameters\n',nConds,nNames);

% Use ndgrid to make the cell array of values{conditions,parameter}. This will
% be written out by rtbWriteConditionsFile
gv = cell(nNames,1);
for ii=1:nNames, gv{ii} = 1:size(params.(names{ii}),1); end

ostr = '[X1'; for ii=2:nNames, ostr = sprintf('%s,X%d',ostr,ii); end
ostr = [ostr,']='];

str = 'ndgrid(gv{1}';
for ii=2:nNames, str = sprintf('%s,gv{%d}',str,ii); end
str = [str,');'];

cmd = [ostr,str];
eval(cmd);

% Combine the Xi matrices into an nConds by nNames matrix
str = 'X = [X1(:)';
for ii=2:nNames, str = sprintf('%s,X%d(:)',str,ii); end
str = [str,'];'];
eval(str);

% Now fill in the nConds x nNames cell array, values, using the ndgrid indices.
for ii=1:nConds
    for jj=1:nNames
        values{ii,jj} = params.(names{jj}){X(ii,jj)};
    end
end

end


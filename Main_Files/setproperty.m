function [OPT OPTset] = setproperty(OPT, vargin)

OPTfields = fieldnames(OPT);


for i = 1: length(OPTfields)
    OPTset.(OPTfields{i}) = false;
end

if ~isempty(vargin)
    
    for i = 1: length(vargin)
        VarNames = find(strcmp(vargin{i},OPTfields));
        
        if ~isempty(VarNames)
          OPTset.(OPTfields{VarNames}) = true; 
          OPT.(OPTfields{VarNames}) = vargin{i+1}; 
        end
    end
    
end

return
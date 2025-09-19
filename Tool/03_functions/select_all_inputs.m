inputNames = {};

for idx = 1:length(l_blocks)
    ss_block = l_blocks{idx};
    inputNames = cat(2, inputNames, ss_block.InputName');
end

inputNames = sort(unique(inputNames)); % display BOTH internal and external outputs
input = inputNames;
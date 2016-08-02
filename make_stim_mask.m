% make_stim_mask.m
%
% Creates a mask of a round stimulus for each stimulus (trial), which is an
% element of stimX, stimY, each is stim_size (either a single value or row
% the same length as stimX and stimY), evaluated at [xx yy], which is
% reshaped into columns, so return is n_stimpts x length(stimX)




function stim_mask = make_stim_mask(stimX,stimY,stim_size,xx,yy)


if length(stim_size)==1
    stim_size = stim_size * ones(length(stimX),1);
end

stim_mask = zeros(length(xx),length(stimX));

for ii = 1:length(stimX)
    
    rr = sqrt((xx-stimX(ii)).^2+(yy-stimY(ii)).^2);
    stim_mask(rr <= stim_size(ii),ii) = 1;
    
end



return
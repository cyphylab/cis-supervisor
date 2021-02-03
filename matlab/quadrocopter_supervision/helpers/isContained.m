function [isContained, idx] = isContained(x0, DSets, check_idx)
idx = zeros(length(DSets),1);
isContained = false;

if (isempty(check_idx))
    for k = 1:length(DSets)
        if (DSets{k}.CIS.contains(x0))
            isContained = true;
            idx(k) = 1;
        end
    end
else
    if (DSets{check_idx}.CIS.contains(x0))
        isContained = true;
        idx(check_idx) = 1;
    end
end
function flag = isInCis(CISs, x0)
flag = 0;
for k = 1:length(CISs)
    if (CISs{k}.CIS.contains(x0))
        flag = 1;
    end
end

end
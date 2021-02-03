% CSV export script

for i = 1 : length(DSets)
nameA = sprintf("cis%d_A.csv", i-1);
nameB = sprintf("cis%d_b.csv", i-1);

writematrix(DSets{i}.CIS.A, nameA);
writematrix(DSets{i}.CIS.b, nameB);
end

writematrix(mdl_dt.Ad, "Ad.csv");
writematrix(mdl_dt.Bd, "Bd.csv");
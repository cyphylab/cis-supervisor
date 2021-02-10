% CSV export script

for i = 1 : length(DSets)
nameA = sprintf("../../data/cis%d_A.csv", i-1);
nameB = sprintf("../../data/cis%d_b.csv", i-1);
writematrix(DSets{i}.CIS.A, nameA);
writematrix(DSets{i}.CIS.b, nameB);

nameA = sprintf("../../data/rcis%d_A.csv", i-1);
nameB = sprintf("../../data/rcis%d_b.csv", i-1);
writematrix(DSets{i}.RCIS.A, nameA);
writematrix(DSets{i}.RCIS.b, nameB);
end

writematrix(mdl_dt.Ad, "../../data/Ad.csv");
writematrix(mdl_dt.Bd, "../../data/Bd.csv");
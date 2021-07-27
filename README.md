# isocorr_peptide
pepcorr.m is the central function, which takes peptide sequence and measured peptide MID, outputs the simulated MID(due to natural abundance only) and
natural abundance corrected MID from optimization using fmincon. fitness function is least square of the error between simulated (natural abundance (*) corrected) and measured
It's optional to add constrains to the upper bounds (0-1) of the corrected MID.

[MID_sim,MID_corr,fval]=pepcorr(pep,MID_measure,constrain)

main.m does the batch processes for all peptides from the input csv sheet.

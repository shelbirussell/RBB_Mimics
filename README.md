# RBB_Mimics
Endosymbiont-host protein mimic detection via reciprocal best BLAST 

Run the following scripts in the given order:

1. Make and submit reciprocal blastp jobs for each sequence pair in a directory.
make_all_by_all_blastp.pl

2. Filter each blast hit file, retaining best hits >100 bp and >30% identity.
blast_besthit_filter4orthology.pl

3. Find reciprocal best blast hits 
find_reciprocal_all_by_all.pl

4. Extract sequences for each rbb hit
make_orthology_multifasta.pl



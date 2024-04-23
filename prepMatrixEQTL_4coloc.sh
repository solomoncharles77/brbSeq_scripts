
#!/bin/bash


Rscript brbSeq_scripts/prepMatrixEQTL_4coloc.R -q tnfa_pilot_BRBseq_cisEQTL_2pc_3pf_Results/tnfa_pilot_BRBseq_cisEQTL_2pc_3pf.csv -f genoFiles/tnfa_pilot_BRBseq.frq

Rscript brbSeq_scripts/prepMatrixEQTL_4coloc.R -q untreated_pilot_BRBseq_cisEQTL_2pc_3pf_Results/untreated_pilot_BRBseq_cisEQTL_2pc_3pf.csv -f genoFiles/untreated_pilot_BRBseq.frq

Rscript brbSeq_scripts/prepMatrixEQTL_4coloc.R -q huvec50_cisEQTL_2pc_3pf_Results/huvec50_cisEQTL_2pc_3pf.csv -f genoFiles/huvec50.frq

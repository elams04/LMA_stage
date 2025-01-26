#!/bin/bash

jid1=$(sbatch prepro.sh)
jid1=${jid1##* }  # Extract the job ID from the submission output

jid2=$(sbatch --dependency=afterok:$jid1 rf.sh)
jid2=${jid2##* }

jid3=$(sbatch --dependency=afterok:$jid2 run.sh)
jid3=${jid3##* }

sbatch --dependency=afterok:$jid3 extractcap.sh


#!/bin/bash

doRun=1

leap3="/home/codex/chege/software/pipelines/twoleap/main.nf"
wpath="/home/codex/chege/software/pipelines/twoleap/test"
hosts="${wpath}/pssh_hosts_list.txt"
datapath="/data/users/lofareor/chege/lptest/L412388"
logsdir="${wpath}/logs"
label="${datapath}/L412388_SAP000_SB*_uv_002_10mins.MS"
mslist=/home/codex/chege/software/pipelines/twoleap/test/mslist_002.txt
nodes="101,102,103,105,114"
split_msout="L412388_SAP000_SB000_uv_002"
ditodd_avg_msout="L412388_SAP000_SB000_uv_003"

if [ ${doRun} == 1 ]; then
    nextflow run ${leap3} --data.path ${datapath} --data.ms "${label}" --ddecal.bp.solint 300 --out.logs ${logsdir} --data.mslist ${mslist} --data.nodes ${nodes} --split.ntimes 120 --split.msout ${split_msout} --ddecal.di.solint 15 --average.ditodd.msout ${ditodd_avg_msout}
fi
#  --data.hosts ${hosts} 
# for node in $(cat hosts.txt); do ssh "$host" "$command" >"output.$host"; done

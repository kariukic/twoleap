#!/bin/bash
doMakeTestDataset=0
doProcess=1

sif=~/software/containers/leor_tools_200525.sif

obsID='L253456'
frompath="/data/users/lofareor/chege/3C196/L253456/zbin2/try2/"
raw_glob="${frompath}/L253456_SAP000_SB*_uv.MS.5ch4s.dppp"
processingDir="/data/users/lofareor/chege/twoleap/${obsID}"

if [ ${doMakeTestDataset} == 1 ]; then

    starttimeslot=2400
    ntimes=120

    for node in {301..323}; do
        outdir="/net/node${node}/${processingDir}"
        mkdir -p "${outdir}"
        for MSi in /net/node${node}/${raw_glob}; do
            msout=$(basename ${MSi})
            msout="${outdir}/${msout}"
            com="singularity exec --bind /data:/data,/net:/net -e -C ${sif} DP3 msin=${MSi} steps=[] msin.starttimeslot=${starttimeslot} msin.ntimes=${ntimes} msin.datacolumn=DATA msout=${msout}"
            echo "${com}"
            ${com}
        done
    done

fi

nodes="301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323"
splitNodes="301,302,303,304,305"

## Paths to directories
project="3C196"
projectDir="/home/codex/chege/projects/${project}/"

workDir="/home/codex/chege/software/pipelines/twoleap/test/"
logsDir="${workDir}/logs/"

##Sky models
skyModelDI="${projectDir}/skymodel/dppp/L253456_SAP000-3.9deg-DI.txt"
# skyModelDD="${projectDir}/skymodel/dppp/L253456_SAP000-3.9deg-clustered-CasA_TauA_lowres.txt"
skyModelDD="/home/codex/chege/software/pipelines/twoleap/test/ddmodel_3clusters.txt"

##Parsets
parsetDISmoothCal="${projectDir}/parsets/DI_smooth_cal.parset"
parsetDISmoothCalApply="${projectDir}/parsets/applycal_DI_smooth.parset"
parsetDIBandpassCal="${projectDir}/parsets/DI_bandpass_cal.parset"
parsetDIBandpassCalApply="${projectDir}/parsets/applycal_DI_bandpass_cal.parset"
parsetDDSmoothCal="${projectDir}/parsets/DD_smooth_cal_minimal.parset"
parsetDDSubtract="${projectDir}/parsets/DD_subtract.parset"

if [ ${doProcess} == 1 ]; then

    raw_glob="${processingDir}/${obsID}_SAP000_SB???_uv.MS.5ch4s.dppp"

    nextflow run ~/software/pipelines/twoleap/main.nf \
        --data.nodes ${nodes} \
        --data.raw_ms_glob "${raw_glob}" \
        --data.path ${processingDir} \
        --ssplit.nodes ${splitNodes} \
        --average.lta_to_di.freqstep 1 \
        --ssplit.di.ntimes 2 \
        --ssplit.di.mses_per_node 3 \
        --ssplit.dd.ntimes 2 \
        --ssplit.dd.mses_per_node 3 \
        --ssplit.di.ms_prefix "fullband_di_smooth_data.MS" \
        --ddecal.di.sourcedb ${skyModelDI} \
        --ddecal.di.parset ${parsetDISmoothCal} \
        --ddecal.di.apply.parset ${parsetDISmoothCalApply} \
        --ddecal.bp.sourcedb ${skyModelDI} \
        --ddecal.bp.parset ${parsetDIBandpassCal} \
        --ddecal.bp.apply.parset ${parsetDIBandpassCalApply} \
        --ddecal.dd.solint 1 \
        --ddecal.dd.parset ${parsetDDSmoothCal} \
        --ddecal.dd.sourcedb ${skyModelDD} \
        --ddecal.dd.subtract.parset ${parsetDDSubtract} \
        --out.logs ${logsDir} \
        --start_at 'postdd'
    # --stop_at 'preprocess'
fi

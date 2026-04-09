#!/usr/bin/env nextflow

include {
    H5ParmCollect ;
    AOqualityCombine ;
    WScleanImage ;
    MergeChansSplitTime ;
    ReadTxtLinesandAppend ;
    WriteMSlist ;
    writeHosts ;
    parseNodes ;
    makeDirectory ;
    readTxtIntoString ;
    readTxtAndAppendString ;
    GetTimeChunksPerSubband ;
    GetData
} from "./processes.nf"


include {
    AddRevision ;
    RunPSPIPE
} from "./makeps.nf"

// Helper function to check if step should run
def shouldRun(step) {

    // Define execution order
    def steps = [
        'untar_ms',
        'preprocess',
        'make_di_ms',
        'di_cal',
        'split_freq',
        'concat_time',
        'bp_cal',
        'average',
        'make_dd_ms',
        'dd_cal',
        'postdd',
        'pspec',
        'image',
    ]

    if (!steps.contains(params.start_at)) {
        error("Invalid starting step.\n start_at: ${params.start_at}.\n Valid steps: ${steps}")
    }

    if (!steps.contains(params.stop_at)) {
        error("Invalid stopping step.\n stop_at: ${params.stop_at}.\n Valid steps: ${steps}.")
    }

    def startIdx: Integer = steps.indexOf(params.start_at)
    def stopIdx: Integer = steps.indexOf(params.stop_at)
    def stepIdx: Integer = steps.indexOf(step)

    if (startIdx > stopIdx) {
        error("start_at ('${params.start_at}') must come before stop_at ('${params.stop_at}') in pipeline sequence:\n ${steps}")
    }

    stepIdx >= startIdx && stepIdx <= stopIdx
}


workflow {
    ext_ch = shouldRun('untar_ms') ? EXTRACT(true) : true

    pre_process_ch = shouldRun('preprocess') ? PreProcess(ext_ch) : true

    split_ch1 = shouldRun('make_di_ms') ? MakeFullBandDITimeChunks(pre_process_ch) : true

    di_ch = shouldRun('di_cal') ? RunDISmooth(split_ch1) : true

    sbs_ch = shouldRun('split_freq') ? SplitMSChannels(di_ch) : true

    sbs_ch2 = shouldRun('concat_time') ? ConcatTimeChunks(sbs_ch) : true

    bp_ch = shouldRun('bp_cal') ? RunDIBandpass(sbs_ch2) : true

    avg_ch = shouldRun('average') ? Average(bp_ch) : true

    split_ch2 = shouldRun('make_dd_ms') ? MakeFullBandDDTimeChunks(avg_ch) : true

    dd_ch = shouldRun('dd_cal') ? Run_DD(split_ch2) : true

    pdd_ch = shouldRun('postdd') ? Run_PostDD(dd_ch) : true

    ps_ch = shouldRun('pspec') ? PowerSpectrum(pdd_ch) : true

    shouldRun('image') ? Image(ps_ch) : null
}


workflow ExtractDATA {
    EXTRACT(true)
}


process InitParams {
    debug true
    publishDir params.out.logs, mode: 'copy'

    input:
    val ready
    val stage

    output:
    path "${stage}_params.json", emit: params_file
    val true, emit: params_standby

    script:
    makeDirectory(params.out.logs)

    def tasks_after_split = ["DI", "DD", "PDD", "SB"]
    //, "AVG"

    if (tasks_after_split.contains(stage)) {
        nodes_list = parseNodes(params.ssplit.nodes)
    }
    else {
        nodes_list = parseNodes(params.data.nodes)
    }

    writeHosts(nodes_list, params.data.hosts)

    """
        echo '${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))}' > ${stage}_params.json
        """
}


process Distribute {
    debug true

    input:
    val ready
    val ch_in
    val entry
    val params_file

    output:
    val true

    script:
    """
        mkdir -p ${params.out.logs}
        pssh -v -i -h ${launchDir}/${params.data.hosts} -t 0 -x "cd ${params.data.path}; bash" /home/codex/chege/software/nextflow run ${params.stagelib} --stage ${entry} --ch_in ${ch_in} -params-file ${params_file}  -with-singularity ${params.image} --singularity_bind_path ${params.binds} > ${params.out.logs}/${entry}.log 2>&1
        """
}


workflow PreProcess {
    take:
    ready

    main:
    stage_ch = channel.of('FCAB')

    stage_params_ch = InitParams(ready, stage_ch)

    cal_ch = Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    nodesList = params.data.nodes?.split(',') as List
    nodes_ch = channel.fromList(nodesList).collect { it }

    mslist_ch = WriteMSlist(cal_ch, nodes_ch, "${params.data.path}/*.noInter.flagged.di_averaged", params.data.di_mslist)

    AOqualityCombine(cal_ch, mslist_ch.per_line_mslist, "raw_aoq_stats_averaged")

    emit:

    AOqualityCombine.out.done
}


workflow MakeFullBandDITimeChunks {
    take:
    ready

    main:

    def mslist_file = file(params.data.di_mslist)
    if (mslist_file.exists()) {
        mses = Channel.value(mslist_file).splitText().map { file(it.trim()) }.collect().map { it.join(' ') }
    }
    else {
        nodesList = params.data.nodes?.split(',') as List
        nodes_ch = channel.fromList(nodesList).collect { it }
        mses = WriteMSlist(ready, nodes_ch, "${params.data.path}/*noInter.flagged.di_averaged", params.data.di_mslist).per_line_mslist.splitText().map { file(it.trim()) }.collect().map { it.join(' ') }
    }

    nodesList = params.ssplit.nodes?.split(',') as List
    nodes_ch = channel.fromList(nodesList).collect { it }

    output_mslist = file(params.out.logs).resolve(params.data.di_mslist)

    MergeChansSplitTime(ready, mses, nodes_ch, params.ssplit.di.ntimes, params.ddecal.di.incol, params.data.di_ms_glob.replace('_T???', '.MS'), params.ssplit.di.mses_per_node, output_mslist)

    emit:

    MergeChansSplitTime.out
}


workflow RunDISmooth {
    take:
    ready

    main:

    stage_ch = channel.of('DI')

    stage_params_ch = InitParams(ready, stage_ch)

    cal_ch = Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    mses_sols_ch = ReadTxtLinesandAppend(cal_ch, params.out.logs, params.data.di_mslist, "/${params.ddecal.di.sols}")

    sols_collect_ch = H5ParmCollect(cal_ch, mses_sols_ch.list_postfix_str, "di_smooth_solutions")

    aoq_comb_ch = AOqualityCombine(sols_collect_ch.done, "${params.out.logs}/${params.data.di_mslist}", "di_smooth_aoq_stats")
    //

    ps_dir = "${params.data.path}/${params.out.results}/${params.pspipe.dir}/di_smooth"
    ps_logs_dir = "${params.out.logs}/power_spectrum/di"

    rev_ch = AddRevision(aoq_comb_ch.done, params.pspipe.obsid, params.ddecal.di.outcol, params.data.path, ps_dir, params.pspipe.node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, params.pspipe.time_start_index, params.pspipe.time_end_index)

    RunPSPIPE(ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.out.logs}/${params.data.di_mslist}.ps", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, ps_logs_dir)

    emit:
    RunPSPIPE.out.ready
}


workflow SplitMSChannels {
    take:
    ready

    main:

    stage_ch = channel.of('SB')

    stage_params_ch = InitParams(ready, stage_ch)

    Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    emit:
    Distribute.out
}


workflow ConcatTimeChunks {
    take:
    ready

    main:

    stage_ch = channel.of('CT')

    stage_params_ch = InitParams(ready, stage_ch)

    from_nodes_list = params.ssplit.nodes?.split(',') as List
    from_nodes_ch = channel.fromList(from_nodes_list).collect { it }

    to_nodes_list = params.data.nodes?.split(',') as List
    to_nodes_ch = channel.fromList(to_nodes_list).collect { it }

    per_sb_chunks_ch = GetTimeChunksPerSubband(stage_params_ch.params_standby, params.data.di_ms_glob.replace('_T???.MS', ''), params.data.bp_subbands_per_node, from_nodes_ch, to_nodes_ch)
    // "fullband_di_smooth_data" 

    Distribute(stage_params_ch.params_standby, per_sb_chunks_ch, stage_ch, stage_params_ch.params_file)

    emit:
    Distribute.out
}


workflow RunDIBandpass {
    take:
    ready

    main:

    stage_ch = channel.of('BP')

    stage_params_ch = InitParams(ready, stage_ch)

    cal_ch = Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    nodesList = params.data.nodes?.split(',') as List
    nodes_ch = channel.fromList(nodesList).collect { it }
    mses_ch = WriteMSlist(cal_ch, nodes_ch, "${params.data.path}/${params.data.bp_ms_glob}", params.data.bp_mslist)

    mses_sols_ch = ReadTxtLinesandAppend(mses_ch.per_line_mslist, params.out.logs, params.data.bp_mslist, "/${params.ddecal.bp.sols}")

    sols_collect_ch = H5ParmCollect(cal_ch, mses_sols_ch.list_postfix_str, "di_bandpass_solutions")

    aoq_comb_ch = AOqualityCombine(sols_collect_ch.done, mses_ch.per_line_mslist, "di_bandpass_aoq_stats")

    ps_dir = "/net/${params.master_node}/${params.data.path}/${params.out.results}/${params.pspipe.dir}/di_bandpass"

    rev_ch = AddRevision(aoq_comb_ch.done, params.pspipe.obsid, params.ddecal.bp.outcol, params.data.path, ps_dir, params.master_node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, 0, 0)

    RunPSPIPE(ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.out.logs}/${params.data.bp_mslist}", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, params.out.logs)

    emit:
    RunPSPIPE.out.ready
}


workflow Average {
    take:
    ready

    main:
    stage_ch = channel.of('AVG')

    stage_params_ch = InitParams(ready, stage_ch)

    avg_ch = Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    nodesList = params.data.nodes?.split(',') as List
    nodes_ch = channel.fromList(nodesList).collect { it }

    WriteMSlist(avg_ch, nodes_ch, "${params.data.path}/${params.data.bp_ms_glob}.dd_averaged", params.data.dd_sb_mslist)

    emit:

    WriteMSlist.out.single_line_mslist
}


workflow MakeFullBandDDTimeChunks {
    take:
    ready

    main:

    mses_ch = ReadTxtLinesandAppend(ready, params.out.logs, params.data.dd_sb_mslist, "/nan")

    nodesList = params.ssplit.nodes?.split(',') as List
    nodes_ch = channel.fromList(nodesList).collect { it }

    output_mslist = file(params.out.logs).resolve(params.data.dd_mslist)

    MergeChansSplitTime(ready, mses_ch.list_str, nodes_ch, params.ssplit.dd.ntimes, 'DATA', params.data.dd_ms_glob.replace('_T???.MS', '.MS'), params.ssplit.dd.mses_per_node, output_mslist)

    emit:

    MergeChansSplitTime.out
}


workflow Run_DD {
    take:
    ready

    main:

    stage_ch = channel.of('DD')

    stage_params_ch = InitParams(ready, stage_ch)

    cal_ch = Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    mses_sols_ch = ReadTxtLinesandAppend(cal_ch, params.out.logs, params.data.dd_mslist, "/${params.ddecal.dd.sols}")

    sols_collect_ch = H5ParmCollect(cal_ch, mses_sols_ch.list_postfix_str, "dd_smooth_solutions")

    aoq_comb_ch = AOqualityCombine(sols_collect_ch.done, "${params.out.logs}/${params.data.dd_mslist}", "dd_smooth_aoq_stats")

    ps_dir = "/net/${params.master_node}/${params.data.path}/${params.out.results}/${params.pspipe.dir}/dd_smooth"

    rev_ch = AddRevision(aoq_comb_ch.done, params.pspipe.obsid, params.ddecal.dd.outcol, params.data.path, ps_dir, params.master_node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, 0, 0)

    RunPSPIPE(ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.out.logs}/${params.data.dd_mslist}.ps", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, params.out.logs)

    emit:
    // AOqualityCombine.out.done
    // WScleanImage.out.done
    RunPSPIPE.out.ready
}


workflow Run_PostDD {
    take:
    ready

    main:
    stage_ch = channel.of('PDD')

    stage_params_ch = InitParams(ready, stage_ch)

    cal_ch = Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    nodesList = params.data.nodes?.split(',') as List
    nodes_ch = channel.fromList(nodesList).collect { it }
    mses_ch = WriteMSlist(cal_ch, nodes_ch, "${params.data.path}/${params.data.dd_ms_glob}.l${params.postdd.uvlambdamin}to${params.postdd.uvlambdamax}", params.data.pdd_mslist)

    AOqualityCombine(ready, mses_ch.per_line_mslist, "pdd_smooth_aoq_stats")

    emit:
    // WriteMSlist.out.single_line_mslist
    AOqualityCombine.out.done
}


workflow PowerSpectrum {
    take:
    ready

    main:
    // def mslist_file = file("${params.out.logs}/${params.data.pdd_mslist}")
    // if (mslist_file.exists()) {
    //     mses_ch = channel.fromPath( "${params.out.logs}/${params.data.pdd_mslist}", checkIfExists: true, type: 'file' )
    // } else {
    //     nodesList = params.data.nodes?.split(',') as List
    //     nodes_ch = channel.fromList( nodesList ).collect { it }
    //     mses_ch = WriteMSlist( ready, nodes_ch, "${params.data.path}/${params.data.dd_ms_glob}.l${params.postdd.uvlambdamin}to${params.postdd.uvlambdamax}", params.data.pdd_mslist).per_line_mslist
    // }

    // aoq_comb_ch = AOqualityCombine( ready, mses_ch, "pdd_smooth_aoq_stats" )

    ps_dir = "/net/${params.master_node}/${params.data.path}/${params.out.results}/${params.pspipe.dir}/post_dd"

    rev_ch = AddRevision(ready, params.pspipe.obsid, 'DATA', params.data.path, ps_dir, params.master_node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, 0, 0)

    RunPSPIPE(ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.out.logs}/${params.data.pdd_mslist}.ps", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, params.out.logs)

    emit:
    RunPSPIPE.out.ready
}


workflow Image {
    take:
    ready

    main:

    mses_sols_ch = ReadTxtLinesandAppend(ready, params.out.logs, params.data.pdd_mslist, "/nan")

    mses_and_imname_ch = mses_sols_ch.list_str.combine(channel.of("postdd_corrected"))

    // WScleanImage ( true, mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.postdd.beam.outcol )

    WScleanImage(true, mses_and_imname_ch, 2560, '2.0amin', params.wsclean.niter, 'iquv', 1, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, 'DATA')

    emit:

    WScleanImage.out.done
}


workflow EXTRACT {
    take:
    ready

    main:

    stage_ch = channel.of('TAR')

    stage_params_ch = InitParams(ready, stage_ch)

    Distribute(stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file)

    emit:
    Distribute.out
}

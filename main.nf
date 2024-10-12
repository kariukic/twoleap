#!/usr/bin/env nextflow

include {
    H5ParmCollect;
    AOqualityCombine;
    WScleanImage;
    ConcatFrequencySplitTime;
    ReadTxtLinesandAppend;
    WriteMSlist;
    writeHosts;
    parseNodes;
    makeDirectory;
    readTxtIntoString;
    readTxtAndAppendString;
} from "./processes.nf"


include {
    AddRevision;
    RunPSPIPE
} from "./makeps.nf"

/*
steps to run
DIcal
BPcal
// ApplyBeam
// PreDD
// AT
DD
PS
*/
workflow {
    // ext_ch = EXTRACT( true )
    fcab_ch = FCAB ( true )
    bp_ch = Run_BP( fcab_ch )
    split_ch = Split( bp_ch )
    di_ch = Run_DI ( split_ch )
    avg_ch = Average ( di_ch )
    dd_ch = Run_DD ( avg_ch )
    ps_ch = PowerSpectrum( dd_ch )
}
// Run_WS( dd_ch )

workflow TwoStep {
    bp_ch = Run_BP( true )
    split_ch = Split( bp_ch )
    di_ch = Run_DI ( split_ch )
    avg_ch = Average ( di_ch )
    at_ch = Run_AT( avg_ch )
    dd_ch = Run_DD ( at_ch )
    Run_WS( dd_ch )
}


workflow ExtractDATA {
    ext_ch = EXTRACT( true )
}

import groovy.json.JsonOutput
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
        makeDirectory( params.out.logs )

        def tasks_after_split = [ "DI", "AVG", "DD", "AT" ]

         if ( tasks_after_split.contains( stage ) ){
            nodes_list  = parseNodes( params.split.nodes )
         }

        else {
            nodes_list  = parseNodes( params.data.nodes )  
        }

        writeHosts ( nodes_list, params.data.hosts )

        """
        echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > ${stage}_params.json
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
        pssh -v -i -h ${launchDir}/${params.data.hosts} -t 0 -x "cd ${params.data.path}; bash" nextflow run ${params.stagelib} --stage ${entry} --ch_in ${ch_in} -params-file ${params_file} > ${params.out.logs}/${entry}.log 2>&1
        """
}


workflow FCAB {

    take:
        ready

    main:
        stage_ch = channel.of ( 'FCAB' )

        stage_params_ch = InitParams( ready, stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        // nodesList = params.data.nodes?.split(',') as List
        // nodes_ch = channel.fromList( nodesList ).collect {it}    
        // mslist_ch = WriteMSlist( cal_ch, nodes_ch, params.data.ms_files.raw, params.data.raw_mslist )
        // output_mslist = file(params.out.logs).resolve( params.data.raw_mslist )

        mses = readTxtIntoString ( params.data.raw_mslist )

        AOqualityCombine( cal_ch, mses, "aoqstats_raw" )

    emit:

        AOqualityCombine.out.qstats

}



workflow Run_BP {

    take:

        ready

    main:

        stage_ch = channel.of ( 'BP' )

        stage_params_ch = InitParams( ready, stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        mses = readTxtIntoString ( params.data.bp_mslist )

        // solution_files = readTxtAndAppendString( params.data.bp_mslist, "/${params.ddecal.bp.sols}" )
        // sols_collect_ch = H5ParmCollect( cal_ch, solution_files, "bp_combined_solutions" )

        aoq_comb_ch = AOqualityCombine( cal_ch, mses, "aoqstats_bp" )  // sols_collect_ch.combined_sols, mses

        mses_and_imname_ch = channel.of ( mses ).combine( channel.of ( "bp_corrected" ) )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.bp.outcol )

    emit:

        WScleanImage.out.done
        // AOqualityCombine.out.qstats
}


workflow Split {
    take:

        ready
    
    main:
        mses = readTxtIntoString ( params.data.bp_mslist )

        nodesList = params.split.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        output_mslist = file(params.out.logs).resolve( params.data.di_mslist )

        ConcatFrequencySplitTime ( ready, mses, nodes_ch, params.split.ntimes, params.ddecal.bp.outcol, params.split.ms_prefix, params.split.mses_per_node, output_mslist )
        
    emit:

       ConcatFrequencySplitTime.out

}


workflow Run_DI {

    take:

        ready

    main:

        stage_ch = channel.of ( 'DI' )

        stage_params_ch = InitParams( ready,  stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        mses_sols_ch = ReadTxtLinesandAppend( cal_ch, params.out.logs, params.data.di_mslist, "/${params.ddecal.di.sols}" )

        // sols_collect_ch = H5ParmCollect( true, mses_sols_ch.list_postfix_str, "di_combined_solutions" )

        aoq_comb_ch = AOqualityCombine( cal_ch, mses_sols_ch.list_str, "aoqstats_di" ) // sols_collect_ch.combined_sols

        mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "di_beam_corrected" ) )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.di.beam.outcol )

    emit:

        WScleanImage.out.done
        // AOqualityCombine.out.qstats

}


workflow Average {

    take:
        ready

    main:
        stage_ch = channel.of ( 'AVG' )

        stage_params_ch = InitParams( ready,  stage_ch )

        Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

    emit:

        Distribute.out
}


workflow Run_AT {

    take:
        ready

    main:

        stage_ch = channel.of ( 'AT' )

        stage_params_ch = InitParams( ready,  stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        nodesList = params.split.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        mses_ch = WriteMSlist( cal_ch, nodes_ch, params.data.ms_files.dd, params.data.dd_mslist)
        // mses_ch = WriteDDMSlist( cal_ch, nodes_ch)

        mses_sols_ch = ReadTxtLinesandAppend( mses_ch.per_line_mslist, params.out.logs, params.data.dd_mslist, "/${params.ddecal.ateams.sols}" )

        // sols_collect_ch = H5ParmCollect( true, mses_sols_ch.list_postfix_str, "dd_combined_solutions" )

        AOqualityCombine( true, mses_sols_ch.list_str, "aoqstats_ateams" ) // sols_collect_ch.combined_sols at p1

        // mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "dd_corrected" ) )

        // WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.polfit, params.ddecal.dd.outcol )

    emit:

        // WScleanImage.out.done
        AOqualityCombine.out.qstats
}


workflow Run_DD {

    take:
        ready

    main:

        stage_ch = channel.of ( 'DD' )

        stage_params_ch = InitParams( ready,  stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        nodesList = params.split.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        mses_ch = WriteMSlist( cal_ch, nodes_ch, params.data.ms_files.dd, params.data.dd_mslist)

        mses_sols_ch = ReadTxtLinesandAppend( mses_ch.per_line_mslist, params.out.logs, params.data.dd_mslist, "/${params.ddecal.dd.sols}" )

        // sols_collect_ch = H5ParmCollect( true, mses_sols_ch.list_postfix_str, "dd_combined_solutions" )

        AOqualityCombine( true, mses_sols_ch.list_str, "aoqstats_dd" ) // sols_collect_ch.combined_sols at p1 //aoq_comb_ch = 

        // mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "dd_corrected" ) )

        // WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.dd.outcol )

    emit:

        // WScleanImage.out.done
        AOqualityCombine.out.qstats
}


workflow Run_WS {

    take:
        ready

    main:

        stage_ch = channel.of ( 'WS' )

        stage_params_ch = InitParams( ready,  stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

}

workflow EXTRACT {

    take:
        ready

    main:

        stage_ch = channel.of ( 'TAR' )

        stage_params_ch = InitParams( ready,  stage_ch )

        Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

    emit:
        Distribute.out

}


workflow PowerSpectrum {
    take:
        ready

    main:
        nodes_list = params.split.nodes.split(',').collect{"node${it}"} as List

        nodesList = params.split.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        mses_ch = WriteMSlist( ready, nodes_ch, params.data.ms_files.dd, params.data.dd_mslist)
        // single_line_mslist= channel.of ("/home/codex/chege/projects/NCP2024/process/redshift2/L246309/logs/dd_mses.txt.ps_without_T017_and_T051")
        // single_line_mslist = channel.of ("/home/codex/chege/projects/NCP2024/process/redshift2/L246309/logs/di_mses.txt.ps")

        ps_dir = "/net/${nodes_list[0]}/${params.data.path}/${params.out.results}/${params.pspipe.dir}"
        rev_ch = AddRevision(mses_ch.single_line_mslist, params.pspipe.obsid, params.ddecal.dd.outcol, params.data.path, ps_dir, nodes_list[0], params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, 0, 0) // TODO: stop using only the final node
        ps_ch = RunPSPIPE(ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.out.logs}/${params.data.dd_mslist}.ps", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.gpr, params.pspipe.ml_gpr, params.out.logs)
}

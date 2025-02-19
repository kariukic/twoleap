#!/usr/bin/env nextflow

include {
    H5ParmCollect;
    AOqualityCombine;
    WScleanImage;
    MergeChansSplitTime;
    ReadTxtLinesandAppend;
    WriteMSlist;
    writeHosts;
    parseNodes;
    makeDirectory;
    readTxtIntoString;
    readTxtAndAppendString;
    GetTimeChunksPerSubband;
} from "./processes.nf"


include {
    AddRevision;
    RunPSPIPE
} from "./makeps.nf"


// workflow {
//     // ext_ch = EXTRACT( true )
//     fcab_ch = FCAB ( true )
//     bp_ch = RunDIBandpass( fcab_ch )
//     split_ch = MakeFullBandTimeChunks( bp_ch )
//     di_ch = RunDISmooth ( split_ch )
//     avg_ch = Average ( di_ch )
//     dd_ch = Run_DD ( avg_ch )
//     ps_ch = PowerSpectrum( dd_ch )
// }


// workflow TwoStep {
//     bp_ch = RunDIBandpass( true )
//     split_ch = MakeFullBandTimeChunks( bp_ch )
//     di_ch = RunDISmooth ( split_ch )
//     avg_ch = Average ( di_ch )
//     at_ch = Run_AT( avg_ch )
//     dd_ch = Run_DD ( at_ch )
//     Run_WS( dd_ch )
// }


workflow {
    // ext_ch = EXTRACT( true )
    fcab_ch = FCAB ( true ) //ext_ch )

    split_ch1 = MakeFullBandDITimeChunks( fcab_ch )
    di_ch = RunDISmooth ( split_ch1 )

    sbs_ch = SplitMSChannels( di_ch )
    sbs_ch2 = ConcatTimeChunks( sbs_ch )
    bp_ch = RunDIBandpass( sbs_ch2 )

    avg_ch = Average ( bp_ch )
    split_ch2 = MakeFullBandDDTimeChunks( true ) //avg_ch )
    dd_ch = Run_DD (  true ) //split_ch2 )
    ps_ch = PowerSpectrum( dd_ch )
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

        def tasks_after_split = [ "DI", "DD", "AT", "SB" ] //, "AVG"

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
        // mses_ch = WriteMSlist( cal_ch, nodes_ch, "${params.data.ms_files.raw}".replace( "_001", "_002" ), "di_subbands.txt" )

        mses = readTxtIntoString ( params.data.raw_mslist )

        AOqualityCombine( cal_ch, mses, "aoqstats_raw" )

    emit:

        AOqualityCombine.out.qstats

}


workflow MakeFullBandDITimeChunks {
    take:

        ready
    
    main:
        mses = readTxtIntoString ( params.data.di_mslist ) //TODO: check this mslist

        nodesList = params.split.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        output_mslist = file(params.out.logs).resolve( "di_mses.txt" )

        MergeChansSplitTime ( ready, mses, nodes_ch, params.split.di.ntimes, params.ddecal.di.incol, params.split.di.ms_prefix, params.split.mses_per_node, output_mslist )
        
    emit:

       MergeChansSplitTime.out

}


workflow RunDISmooth {

    take:

        ready

    main:

        stage_ch = channel.of ( 'DI' )

        stage_params_ch = InitParams( ready,  stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        mses_sols_ch = ReadTxtLinesandAppend( cal_ch, params.out.logs, "di_mses.txt", "/${params.ddecal.di.sols}" )

        AOqualityCombine( cal_ch, mses_sols_ch.list_str, "aoqstats_di_smooth" ) //aoq_comb_ch

        // mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "di_smooth_corrected" ) )

        // WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.di.outcol )

    emit:

        // WScleanImage.out.done
        AOqualityCombine.out.qstats

}


workflow SplitMSChannels {

    take:
        ready

    main:

        stage_ch = channel.of ( 'SB' )

        stage_params_ch = InitParams( ready,  stage_ch )

        Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

    emit:
        Distribute.out

}


workflow ConcatTimeChunks {

    take:
        ready

    main:

        stage_ch = channel.of ( 'CT' )

        stage_params_ch = InitParams( ready,  stage_ch )

        from_nodes_list = params.split.nodes?.split(',') as List
        from_nodes_ch = channel.fromList( from_nodes_list ).collect {it}

        to_nodes_list = params.data.nodes?.split(',') as List
        to_nodes_ch = channel.fromList( to_nodes_list ).collect {it}

        per_sb_chunks_ch = GetTimeChunksPerSubband( stage_params_ch.params_standby, params.split.di.ms_prefix, params.data.bp_subbands_per_node, from_nodes_ch, to_nodes_ch)

        Distribute ( stage_params_ch.params_standby, per_sb_chunks_ch, stage_ch, stage_params_ch.params_file )

    emit:
        Distribute.out

}



workflow RunDIBandpass {

    take:

        ready

    main:

        stage_ch = channel.of ( 'BP' )

        stage_params_ch = InitParams( ready, stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        nodesList = params.data.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}
        mses_ch = WriteMSlist( cal_ch, nodes_ch, params.data.ms_files.bp, "di_bandpass_subbands_mslist.txt")

        mses_sols_ch = ReadTxtLinesandAppend( mses_ch.per_line_mslist, params.out.logs, "di_bandpass_subbands_mslist.txt", "/${params.ddecal.bp.sols}" )

        aoq_comb_ch = AOqualityCombine( true, mses_sols_ch.list_str, "aoqstats_di_bandpass" ) 

        mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of ( "di_bandpass_beam_corrected" ) )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.beam.outcol )

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

        avg_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        nodesList = params.data.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        WriteMSlist( avg_ch, nodes_ch, "${params.data.ms_files.bp}".replace( "_002", "_003" ), "dd_subband_mslist.txt" )

    emit:

        WriteMSlist.out.single_line_mslist
}


workflow MakeFullBandDDTimeChunks {
    take:

        ready
    
    main:

        mses_ch = ReadTxtLinesandAppend( ready, params.out.logs, "dd_subband_mslist.txt", "/nan" )

        nodesList = params.split.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        output_mslist = file(params.out.logs).resolve( "dd_fullband_time_chunks_mslist.txt" )

        MergeChansSplitTime ( true, mses_ch.list_str, nodes_ch, params.split.dd.ntimes, 'DATA', params.split.dd.ms_prefix, params.split.mses_per_node, output_mslist )
        
    emit:

       MergeChansSplitTime.out

}


workflow Run_DD {

    take:
        ready

    main:

        stage_ch = channel.of ( 'DD' )

        stage_params_ch = InitParams( ready,  stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        mses_sols_ch = ReadTxtLinesandAppend( cal_ch, params.out.logs, "dd_fullband_time_chunks_mslist.txt", "/${params.ddecal.dd.sols}" )

        aoq_comb_ch = AOqualityCombine( true, mses_sols_ch.list_str, "aoqstats_dd" )

        mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "dd_corrected" ) )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.dd.outcol )

    emit:

        WScleanImage.out.done
        // AOqualityCombine.out.qstats
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


// Add this step after DI, BP and DD
workflow PowerSpectrum {
    take:
        ready

    main:
        nodes_list = params.split.nodes.split(',').collect{"node${it}"} as List

        // nodesList = params.split.nodes?.split(',') as List
        // nodes_ch = channel.fromList( nodesList ).collect {it}

        // mses_ch = WriteMSlist( ready, nodes_ch, params.data.ms_files.dd, "dd_fullband_time_chunks_mslist.txt")
        // single_line_mslist= channel.of ("/home/codex/chege/projects/NCP2024/process/redshift2/L246309/logs_sage_order_trial3/dd_mses_without_T17495051.txt.ps")
        // single_line_mslist = channel.of ("/home/codex/chege/projects/NCP2024/process/redshift2/L246309/logs/di_mses.txt.ps")

        ps_dir = "/net/${nodes_list[0]}/${params.data.path}/${params.out.results}/${params.pspipe.dir}"
        rev_ch = AddRevision(ready, params.pspipe.obsid, 'DATA', params.data.path, ps_dir, nodes_list[0], params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, 0, 0) // TODO: stop using only the final node @1 mses_ch.single_line_mslist //params.ddecal.dd.outcol

        // listed_bp_mses_ch = channel.of ("/home/codex/chege/projects/NCP2024/process/redshift2/L246309/bp_mslist.txt3") #/home/codex/chege/projects/modelling/L246309/di_cal_logs/di_bandpass_subbands_mslist.txt.ps_r2
        ps_ch = RunPSPIPE(ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.out.logs}/di_mses.txt", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.gpr, params.pspipe.ml_gpr, params.out.logs) //dd_fullband_time_chunks_mslist.txt.ps"

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

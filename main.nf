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


workflow FullPipeline {
    ext_ch = EXTRACT( true )
    pre_process_ch = PreProcess ( ext_ch )

    split_ch1 = MakeFullBandDITimeChunks( pre_process_ch )
    di_ch = RunDISmooth ( split_ch1 )

    sbs_ch = SplitMSChannels( di_ch )   
    sbs_ch2 = ConcatTimeChunks( sbs_ch )
    bp_ch = RunDIBandpass( sbs_ch2 )

    avg_ch = Average ( bp_ch )
    split_ch2 = MakeFullBandDDTimeChunks( avg_ch )
    dd_ch = Run_DD ( split_ch2 )
    pdd_ch = Run_PostDD ( dd_ch )
    ps_ch = PowerSpectrum( pdd_ch )
    fi_ch = FinalImages( ps_ch )
}

workflow {

    validStartPoints = ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd', 'postdd', 'powerspec', 'finalimages']

    // Validate the starting point parameter
    if (!(params.startFrom in validStartPoints)) {
        exit 1, "Invalid start point '${params.startFrom}'. Valid options: ${validStartPoints.join(', ')}"
    }

    // Define the workflow segments with conditional execution
    if (params.startFrom == 'extract') {
        ext_ch = EXTRACT(true)
    }
    
    if (params.startFrom == 'extract' || params.startFrom == 'pre_process') {
        pre_process_ch = params.startFrom == 'extract' ? PreProcess(ext_ch) : PreProcess(true)
    }
    
    if (params.startFrom == 'extract' || params.startFrom == 'pre_process' || params.startFrom == 'dismooth') {
        split_ch1 = params.startFrom in ['extract', 'pre_process'] ? 
                   MakeFullBandDITimeChunks(pre_process_ch) : 
                   MakeFullBandDITimeChunks(true)
        di_ch = RunDISmooth(split_ch1)
    }
    
    if (params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms']) {
        sbs_ch = params.startFrom in ['extract', 'pre_process', 'dismooth'] ? 
                SplitMSChannels(di_ch) : 
                SplitMSChannels(true)
        sbs_ch2 = ConcatTimeChunks(sbs_ch)
    }
    
    if (params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass']) {
        bp_ch = params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms'] ? 
               RunDIBandpass(sbs_ch2) : 
               RunDIBandpass(true)
    }
    
    if (params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average']) {
        avg_ch = params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass'] ? 
                Average(bp_ch) : 
                Average(true)
    }
    
    if (params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd']) {
        split_ch2 = params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average'] ? 
                   MakeFullBandDDTimeChunks(avg_ch) : 
                   MakeFullBandDDTimeChunks(true)
        dd_ch = Run_DD(split_ch2)
    }
    
    if (params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd', 'postdd']) {
        pdd_ch = params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd'] ? 
                Run_PostDD(dd_ch) : 
                Run_PostDD(true)
    }
    
    if (params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd', 'postdd', 'powerspec']) {
        ps_ch = params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd', 'postdd'] ? 
               PowerSpectrum(pdd_ch) : 
               PowerSpectrum(true)
    }
    
    // if (params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd', 'postdd', 'powerspec', 'finalimages']) {
    //     fi_ch = params.startFrom in ['extract', 'pre_process', 'dismooth', 'splitms', 'bandpass', 'average', 'dd', 'postdd', 'powerspec'] ? 
    //            FinalImages(ps_ch) : 
    //            FinalImages(true)
    // }
}


// workflow {

//     validSteps = [
//         'extract', 
//         'pre_process', 
//         'dismooth', 
//         'splitms', 
//         'bandpass',
//         'average', 
//         'dd', 
//         'postdd', 
//         'powerspec', 
//         'finalimages'
//     ]

//     // Validate parameters
//     if (!(params.startFrom in validSteps)) {
//         exit 1, "Invalid start point '${params.startFrom}'. Valid options: ${validSteps.join(', ')}"
//     }
//     if (!(params.stopAfter in validSteps)) {
//         exit 1, "Invalid stop point '${params.stopAfter}'. Valid options: ${validSteps.join(', ')}"
//     }
//     if (validSteps.indexOf(params.startFrom) > validSteps.indexOf(params.stopAfter)) {
//         exit 1, "Start point '${params.startFrom}' cannot be after stop point '${params.stopAfter}'"
//     }

//     // Helper function to check if step should run
//     def shouldRun = { step ->
//         def stepIdx = validSteps.indexOf(step)
//         def startIdx = validSteps.indexOf(params.startFrom)
//         def stopIdx = validSteps.indexOf(params.stopAfter)
//         return (stepIdx >= startIdx && stepIdx <= stopIdx)
//     }

//     // Helper function to get input for a step
//     def getInput = { step, normalInput ->
//         return shouldRun(step) && validSteps.indexOf(step) > validSteps.indexOf(params.startFrom) ? 
//                normalInput : 
//                true
//     }

//     // Workflow steps
//     if (shouldRun('extract')) {
//         ext_ch = EXTRACT(true)
//     }
    
//     if (shouldRun('pre_process')) {
//         pre_process_ch = getInput('pre_process', ext_ch) | PreProcess
//     }
    
//     if (shouldRun('dismooth')) {
//         split_ch1 = getInput('dismooth', pre_process_ch) | MakeFullBandDITimeChunks
//         di_ch = RunDISmooth(split_ch1)
//     }
    
//     if (shouldRun('splitms')) {
//         sbs_ch = getInput('splitms', di_ch) | SplitMSChannels
//         sbs_ch2 = ConcatTimeChunks(sbs_ch)
//     }
    
//     if (shouldRun('bandpass')) {
//         bp_ch = getInput('bandpass', sbs_ch2) | RunDIBandpass
//     }
    
//     if (shouldRun('average')) {
//         avg_ch = getInput('average', bp_ch) | Average
//     }
    
//     if (shouldRun('dd')) {
//         split_ch2 = getInput('dd', avg_ch) | MakeFullBandDDTimeChunks
//         dd_ch = Run_DD(split_ch2)
//     }
    
//     if (shouldRun('postdd')) {
//         pdd_ch = getInput('postdd', dd_ch) | Run_PostDD
//     }
    
//     if (shouldRun('powerspec')) {
//         ps_ch = getInput('powerspec', pdd_ch) | PowerSpectrum
//     }
    
//     if (shouldRun('finalimages')) {
//         fi_ch = getInput('finalimages', ps_ch) | FinalImages
//     }
// }



workflow ExtractDATA {
    EXTRACT( true )
}

// import groovy.json.JsonOutput
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

        def tasks_after_split = [ "DI", "DD", "PDD", "AT", "SB" ] //, "AVG"

         if ( tasks_after_split.contains( stage ) ){
            nodes_list  = parseNodes( params.split.nodes )
         }

        else {
            nodes_list  = parseNodes( params.data.nodes )  
        }

        writeHosts ( nodes_list, params.data.hosts )

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
        pssh -v -i -h ${launchDir}/${params.data.hosts} -t 0 -x "cd ${params.data.path}; bash" /home/codex/chege/software/nextflow run ${params.stagelib} --stage ${entry} --ch_in ${ch_in} -params-file ${params_file}  -with-singularity ${params.image} --singularity_bind_path ${params.binds} > ${params.out.logs}/${entry}.log 2>&1
        """
}


workflow PreProcess {

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

        MergeChansSplitTime ( ready, mses, nodes_ch, params.split.di.ntimes, params.ddecal.di.incol, params.split.di.ms_prefix, params.split.di.mses_per_node, output_mslist )
        
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

        mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of ( "di_bandpass_corrected" ) )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.bp.outcol )

    emit:

        WScleanImage.out.done
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

        MergeChansSplitTime ( true, mses_ch.list_str, nodes_ch, params.split.dd.ntimes, 'DATA', params.split.dd.ms_prefix, params.split.dd.mses_per_node, output_mslist )
        
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

        AOqualityCombine( true, mses_sols_ch.list_str, "aoqstats_dd" )
    emit:

        AOqualityCombine.out.qstats
}


workflow Run_PostDD {
    take:
        ready

    main:
        stage_ch = channel.of ( 'PDD' )

        stage_params_ch = InitParams( ready,  stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        mses_sols_ch = ReadTxtLinesandAppend( cal_ch, params.out.logs, "dd_fullband_time_chunks_mslist.txt", "/nan" )

        AOqualityCombine( true, mses_sols_ch.list_str, "aoqstats_postdd" )


    emit:

        AOqualityCombine.out.qstats
}


// Add this step after DI, BP and DD
workflow PowerSpectrum {
    take:
        ready

    main:

        nodesList = params.split.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}
        mses_ch = WriteMSlist( ready, nodes_ch, params.data.ms_files.dd, "dd_fullband_time_chunks_mslist.txt" )

        ps_dir = "/net/${params.master_node}/${params.data.path}/${params.out.results}/${params.pspipe.dir}"

        rev_ch = AddRevision( mses_ch.single_line_mslist, params.data.obsID, params.postdd.beam.outcol, params.data.path, ps_dir, params.master_node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, 0, 0 )

        RunPSPIPE( ps_dir, rev_ch.toml_file, params.data.obsID, "${params.out.logs}/dd_fullband_time_chunks_mslist.txt.ps_high_el", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, params.out.logs )

        // RunPSPIPE( ps_dir, rev_ch.toml_file, params.data.obsID, "${params.out.logs}/dd_fullband_time_chunks_mslist.txt.ps", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, params.out.logs )

    emit:
       RunPSPIPE.out.ready

}


workflow FinalImages {

    take:
        ready

    main:

        mses_sols_ch = ReadTxtLinesandAppend( ready, params.out.logs, "dd_fullband_time_chunks_mslist.txt", "/nan" )

        mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "postdd_corrected" ) )

        WScleanImage ( true, mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.postdd.beam.outcol )

    emit:

        WScleanImage.out.done

}


// workflow Run_WS {

//     take:
//         ready

//     main:

//         stage_ch = channel.of ( 'WS' )

//         stage_params_ch = InitParams( ready,  stage_ch )

//         Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

// }

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
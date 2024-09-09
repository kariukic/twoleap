#!/usr/bin/env nextflow

include {
    readTxtIntoString;
    readTxtAndAppendString;
    H5ParmCollect;
    AOqualityCombine;
    WScleanImage;
    ConcatFrequencySplitTime;
    ReadTxtLinesandAppend;
    AverageDItoDDMS;
    WriteDDMSlist;
} from "./processes.nf"

/*
steps to run
DIcal
BPcal
// ApplyBeam
// PreDD
DD
PS
*/
workflow {
    // // bp_ch = Run_BP( true )
    // // split_ch = Split( true ) //bp_ch )
    // di_ch = Run_DI ( true) //split_ch )
    // avg_ch = Average ( true ) //di_ch )
    dd_ch = Run_DD ( true ) //avg_ch )
    Run_WS( dd_ch )
}



import groovy.json.JsonOutput
process GetParams {
    debug true
    publishDir params.out.logs, mode: 'copy'

    input:
        val stage

    output:
        path "${stage}_params.json", emit: params_file
        val true, emit: params_standby

    script:
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
        pssh -v -i -h ${params.data.hosts} -t 0 -x "cd ${params.data.path}; bash" nextflow run ${params.stagelib} --stage ${entry} --ch_in ${ch_in} -params-file ${params_file} > ${params.out.logs}/${entry}.log 2>&1
        """
}


workflow Run_BP {

    take:

        stage_ready

    main:

        stage_ch = channel.of ( 'BP' )

        stage_params_ch = GetParams( stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, stage_ready, stage_ch, stage_params_ch.params_file )

        mses = readTxtIntoString ( params.data.mslist )

        solution_files = readTxtAndAppendString( params.data.mslist, "/${params.ddecal.bp.sols}" )

        // sols_collect_ch = H5ParmCollect( cal_ch, solution_files, "bp_combined_solutions" )

        aoq_comb_ch = AOqualityCombine( cal_ch, mses, "aoqstats_bp" )  // sols_collect_ch.combined_sols, mses

        mses_and_imname_ch = channel.of ( mses ).combine( channel.of ( "bp_corrected" ) )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.pol_fit, params.ddecal.bp.outcol )

    emit:

        WScleanImage.out.done
        // AOqualityCombine.out.qstats
}


workflow Split {
    take:

        ready
    
    main:
        mses = readTxtIntoString ( params.data.mslist )

        nodesList = params.data.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        output_mslist = file(params.out.logs).resolve( "di_mses.txt" )

        ConcatFrequencySplitTime ( ready, mses, nodes_ch, params.split.ntimes, params.ddecal.bp.outcol, params.split.msout, output_mslist )
        
    emit:

       ConcatFrequencySplitTime.out

}


workflow Run_DI {

    take:

        ready

    main:

        stage_ch = channel.of ( 'DI' )

        stage_params_ch = GetParams( stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        mses_sols_ch = ReadTxtLinesandAppend( cal_ch, params.out.logs, "di_mses.txt", "/${params.ddecal.di.sols}" )

        // sols_collect_ch = H5ParmCollect( true, mses_sols_ch.list_postfix_str, "di_combined_solutions" )

        aoq_comb_ch = AOqualityCombine( cal_ch, mses_sols_ch.list_str, "aoqstats_di" ) // sols_collect_ch.combined_sols

        mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "di_beam_corrected" ) )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.pol_fit, params.ddecal.di.beam.outcol )

    emit:

        WScleanImage.out.done
        // AOqualityCombine.out.qstats

}


workflow Average {

    take:
        ready

    main:
        stage_ch = channel.of ( 'AVG' )

        stage_params_ch = GetParams( stage_ch )

        Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

    emit:

        Distribute.out
}


workflow Run_DD {

    take:
        ready

    main:

        stage_ch = channel.of ( 'DD' )

        stage_params_ch = GetParams( stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

        nodesList = params.data.nodes?.split(',') as List
        nodes_ch = channel.fromList( nodesList ).collect {it}

        mses_ch = WriteDDMSlist( cal_ch, nodes_ch)

        mses_sols_ch = ReadTxtLinesandAppend( mses_ch, params.out.logs, "dd_mses.txt", "/${params.ddecal.dd.sols}" )

        // sols_collect_ch = H5ParmCollect( true, mses_sols_ch.list_postfix_str, "dd_combined_solutions" )

        AOqualityCombine( true, mses_sols_ch.list_str, "aoqstats_dd" ) // sols_collect_ch.combined_sols at p1

        // mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "dd_corrected" ) )

        // WScleanImage ( aoq_comb_ch.qstats.collect(), mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.pol_fit, params.ddecal.dd.outcol )

    emit:

        // WScleanImage.out.done
        AOqualityCombine.out.qstats
}


workflow Run_WS {


    take:
        ready

    main:

        stage_ch = channel.of ( 'WS' )

        stage_params_ch = GetParams( stage_ch )

        cal_ch = Distribute ( stage_params_ch.params_standby, ready, stage_ch, stage_params_ch.params_file )

}
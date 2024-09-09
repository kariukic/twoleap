#!/usr/bin/env nextflow

include {
    ScaleData;
    ClipData;
    DP3Calibrate;
    ApplyGains;
    AOqualityCollect;
    ApplyBEAM;
    AverageDItoDDMS;
    MakeDP3ClustersListFile;
    SubtractSources;
    WScleanImage;
    readTxtIntoString;
} from './processes.nf'

workflow {

    if ( params.stage == "BP" ) {
        
        BP ( params.ch_in )

    }

    if ( params.stage == "DI" ) {
        
        DI ( params.ch_in )

    }

    if ( params.stage == "AVG" ) {
        
        AVG ( params.ch_in )

    }

    if ( params.stage == "DD" ) {
        
        DD ( params.ch_in )

    }

    if ( params.stage == "WS" ) {
        
        WS ( params.ch_in )

    }

}


workflow BP {

    take:
        start_ch

    main:

        // String mspattern= params.ms.split(',').collect{"${it}"}.join("*")
        mset_ch = channel.fromPath( params.data.ms, glob: true, checkIfExists: true, type: 'dir' )

        // scale_data_ch = ScaleData( mset_ch )

        clip_ch = ClipData( mset_ch ) //scale_data_ch 

        sols_ch = DP3Calibrate( true, clip_ch, params.ddecal.bp.parset, params.ddecal.bp.sourcedb, params.ddecal.bp.sols, params.ddecal.bp.incol, params.ddecal.bp.solint ) // at postion 2

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.bp.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.bp.apply.parset, params.ddecal.bp.incol, params.ddecal.bp.outcol ) //  at position 1

        AOqualityCollect( true, apply_gains_ch, params.ddecal.bp.outcol )
    
    emit:

        AOqualityCollect.out

}


workflow DI {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath( "${params.data.path}/${params.split.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        sols_ch = DP3Calibrate( start_ch, mset_ch, params.ddecal.di.parset, params.ddecal.di.sourcedb, params.ddecal.di.sols, params.ddecal.di.incol, params.ddecal.di.solint )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.di.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.di.apply.parset, params.ddecal.di.incol, params.ddecal.di.outcol) //

        apply_elbeam_ch = ApplyBEAM(true, apply_gains_ch, params.ddecal.di.beam.parset, params.ddecal.di.outcol, params.ddecal.di.beam.outcol )

        AOqualityCollect( true, apply_elbeam_ch, params.ddecal.di.beam.outcol )
    
    emit:

        AOqualityCollect.out

}


workflow AVG {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath( "${params.data.path}/${params.split.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        AverageDItoDDMS ( start_ch, mset_ch, params.ddecal.di.beam.outcol, params.average.ditodd.timestep, params.average.ditodd.freqstep )

    emit:

        AverageDItoDDMS.out.done_averaging

}


workflow DD {
    take:
        start_ch

    main:

        mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        // cluster_ch = MakeClusters( source_select_ch, params.number_of_clusters, "l3_ncp_clusters.ao" )
        // OR
        // add a clustrsfile in the nextflow.config file

        clip_ch = ClipData( mset_ch )

        calibrate_ch = DP3Calibrate( start_ch, clip_ch, params.ddecal.dd.parset, params.ddecal.dd.sourcedb, params.ddecal.dd.sols, params.ddecal.dd.incol, params.ddecal.dd.solint )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.dd.sols}" }

        // clusters_ch  = MakeDP3ClustersListFile( calibrate_ch.collect(), params.number_of_clusters, "clusters_list.txt" )
        clusters_ch = channel.of( params.ddecal.dd.subtract.clusters )

        mset_and_sourcedb_ch = mset_ch.flatten().combine( channel.of( params.ddecal.dd.sourcedb ) )

        mset_sourcedb_solutions_and_clusters_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).combine( clusters_ch )

        subtract_ncp_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_and_clusters_ch, params.ddecal.dd.subtract.parset, params.ddecal.dd.incol, params.ddecal.dd.outcol )

        AOqualityCollect( true, subtract_ncp_ch, params.ddecal.dd.outcol )

    emit:

        AOqualityCollect.out

}


workflow WS {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        // mset_ch = channel.fromPath( "${params.data.path}/${params.split.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_dd_corrected_wide" }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( start_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.pol_fit, params.ddecal.dd.outcol )  // params.ddecal.di.beam.outcol ) //

}
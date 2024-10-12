#!/usr/bin/env nextflow

include {
    ScaleData;
    ClipData;
    DP3CalibrateDI;
    DP3CalibrateDD;
    ApplyGains;
    AOqualityCollect;
    ApplyBEAM;
    MakeDP3ClustersListFile;
    SubtractSources;
    WScleanImage;
    Flag;
    Compress;
    Average;
    UnpackMSTarball;
    readTxtIntoString;
} from './processes.nf'

workflow {

    if ( params.stage == "TAR" ) {
        
        TAR ( params.ch_in )

    }

    if ( params.stage == "FCAB" ) {
        
        FCAB ( params.ch_in )

    }

    if ( params.stage == "BP" ) {
        
        BP ( params.ch_in )

    }

    if ( params.stage == "DI" ) {
        
        DI ( params.ch_in )

    }

    if ( params.stage == "AVG" ) {
        
        AVG ( params.ch_in )

    }

    if ( params.stage == "AT" ) {
        
        ATEAMS ( params.ch_in )

    }

    if ( params.stage == "DD" ) {
        
        DD ( params.ch_in )

    }

    if ( params.stage == "WS" ) {
        
        WS ( params.ch_in )

    }

}


//FlagCompressBackupAverage
workflow FCAB {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath( params.data.ms_files.raw, glob: true, checkIfExists: true, type: 'dir' )

        compress_ch = Compress( mset_ch )

        compressed_mset_ch = mset_ch.collect { "${params.data.path}/" + it.getName().replace( ".MS", ".DCMS" ) }

        averaged_msnames_ch = mset_ch.collect { it.getName().replace( "_001", "_002" ) } //"15ch2s", "1ch4s" ) } //TODO: replace thes numbers with label params
        
        all_msets_and_averaged_msnames_ch = compressed_mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        avg_ch = Average ( compress_ch.collect(), all_msets_and_averaged_msnames_ch, params.average.lta_to_di.column, params.average.lta_to_di.timestep, params.average.lta_to_di.freqstep )

        AOqualityCollect( avg_ch.done_averaging.collect(), mset_ch, params.average.lta_to_di.column )
        
        // BackUP ( )

    emit:

        AOqualityCollect.out

}


workflow BP {

    take:
        start_ch

    main:

        // String mspattern= params.ms.split(',').collect{"${it}"}.join("*")
        mset_ch = channel.fromPath( params.data.ms_files.bp, glob: true, checkIfExists: true, type: 'dir' )

        scale_ch = ScaleData( mset_ch )

        clip_ch = ClipData( scale_ch )

        sols_ch = DP3CalibrateDI( true, clip_ch, params.ddecal.bp.parset, params.ddecal.bp.sourcedb, params.ddecal.bp.sols, params.ddecal.bp.incol, params.ddecal.bp.solint )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.bp.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.bp.apply.parset, params.ddecal.bp.incol, params.ddecal.bp.outcol )

        aoq_ch = AOqualityCollect( true, apply_gains_ch, params.ddecal.bp.outcol )
    
        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, 3, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.bp.outcol )

    emit:

        WScleanImage.out.done


}


workflow DI {

    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "${params.data.path}/${params.split.ms_prefix}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.di, glob: true, checkIfExists: true, type: 'dir' )

        sols_ch = DP3CalibrateDI( start_ch, mset_ch, params.ddecal.di.parset, params.ddecal.di.sourcedb, params.ddecal.di.sols, params.ddecal.di.incol, params.ddecal.di.solint )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.di.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.di.apply.parset, params.ddecal.di.incol, params.ddecal.di.outcol) //

        apply_elbeam_ch = ApplyBEAM(true, apply_gains_ch, params.ddecal.di.beam.parset, params.ddecal.di.outcol, params.ddecal.di.beam.outcol )

        aoq_ch = AOqualityCollect( true, apply_elbeam_ch, params.ddecal.di.beam.outcol )
    
        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, 3, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.di.beam.outcol )

    emit:

        WScleanImage.out.done

}


workflow AVG {

    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "${params.data.path}/${params.split.ms_prefix}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.di, glob: true, checkIfExists: true, type: 'dir' )

        averaged_msnames_ch = mset_ch.collect { it.getName().replace( "_002", "_003" ) }
        
        all_msets_and_averaged_msnames_ch = mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        Average ( true, all_msets_and_averaged_msnames_ch, params.ddecal.di.beam.outcol, params.average.ditodd.timestep, params.average.ditodd.freqstep )

    emit:

        Average.out.done_averaging

}

workflow DD {
    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.dd, glob: true, checkIfExists: true, type: 'dir' )

        // cluster_ch = MakeClusters( source_select_ch, params.number_of_clusters, "l3_ncp_clusters.ao" )
        // OR
        // add a clustrsfile in the nextflow.config file

        clip_ch = ClipData( mset_ch )

        calibrate_ch = DP3CalibrateDD( start_ch, clip_ch, params.ddecal.dd.parset, params.ddecal.dd.sourcedb, params.ddecal.dd.sols, params.ddecal.dd.incol, params.ddecal.dd.solint )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.dd.sols}" }

        // clusters_ch  = MakeDP3ClustersListFile( calibrate_ch.collect(), params.number_of_clusters, "clusters_list.txt" )
        clusters_ch = channel.of( params.ddecal.dd.subtract.clusters )

        mset_and_sourcedb_ch = mset_ch.flatten().combine( channel.of( params.ddecal.dd.sourcedb ) )

        mset_sourcedb_solutions_and_clusters_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).combine( clusters_ch )

        subtract_ncp_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_and_clusters_ch, params.ddecal.dd.subtract.parset, params.ddecal.dd.incol, params.ddecal.dd.outcol )

        aoq_ch = AOqualityCollect( true, subtract_ncp_ch, params.ddecal.dd.outcol )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, 3, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.dd.outcol )

    emit:

        WScleanImage.out.done

}


workflow ATEAMS {
    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.dd, glob: true, checkIfExists: true, type: 'dir' )

        // flag_ch = Flag (mset_ch, params.ddecal.ateams.outcol)

        clip_ch = ClipData( mset_ch )

        calibrate_ch = DP3CalibrateDD( start_ch, clip_ch, params.ddecal.ateams.parset, params.ddecal.ateams.sourcedb, params.ddecal.ateams.sols, params.ddecal.ateams.incol, params.ddecal.ateams.solint )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.ateams.sols}" }

        // clusters_ch  = MakeDP3ClustersListFile( calibrate_ch.collect(), params.number_of_clusters, "clusters_list.txt" )
        clusters_ch = channel.of( params.ddecal.ateams.subtract.clusters )

        mset_and_sourcedb_ch = mset_ch.flatten().combine( channel.of( params.ddecal.ateams.sourcedb ) )

        mset_sourcedb_solutions_and_clusters_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).combine( clusters_ch )

        subtract_ncp_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_and_clusters_ch, params.ddecal.ateams.subtract.parset, params.ddecal.ateams.incol, params.ddecal.ateams.outcol )

        aoq_ch = AOqualityCollect( true, subtract_ncp_ch, params.ddecal.ateams.outcol )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, 3, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.ateams.outcol )

    emit:

        WScleanImage.out.done

}


workflow WS {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*flagged.MS", glob: true, checkIfExists: true, type: 'dir' )
        // mset_ch = channel.fromPath( "${params.data.path}/${params.split.ms_prefix}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( start_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.wsclean.column )  //params.ddecal.di.beam.outcol ) // params.ddecal.dd.outcol

}


workflow TAR {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath( params.data.ms_files.raw, glob: true, checkIfExists: true, type: 'file' )

        untar_ch = UnpackMSTarball( mset_ch, params.data.label )

}
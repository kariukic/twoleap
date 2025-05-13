#!/usr/bin/env nextflow

include {
    ScaleData;
    ClipData;
    FlagStations;
    DP3CalibrateDI;
    DP3GainCalDI;
    DP3CalibrateDD;
    ApplyGains;
    AOqualityCollect;
    ApplyBEAM;
    MakeDP3ClustersListFile;
    SubtractSources;
    WScleanImage;
    AOFlag;
    FlagIntra;
    UVWFlag;
    Compress;
    FilterInter;
    Average;
    Demix;
    UnpackMSTarball;
    SplitMSToSubbands;
    ConcatMSesinTime;
    FitBpol;
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
        
        DIBandpass ( params.ch_in )

    }

    if ( params.stage == "DI" ) {
        
        DISmooth ( params.ch_in )

    }

    if ( params.stage == "SB" ) {
        
        SB ( params.ch_in )

    }

    if ( params.stage == "CT" ) {
        
        CT ( params.ch_in )

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

    if ( params.stage == "PDD" ) {

        PostDD ( params.ch_in  )

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
        flag_intra_ch = FlagIntra ( start_ch, mset_ch )

        filter_ch = FilterInter ( flag_intra_ch.collect(), mset_ch )
        
        filtered_mset_ch = filter_ch.collect { "${params.data.path}/" + it.getName().replace( '.FMS.5ch4s.dppp', '.FMS.5ch4s.dppp' ) } //15ch2s replace!!!!!!!!!!!!!

        // demix_ch = Demix(params.demix.enabled, filtered_mset_ch.flatten(), params.demix.parset, params.demix.sourcedb)
        // demixed_mset_ch = filter_ch.collect { "${params.data.path}/" + it.getName().replace( '.FDMS', '.FDMS' ) } //demix_ch.collect

        flag_ch = AOFlag ( true, filtered_mset_ch.flatten(), params.average.lta_to_di.column, params.average.lta_to_di.aoflagger_strategy, 1 ) //demixed_mset_ch.flatten() @ 2

        averaged_msnames_ch = filter_ch.collect { it.getName().replace( ".FMS.5ch4s.dppp", "_002_3c196.MS" ) } // "_001", "_002" ) } //TODO: replace thes numbers with label params
        all_msets_and_averaged_msnames_ch = filtered_mset_ch.flatten().merge( averaged_msnames_ch.flatten() ) //demixed_mset_ch.flatten().merge ...

        avg_ch = Average ( flag_ch.collect(), all_msets_and_averaged_msnames_ch, params.average.lta_to_di.column, params.average.lta_to_di.timestep, params.average.lta_to_di.freqstep  )

        flagged_filtered_averaged_mset_ch = mset_ch.collect { "${params.data.path}/" + it.getName().replace(  ".MS.5ch4s.dppp", "_002_3c196.MS"  ) }

        AOqualityCollect( avg_ch.done_averaging.collect(), flagged_filtered_averaged_mset_ch.flatten(), params.average.lta_to_di.column )

    emit:

        AOqualityCollect.out 

}


workflow DISmooth {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath( params.data.ms_files.di, glob: true, checkIfExists: true, type: 'dir' )

        sols_ch = DP3CalibrateDI( start_ch, mset_ch, params.ddecal.di.parset, params.ddecal.di.sourcedb, params.ddecal.di.sols, params.ddecal.di.incol, params.ddecal.di.solint, params.ddecal.di.uvlambdamin, params.ddecal.di.uvlambdamax, params.ddecal.di.nchan )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.di.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.di.apply.parset, params.ddecal.di.incol, params.ddecal.di.outcol)

        aoq_ch = AOqualityCollect( true, apply_gains_ch, params.ddecal.di.outcol )
    
        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.di.outcol )

    emit:

        WScleanImage.out.done

}


workflow SB {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath( params.data.ms_files.di, glob: true, checkIfExists: true, type: 'dir' )

        SplitMSToSubbands( start_ch, mset_ch, params.average.nchans_per_subband_after_averaging, params.ddecal.di.outcol )

    emit:

        SplitMSToSubbands.out

}

workflow CT {

    take:
        start_ch

    main:
        
        txts_ch = channel.fromPath( "${params.data.path}/${params.split.di.ms_prefix}_*.txt", glob: true, checkIfExists: true )

        msouts_ch =  txts_ch.collect { it.getSimpleName() + ".MS" }

        txts_and_msouts_ch = txts_ch.merge( msouts_ch.flatten() )

        ConcatMSesinTime( start_ch, txts_and_msouts_ch )

    emit:

        ConcatMSesinTime.out.done

}



workflow DIBandpass {

    take:
        start_ch

    main:

        // String mspattern= params.ms.split(',').collect{"${it}"}.join("*")
        // mset_ch = channel.fromPath( params.data.ms_files.bp, glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.bp, glob: true, checkIfExists: true, type: 'dir' )

        sols_ch = DP3GainCalDI( start_ch, mset_ch, params.ddecal.bp.parset, params.ddecal.bp.sourcedb, params.ddecal.bp.sols, params.ddecal.bp.incol, params.ddecal.bp.solint, params.ddecal.bp.uvlambdamin, params.ddecal.bp.uvlambdamax, params.ddecal.bp.nchan )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.bp.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.bp.apply.parset, params.ddecal.bp.incol, params.ddecal.bp.outcol )

        flag_ch = AOFlag ( true, apply_gains_ch, params.ddecal.bp.outcol, params.ddecal.bp.aoflagger_strategy, 0 )

        aoq_ch = AOqualityCollect( flag_ch, apply_gains_ch, params.ddecal.bp.outcol) //params.ddecal.beam.outcol ) // @ 1
    
        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_subband, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.bp.outcol ) //params.ddecal.beam.outcol )

    emit:

        WScleanImage.out.done

}


workflow AVG {

    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "${params.data.path}/${params.split.ms_prefix}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.bp, glob: true, checkIfExists: true, type: 'dir' )

        averaged_msnames_ch = mset_ch.collect { it.getName().replace( "_002", "_003" ) }
        
        all_msets_and_averaged_msnames_ch = mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        Average ( start_ch, all_msets_and_averaged_msnames_ch, params.ddecal.bp.outcol, params.average.ditodd.timestep, params.average.ditodd.freqstep ) //params.ddecal.beam.outcol

    emit:

        Average.out.done_averaging

}


workflow DD {
    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.dd, glob: true, checkIfExists: true, type: 'dir' )

        calibrate_ch = DP3CalibrateDD( start_ch, mset_ch, params.ddecal.dd.parset, params.ddecal.dd.sourcedb, params.ddecal.dd.sols, params.ddecal.dd.incol, params.ddecal.dd.solint, params.ddecal.dd.uvlambdamin, params.ddecal.dd.uvlambdamax, params.ddecal.dd.nchan, params.ddecal.dd.flagstations )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.dd.sols}" }

        mset_and_sourcedb_ch = mset_ch.flatten().combine( channel.of( params.ddecal.dd.sourcedb ) )

        mset_sourcedb_solutions_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() )

        subtract_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_ch, params.ddecal.dd.subtract.parset, params.ddecal.dd.incol, params.ddecal.dd.outcol )

        aoq_ch = AOqualityCollect( true, subtract_ch, params.ddecal.dd.outcol )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.dd.outcol )

    emit:

        WScleanImage.out.done

}


workflow PostDD {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath( params.data.ms_files.dd, glob: true, checkIfExists: true, type: 'dir' )

        // flag_intra_ch = FlagIntra ( start_ch, mset_ch )

        aoflagger_ch = AOFlag ( start_ch, mset_ch, params.ddecal.dd.outcol, params.postdd.aoflagger_strategy, 1 )

        uvwflag_ch = UVWFlag (aoflagger_ch.collect(), mset_ch, params.ddecal.dd.outcol, params.postdd.uvlambdamin, params.postdd.uvlambdamax)

        beam_ch = ApplyBEAM( uvwflag_ch.collect(), mset_ch, params.postdd.beam.parset, params.ddecal.dd.outcol, params.postdd.beam.outcol )

        aoq_ch = AOqualityCollect( true, beam_ch, params.postdd.beam.outcol )

        im_names_ch =  mset_ch.collect { it.getSimpleName()  + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.postdd.beam.outcol )

    emit:

        WScleanImage.out.done

}


workflow ATEAMS {
    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.ms_files.dd, glob: true, checkIfExists: true, type: 'dir' )

        // flag_ch = AOFlag (mset_ch, params.ddecal.ateams.outcol)

        clip_ch = ClipData( mset_ch )

        calibrate_ch = DP3CalibrateDD( start_ch, clip_ch, params.ddecal.ateams.parset, params.ddecal.ateams.sourcedb, params.ddecal.ateams.sols, params.ddecal.ateams.incol, params.ddecal.ateams.solint, params.ddecal.dd.uvlambdamin, params.ddecal.dd.uvlambdamax, params.ddecal.dd.nchan, params.ddecal.dd.flagstations )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.ateams.sols}" }

        // clusters_ch  = MakeDP3ClustersListFile( calibrate_ch.collect(), params.number_of_clusters, "clusters_list.txt" )
        clusters_ch = channel.of( params.ddecal.ateams.subtract.clusters )

        mset_and_sourcedb_ch = mset_ch.flatten().combine( channel.of( params.ddecal.ateams.sourcedb ) )

        mset_sourcedb_solutions_and_clusters_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).combine( clusters_ch )

        subtract_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_and_clusters_ch, params.ddecal.ateams.subtract.parset, params.ddecal.ateams.incol, params.ddecal.ateams.outcol )

        aoq_ch = AOqualityCollect( true, subtract_ch, params.ddecal.ateams.outcol )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.ateams.outcol )

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

        WScleanImage ( start_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.wsclean.column )  //params.ddecal.beam.outcol ) // params.ddecal.dd.outcol

}


workflow TAR {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath( params.data.ms_files.raw, glob: true, checkIfExists: true, type: 'file' )

        untar_ch = UnpackMSTarball( mset_ch, params.data.label )

}
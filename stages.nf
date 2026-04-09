#!/usr/bin/env nextflow

include {
    ScaleData;
    ClipData;
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
    FlagBaselines;
    readTxtIntoString;
    GetMSColumn;
} from './processes.nf'

    // FlagStations;

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

        mset_ch = channel.fromPath( params.data.raw_ms_glob, glob: true, checkIfExists: true, type: 'dir' )
        flag_intra_ch = FlagIntra ( start_ch, mset_ch )

        filter_ch = FilterInter ( flag_intra_ch.collect(), mset_ch )

        filtered_mset_ch = filter_ch.collect { mset -> "${params.data.path}/" + mset.getName() }
        
        flag_ch = AOFlag ( true, filtered_mset_ch.flatten(), params.average.lta_to_di.column, params.average.lta_to_di.aoflagger_strategy, 1, 1)

        averaged_msnames_ch = filter_ch.collect { mset -> mset.getName() + '.flagged.di_averaged' }

        all_msets_and_averaged_msnames_ch = filtered_mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        // mset_ch = channel.fromPath( params.data.raw_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

        // filtered_mset_ch = mset_ch.collect { mset -> "${params.data.path}/" + mset.getName()+ '.noInter' }

        // averaged_msnames_ch = mset_ch.collect { mset -> mset.getName() + '.noInter.flagged.di_averaged3' } //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // all_msets_and_averaged_msnames_ch = filtered_mset_ch.flatten().merge( averaged_msnames_ch.flatten() )


        //////////////////////////////////////////////////////////////////////////////////////////

        avg_ch = Average ( flag_ch.collect(), all_msets_and_averaged_msnames_ch, params.average.lta_to_di.column, params.average.lta_to_di.timestep, params.average.lta_to_di.freqstep  )

        flagged_filtered_averaged_mset_ch = mset_ch.collect { mset -> "${params.data.path}/" + mset.getName() + '.noInter.flagged.di_averaged'  }

        AOqualityCollect( avg_ch.done_averaging.collect(), flagged_filtered_averaged_mset_ch.flatten(), params.average.lta_to_di.column )
    emit:
        // Average.out.done_averaging
        AOqualityCollect.out 

}
//TODO:
// turn on aoflager
// return avg_ch
// turn on aoquality collect

workflow DISmooth {

    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "fullband_di_smooth_data_T???.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.di_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

        sols_ch = DP3CalibrateDI( start_ch, mset_ch, params.ddecal.di.parset, params.ddecal.di.sourcedb, params.ddecal.di.sols, params.ddecal.di.incol, params.ddecal.di.solint, params.ddecal.di.uvlambdamin, params.ddecal.di.uvlambdamax, params.ddecal.di.uvmmax, params.ddecal.di.nchan, params.ddecal.di.flagstations, params.ddecal.di.calmode, params.ddecal.di.smoothnessconstraint, params.ddecal.di.propagate_sols, params.ddecal.di.maxiter, params.ddecal.di.beamproximitylimit, params.ddecal.di.usebeam, params.ddecal.di.beammode, params.ddecal.di.propagate_converged_sols_only, params.ddecal.di.scaling_coefficient )

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

        // mset_ch = channel.fromPath( "fullband_di_smooth_data_T???.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.di_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

        SplitMSToSubbands( start_ch, mset_ch, params.average.nchans_per_subband_after_averaging, params.ddecal.di.outcol )

    emit:

        SplitMSToSubbands.out

}

workflow CT {

    take:
        start_ch

    main:
        def txts_glob=params.data.di_ms_glob.replace('_T???.MS', '_SB???.txt')

        def di_ms_stem = params.data.di_ms_glob.replace('_T???.MS', '')
        def bp_ms_stem = params.data.bp_ms_glob.replace('_SB???.MS', '')

        txts_ch = channel.fromPath( "${params.data.path}/${txts_glob}", glob: true, checkIfExists: true ) // fullband_di_smooth_data_SB???.txt"

        msouts_ch =  txts_ch.collect { it.getSimpleName().replace(di_ms_stem, bp_ms_stem ) + '.MS' } // //"fullband_di_smooth", "di_bandpass")

        txts_and_msouts_ch = txts_ch.merge( msouts_ch.flatten() )

        ConcatMSesinTime( start_ch, txts_and_msouts_ch )

    emit:

        ConcatMSesinTime.out.done

}



workflow DIBandpass {

    take:
        start_ch

    main:

        // String mspattern= params.ms.ssplit(',').collect{"${it}"}.join("*")
        mset_ch = channel.fromPath( params.data.bp_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

        sols_ch = DP3GainCalDI( start_ch, mset_ch, params.ddecal.bp.parset, params.ddecal.bp.sourcedb, params.ddecal.bp.sols, params.ddecal.bp.incol, params.ddecal.bp.solint, params.ddecal.bp.uvlambdamin, params.ddecal.bp.uvlambdamax, params.ddecal.bp.uvmmax, params.ddecal.bp.nchan )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.bp.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.bp.apply.parset, params.ddecal.bp.incol, params.ddecal.bp.outcol )

        flag_ch = AOFlag ( true, apply_gains_ch, params.ddecal.bp.outcol, params.ddecal.bp.aoflagger_strategy, 1, 0 )

        aoq_ch = AOqualityCollect( flag_ch, apply_gains_ch, params.ddecal.bp.outcol) //params.ddecal.beam.outcol )
    
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

        // mset_ch = channel.fromPath( "di_bandpass_data_SB???.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.bp_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

        averaged_msnames_ch = mset_ch.collect { it.getName() + '.dd_averaged' } //.replace( "_002", "_003" ) }
        
        all_msets_and_averaged_msnames_ch = mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        Average ( start_ch, all_msets_and_averaged_msnames_ch, params.ddecal.bp.outcol, params.average.ditodd.timestep, params.average.ditodd.freqstep ) //params.ddecal.beam.outcol

    emit:

        Average.out.done_averaging

}


workflow DD {
    take:
        start_ch

    main:

        // mset_ch = channel.fromPath( "fullband_dd_smooth_data_T???.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.dd_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

        calibrate_ch = DP3CalibrateDD( 
            start_ch, mset_ch, params.ddecal.dd.parset, params.ddecal.dd.sourcedb, 
            params.ddecal.dd.sols, params.ddecal.dd.incol, params.ddecal.dd.calmode, 
            params.ddecal.dd.solint, params.ddecal.dd.uvlambdamin, 
            params.ddecal.dd.uvlambdamax, params.ddecal.dd.uvmmax, params.ddecal.dd.nchan, 
            params.ddecal.dd.usebeam, params.ddecal.dd.beammode, 
            params.ddecal.dd.smoothnessconstraint, params.ddecal.dd.truncateksmoothkernel, 
            params.ddecal.dd.robust_reg, params.ddecal.dd.propagate_sols, 
            params.ddecal.dd.maxiter, params.ddecal.dd.beamproximitylimit, 
            params.ddecal.dd.correctfreqsmearing, params.ddecal.dd.flagstations, 
            params.ddecal.dd.propagate_converged_sols_only, 0.15, 1, 'directioniterative'
        )

        // Create tuples with [mset, sourcedb, solutions] for each item
        mset_sourcedb_solutions_ch = mset_ch.map { mset ->
            [
                mset,
                "${mset}/filtered_sky_model.txt",
                "${mset}/${params.ddecal.dd.sols}"
            ]
        }

        subtract_ch = SubtractSources( 
            calibrate_ch.done.collect(), 
            mset_sourcedb_solutions_ch, 
            params.ddecal.dd.subtract.parset, 
            params.ddecal.dd.incol, 
            params.ddecal.dd.outcol, 
            params.ddecal.dd.subtract.exclude_directions,
        )

        aoflagger_ch = AOFlag ( true, subtract_ch, params.ddecal.dd.outcol, params.postdd.aoflagger_strategy, 1, 0 )

        aoq_ch = AOqualityCollect( aoflagger_ch, mset_ch, params.ddecal.dd.outcol ) //

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.dd.outcol )

    emit:
        // AOqualityCollect.out
        WScleanImage.out.done

}



workflow PostDD {
    take:
        start_ch

    main:
        // mset_ch = channel.fromPath( "fullband_dd_smooth_data_T???.MS", glob: true, checkIfExists: true, type: 'dir' )
        mset_ch = channel.fromPath( params.data.dd_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

        filter_ch = FlagBaselines(start_ch, mset_ch, params.postdd.filter_baselines)

        beam_ch = ApplyBEAM( filter_ch.collect(), mset_ch, params.postdd.beam.parset, params.ddecal.dd.outcol, params.postdd.beam.outcol )

        uvwflag_ch = UVWFlag ( beam_ch.collect(), mset_ch, params.postdd.beam.outcol, params.postdd.uvlambdamin, params.postdd.uvlambdamax )

        flagged_msnames_ch = mset_ch.map { "${params.data.path}/" + it.getName().replace( ".MS", ".MS.l${params.postdd.uvlambdamin}to${params.postdd.uvlambdamax}" ) }

        aoq_ch = AOqualityCollect( uvwflag_ch.done.collect(), flagged_msnames_ch.flatten(), 'DATA' )

        // uvwflag_ch = UVWFlag (beam_ch.collect(), mset_ch, params.ddecal.dd.outcol, params.postdd.uvlambdamin, params.postdd.uvlambdamax)
        // beam_ch = ApplyBEAM( uvwflag_ch.done.collect(), flagged_msnames_ch, params.postdd.beam.parset, 'DATA', params.postdd.beam.outcol )
        // // im_names_ch =  flagged_msnames_ch.map { it.getSimpleName()  + "_" + params.wsclean.imname }

        im_names_ch = mset_ch.map { it.getSimpleName().replace( ".MS", ".MS.l${params.postdd.uvlambdamin}to${params.postdd.uvlambdamax}" ) +  "_" + params.wsclean.imname }

        mset_and_names_ch = flagged_msnames_ch.merge( im_names_ch.flatten() )

        WScleanImage ( aoq_ch.collect(), mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  'DATA' )

    emit:

        WScleanImage.out.done
        // UVWFlag.out.done.collect()

}


// workflow PostDD {
//     take:
//         start_ch

//     main:
//         mset_ch = channel.fromPath( params.data.dd_ms_glob, glob: true, checkIfExists: true, type: 'dir' )

//         aoflagger_ch = AOFlag ( start_ch, mset_ch, params.ddecal.dd.outcol, params.postdd.aoflagger_strategy, 1, 1 )

//         uvwflag_ch = UVWFlag (aoflagger_ch.collect(), mset_ch, params.ddecal.dd.outcol, params.postdd.uvlambdamin, params.postdd.uvlambdamax)

//         beam_ch = ApplyBEAM( uvwflag_ch.collect(), mset_ch, params.postdd.beam.parset, params.ddecal.dd.outcol, params.postdd.beam.outcol )

//         aoq_ch = AOqualityCollect( true, beam_ch, params.postdd.beam.outcol )

//         im_names_ch =  mset_ch.collect { it.getSimpleName()  + "_" + params.wsclean.imname }

//         mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

//         WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.postdd.beam.outcol )

//     emit:

//         WScleanImage.out.done

// }


workflow WS {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath( "${params.data.path}/${params.average.ditodd.msout}_T*flagged.MS", glob: true, checkIfExists: true, type: 'dir' )
        // mset_ch = channel.fromPath( "${params.data.path}/${params.ssplit.ms_prefix}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( start_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.wsclean.column )  //params.ddecal.beam.outcol ) // params.ddecal.dd.outcol

}


workflow TAR {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath( params.data.raw_tar_glob, glob: true, checkIfExists: true, type: 'file' )

        untar_ch = UnpackMSTarball( start_ch, mset_ch, params.data.label )


}
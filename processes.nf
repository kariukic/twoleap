#!/usr/bin/env nextflow

include {
    getTime
} from './modules/utils.nf'


process FlagIntra {
    debug true
    label 'sing'

    input:
    val ready
    path ms

    output:
    val true

    script:
    time = getTime()
    """
        python3 ${projectDir}/templates/flag_intrastations.py -i ${ms} > "${ms}/flag_intrastations_${time}.log" 2>&1
        """
}
//         """
//         // #!/usr/bin/env python3
//         // import casacore.tables as tab
//         // tab.taql(
//         //         r'UPDATE ${ms} SET FLAG=True WHERE mscal.baseline("/(.*)HBA0&\1HBA1/")'
//         //     )
//         """
// }

process ClipData {
    debug true
    label 'sing'

    input:
    path ms

    output:
    path "${ms}"

    script:
    time = getTime()
    """
        python3 ${projectDir}/templates/clip_data.py -i ${ms} --flag_intrastations --flag_badbaselines -c DATA  > "${ms}/clip_${time}.log" 2>&1
        """
}


process FlagBaselines {
    debug true
    label 'sing'

    input:
    val ready
    path ms
    val baselines

    output:
    val true

    script:
    time = getTime()
    """
        DP3 steps=[preflagger] msin=${ms} preflagger.baseline="${baselines}" msout=. msout.overwrite=True > "${params.out.logs}/filter_bad_stations.log"
        """
}


process FilterInter {
    debug true
    label 'singDPPP'
    publishDir "${params.data.path}", mode: 'move'

    input:
    val ready
    path ms

    output:
    path "${ms.getName() + '.noInter'}"

    script:
    time = getTime()
    """
        DP3 msin=${ms} steps=[filter] filter.remove=true filter.baseline="[CR]S*&&" msout="${ms.getName() + '.noInter'}" msout.overwrite=True > "${params.out.logs}/filter_intrastations_and_international_stations_${ms}_${time}.log" 2>&1
        """
}


process FitBpol {
    debug true
    publishDir "${ms}", mode: 'copy'

    input:
    val ready
    tuple path(ms), path(solsfile)
    val degree

    output:
    path "${solsfile.getSimpleName()}_degree${degree}_bpol.h5"

    script:
    time = getTime()
    """
        python3 ${projectDir}/templates/fit_bpol.py -s ${solsfile} -d ${degree}   > "${params.out.logs}/fit_degree${degree}_bpol_${ms}_${solsfile}_${time}.log" 2>&1
        """
}



process UnpackMSTarball {
    debug true
    label 'sing'
    publishDir "${params.data.path}", mode: 'move'

    input:
    val ready
    path ms
    val mslabel

    output:
    path "${ms.getSimpleName()}_${mslabel}.MS"

    script:
    time = getTime()
    """
        python3 ${projectDir}/templates/untarLTAData.py -i ${ms} -l ${mslabel}  > "${params.out.logs}/extract_${ms}_${time}.log" 2>&1
        """
}


process GetMSColumn {
    debug true
    label 'sing'

    input:
    val ready
    tuple path(ms), val(msout)
    val column

    output:
    path "${ms}"

    script:
    time = getTime()
    """
        DP3 steps=[] msin=${ms} msin.datacolumn=${column} msout=${msout} msout.overwrite=True > "${params.out.logs}/get_column_${ms}_${column}_${time}.log" 2>&
        """
}


process DP3GainCalDI {
    debug true
    label 'sing'
    // maxForks 4
    publishDir "${ms}", mode: 'copy'

    input:
    val ready
    path ms
    path parset
    path sourcedb
    val solsfile
    val incol
    val solint
    val uvlambdamin
    val uvlambdamax
    val uvmmax
    val nchan

    output:
    path "${solsfile}"

    script:

    time = getTime()

    """
        DP3 ${parset} msin=${ms} msin.datacolumn=${incol} gaincal.sourcedb=${sourcedb} gaincal.parmdb=${solsfile} gaincal.solint=${solint} gaincal.uvlambdamin=${uvlambdamin} gaincal.uvlambdamax=${uvlambdamax} gaincal.uvmmax=${uvmmax} gaincal.nchan=${nchan} > "${ms}/cal_${solsfile}_${time}.log"
        """
}


process DP3CalibrateDI {
    debug true
    label 'sing'
    maxForks 4
    publishDir "${ms}", mode: 'copy'

    input:
    val ready
    path ms
    path parset
    path sourcedb
    val solsfile
    val incol
    val solint
    val uvlambdamin
    val uvlambdamax
    val uvmmax
    val nchan
    val flagstations
    val calmode
    val smoothnessconstraint
    val propagate_sols
    val maxiter
    val beamproximitylimit
    val usebeam
    val beammode
    val propagate_converged_sols_only
    val scaling_coefficient

    output:
    path "${solsfile}"

    script:

    time = getTime()
    if (flagstations) {
        """
            DP3 ${parset} steps=[preflagger,scaledata,ddecal] msin=${ms} msin.datacolumn=${incol} preflagger.baseline="${flagstations}" scaledata.stations=[*] scaledata.coeffs=[${scaling_coefficient}] ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} ddecal.uvlambdamin=${uvlambdamin} ddecal.uvlambdamax=${uvlambdamax} ddecal.uvmmax=${uvmmax} ddecal.nchan=${nchan} ddecal.mode=${calmode} ddecal.smoothnessconstraint=${smoothnessconstraint} ddecal.propagatesolutions=${propagate_sols} ddecal.maxiter=${maxiter} ddecal.beamproximitylimit=${beamproximitylimit} ddecal.usebeammodel=${usebeam} ddecal.beammode=${beammode} ddecal.propagateconvergedonly=${propagate_converged_sols_only} > "${ms}/cal_${solsfile}_${time}.log"
            """
    }
    else {
        """
            DP3 ${parset} steps=[scaledata,ddecal] msin=${ms} msin.datacolumn=${incol} scaledata.stations=[*] scaledata.coeffs=[${scaling_coefficient}] ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} ddecal.uvlambdamin=${uvlambdamin} ddecal.uvlambdamax=${uvlambdamax} ddecal.uvmmax=${uvmmax} ddecal.nchan=${nchan} ddecal.mode=${calmode} ddecal.smoothnessconstraint=${smoothnessconstraint} ddecal.propagatesolutions=${propagate_sols} ddecal.maxiter=${maxiter} ddecal.beamproximitylimit=${beamproximitylimit} ddecal.usebeammodel=${usebeam} ddecal.beammode=${beammode} ddecal.propagateconvergedonly=${propagate_converged_sols_only} > "${ms}/cal_${solsfile}_${time}.log"
            """
    }
}

process DP3CalibrateDD {
    debug true
    label 'sing'
    maxForks 1
    publishDir "${ms}", mode: 'copy'

    input:
    val ready
    path ms
    path parset
    path sourcedb
    val solsfile
    val incol
    val calmode
    val solint
    val uvlambdamin
    val uvlambdamax
    val uvmmax
    val nchan
    val usebeam
    val beammode
    val smoothnessconstraint
    val truncateksmoothkernel
    val robust_reg
    val propagate_sols
    val maxiter
    val beamproximitylimit
    val correctfreqsmearing
    val flagstations
    val propagate_converged_sols_only
    val flux_threshold
    // Add this input
    val smoothness_max_factor
    // Add this input
    val solver

    output:
    path "${solsfile}", emit: solsfile
    path "filtered_sky_model.txt", emit: filtered_model
    path "smoothness_factors.csv"
    val true, emit: done

    script:
    time = getTime()

    // Create filtered sky model and get smoothness factors
    smoothness_command = """
        python3 ${projectDir}/templates/dd_smoothness_factors_v3.py \\
            ${ms} \\
            ${sourcedb} \\
            filtered_sky_model.txt \\
            --flux-threshold ${flux_threshold} \\
            --max-smoothness ${smoothness_max_factor} \\
            --exclude c3c196,cluster1 \\
            --exclude-values 1,1 \\
            --output-factors smoothness_factors.csv > "${ms}/filter_sources_${time}.log"
        """

    // Calculate solutions_per_direction
    // --min_flux 25 --max_flux 50 use this to force an output of all 1s, migth be useful for testing
    //  --min_flux 1 --max_flux 5, use this to get actual solutions per direction, works for Ateams subtraction
    solutions_per_direction_command = """
            python3 ${projectDir}/templates/sols_per_dir.py ${sourcedb} ${ms} --solint ${solint} --exclude_direction ${params.ddecal.dd.subtract.exclude_directions} --flux_threshold ${flux_threshold} --min_flux 25 --max_flux 50 > "${ms}/solperdir.txt"
        """

    // Extract smoothness factors from CSV file
    extract_factors_command = """
        # Read smoothness factors for passing patches only
        python3 -c "
        import csv
        factors = []
        with open('smoothness_factors.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['passes_threshold'] == 'YES':
                    factors.append(row['smoothness_factor'])
        # Print as comma-separated list for DP3
        print(','.join(factors))
        " > smoothness_factors_list.txt
        """

    if (flagstations) {
        """
            #step 0: get number of subsolutions per direction
            ${solutions_per_direction_command}
            solutions_per_direction=\$(cat solutions_per_direction.txt)

            # Step 1: Filter sources and calculate smoothness factors
            ${smoothness_command}
            
            # Step 2: Extract smoothness factors as comma-separated list
            ${extract_factors_command}
            
            # Step 3: Read smoothness factors
            smoothness_factors=\$(cat smoothness_factors_list.txt)
            
            # Step 4: Run calibration with filtered model and smoothness factors
            DP3 ${parset} steps=[preflagger,ddecal] \\
                msin=${ms} \\
                preflagger.baseline="${flagstations}" \\
                msin.datacolumn=${incol} \\
                ddecal.sourcedb=filtered_sky_model.txt \\
                ddecal.smoothness_dd_factors="[\${smoothness_factors}]" \\
                ddecal.solutions_per_direction="\${solutions_per_direction}" \\
                ddecal.h5parm=${solsfile} \\
                ddecal.solint=${solint} \\
                ddecal.uvlambdamin=${uvlambdamin} \\
                ddecal.uvlambdamax=${uvlambdamax} \\
                ddecal.uvmmax=${uvmmax} \\
                ddecal.nchan=${nchan} \\
                ddecal.mode=${calmode} \\
                ddecal.smoothnessconstraint=${smoothnessconstraint} \\
                ddecal.model_weighted_constraints=${robust_reg} \\
                ddecal.propagatesolutions=${propagate_sols} \\
                ddecal.maxiter=${maxiter} \\
                ddecal.beamproximitylimit=${beamproximitylimit} \\
                ddecal.correctfreqsmearing=${correctfreqsmearing} \\
                ddecal.usebeammodel=${usebeam} \\
                ddecal.beammode=${beammode} \\
                ddecal.smoothness_kernel_truncation=${truncateksmoothkernel} \\
                ddecal.solveralgorithm=${solver} \\
                ddecal.propagateconvergedonly=${propagate_converged_sols_only} > "${ms}/cal_${solsfile}_${time}.log"
            """
    }
    else {
        """
            # Step 1: Filter sources and calculate smoothness factors
            ${smoothness_command}
            
            # Step 2: Extract smoothness factors as comma-separated list
            ${extract_factors_command}
            
            # Step 3: Read smoothness factors
            smoothness_factors=\$(cat smoothness_factors_list.txt)
            
            # Step 4: Run calibration with filtered model and smoothness factors
            DP3 ${parset} steps=[ddecal] \\
                msin=${ms} \\
                msin.datacolumn=${incol} \\
                ddecal.sourcedb=filtered_sky_model.txt \\
                ddecal.smoothness_dd_factors="[\${smoothness_factors}]" \\
                ddecal.h5parm=${solsfile} \\
                ddecal.solint=${solint} \\
                ddecal.uvlambdamin=${uvlambdamin} \\
                ddecal.uvlambdamax=${uvlambdamax} \\
                ddecal.uvmmax=${uvmmax} \\
                ddecal.nchan=${nchan} \\
                ddecal.mode=${calmode} \\
                ddecal.smoothnessconstraint=${smoothnessconstraint} \\
                ddecal.solutions_per_direction="\${solutions_per_direction}" \\
                ddecal.model_weighted_constraints=${robust_reg} \\
                ddecal.propagatesolutions=${propagate_sols} \\
                ddecal.maxiter=${maxiter} \\
                ddecal.beamproximitylimit=${beamproximitylimit} \\
                ddecal.correctfreqsmearing=${correctfreqsmearing} \\
                ddecal.usebeammodel=${usebeam} \\
                ddecal.beammode=${beammode} \\
                ddecal.smoothness_kernel_truncation=${truncateksmoothkernel} \\
                ddecal.solveralgorithm=${solver} \\
                ddecal.propagateconvergedonly=${propagate_converged_sols_only} > "${ms}/cal_${solsfile}_${time}.log"
            """
    }
}


process DP3CalibrateDDD {
    debug true
    label 'singDPPP'
    maxForks 2
    publishDir "${ms}", mode: 'copy'

    input:
    val ready
    path ms
    path parset
    path sourcedb
    val solsfile
    val incol
    val calmode
    val solint
    val uvlambdamin
    val uvlambdamax
    val uvmmax
    val nchan
    val usebeam
    val beammode
    val smoothnessconstraint
    val truncateksmoothkernel
    val robust_reg
    val propagate_sols
    val maxiter
    val beamproximitylimit
    val correctfreqsmearing
    val flagstations
    val propagate_converged_sols_only

    output:
    path "${solsfile}"

    script:

    time = getTime()
    // chosen_sourcedb='ddmodel.txt'

    if (flagstations) {

        """
            # python3 /home/codex/chege/projects/3C196/notebooks_v2/ateam_model/add_cyg_to_model.py -m ${ms} > "${ms}/choosing_model_${time}.log"

            DP3 ${parset} steps=[preflagger,ddecal] msin=${ms} preflagger.baseline="${flagstations}" msin.datacolumn=${incol} ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} ddecal.uvlambdamin=${uvlambdamin} ddecal.uvlambdamax=${uvlambdamax} ddecal.uvmmax=${uvmmax} ddecal.nchan=${nchan} ddecal.mode=${calmode} ddecal.smoothnessconstraint=${smoothnessconstraint} ddecal.model_weighted_constraints=${robust_reg} ddecal.propagatesolutions=${propagate_sols} ddecal.maxiter=${maxiter} ddecal.beamproximitylimit=${beamproximitylimit} ddecal.correctfreqsmearing=${correctfreqsmearing} ddecal.usebeammodel=${usebeam} ddecal.beammode=${beammode} ddecal.smoothness_kernel_truncation=${truncateksmoothkernel} ddecal.propagateconvergedonly=${propagate_converged_sols_only} > "${ms}/cal_${solsfile}_${time}.log"
            """
    }
    else {

        """
            # python3 /home/codex/chege/projects/3C196/notebooks_v2/ateam_model/add_cyg_to_model.py -m ${ms} > "${ms}/choosing_model_${time}.log"

            DP3 ${parset} steps=[ddecal] msin=${ms} msin.datacolumn=${incol} ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} ddecal.uvlambdamin=${uvlambdamin} ddecal.uvlambdamax=${uvlambdamax} ddecal.uvmmax=${uvmmax} ddecal.nchan=${nchan} ddecal.mode=${calmode} ddecal.smoothnessconstraint=${smoothnessconstraint} ddecal.model_weighted_constraints=${robust_reg} ddecal.propagatesolutions=${propagate_sols} ddecal.maxiter=${maxiter} ddecal.beamproximitylimit=${beamproximitylimit} ddecal.correctfreqsmearing=${correctfreqsmearing} ddecal.usebeammodel=${usebeam} ddecal.beammode=${beammode} ddecal.smoothness_kernel_truncation=${truncateksmoothkernel} ddecal.propagateconvergedonly=${propagate_converged_sols_only} > "${ms}/cal_${solsfile}_${time}.log"
            """
    }
}


process ApplyGains {
    debug true
    label 'singDPPP'

    input:
    val ready
    tuple path(ms), path(solsfile)
    val parset
    val incol
    val outcol

    output:
    path "${ms}"

    script:

    time = getTime()

    """
        DP3 ${parset} msin=${ms} applycal.parmdb=${solsfile} msin.datacolumn=${incol} msout.datacolumn=${outcol} > "${ms}/applygains_${solsfile}_to_${outcol}_${time}.log" 2>&1
        """
}



//Subtract a sky direction(s). For subtracting specific directions use this in the shell instead of script
// '''
// #directions_to_subtract=$(<!{sources_to_subtract_file})
// #DP3 !{subtraction_parset} msin=!{full_ms_path} sub.applycal.parmdb=!{calibration_solutions_file} sub.sourcedb=!{sourcedb_name} sub.directions=${directions_to_subtract} msin.datacolumn=!{input_datacolumn} msout.datacolumn=!{output_datacolumn}> di_sub.log
// process SubtractSources {
//     label 'singDPPP'
//     input:
//         val ready
//         tuple path(full_ms_path), path(sourcedb_name), path(calibration_solutions_file) //, path(sources_to_subtract_file)
//         path subtraction_parset
//         val input_datacolumn
//         val output_datacolumn

//     output:
//         path "${full_ms_path}"

//     script:
//         time = getTime()
//         // chosen_sourcedb='ddmodel.txt'
//         """
//         # python3 /home/codex/chege/projects/3C196/notebooks_v2/ateam_model/add_cyg_to_model.py -m ${full_ms_path} > "${full_ms_path}/choosing_model_${time}.log"
//         DP3 ${subtraction_parset} msin=${full_ms_path} sub.applycal.parmdb=${calibration_solutions_file} sub.sourcedb=${sourcedb_name} msin.datacolumn=${input_datacolumn} msout.datacolumn=${output_datacolumn} > "${full_ms_path}/dd_subtract_${time}.log" 2>&1
//         """
// }


process SubtractSources {
    label 'singDPPP'

    input:
    val ready
    tuple path(full_ms_path), path(sourcedb_name), path(calibration_solutions_file)
    path subtraction_parset
    val input_datacolumn
    val output_datacolumn
    val exclude_directions

    output:
    path "${full_ms_path}"

    script:
    time = getTime()
    """
        # Get subtraction directions using Python script
        SUBTRACT_DIRECTIONS=\$(python3 ${projectDir}/templates/get_subtract_directions.py "${sourcedb_name}" "${exclude_directions}")
        echo "Subtraction directions: \$SUBTRACT_DIRECTIONS"
        
        # Run DP3 subtraction
        DP3 ${subtraction_parset} \\
            msin=${full_ms_path} \\
            sub.applycal.parmdb=${calibration_solutions_file} \\
            sub.sourcedb=${sourcedb_name} \\
            sub.directions=\${SUBTRACT_DIRECTIONS} \\
            msin.datacolumn=${input_datacolumn} \\
            msout.datacolumn=${output_datacolumn} > "${full_ms_path}/dd_subtract_${time}.log" 2>&1
        """
}


process MakeDP3ClustersListFile {
    input:
    val ready
    val num
    val fname

    output:
    path "${fname}"

    script:
    """
        #!/usr/bin/env python3
        clusters_str = ",".join([f"[cluster{c}]" for c in list(range(1, ${num}))])
        with open("${fname}", "w") as txt:
            txt.write( f"[{clusters_str}]" )
        """
}


process WScleanImage {
    label 'sing'
    publishDir "${params.data.path}/${params.out.results}/wsclean/${datacol}", pattern: "*.fits", mode: "move", overwrite: true
    publishDir "${params.data.path}/${params.out.results}/wsclean/${datacol}/plots", pattern: "*.png", mode: "move", overwrite: true

    input:
    val ready
    tuple val(mses), val(imname)
    val size
    val scale
    val niter
    val pol
    val chansout
    val minuvl
    val maxuvl
    val weight
    val polfit
    val datacol

    output:
    path "*.fits"
    val true, emit: done

    script:
    if (chansout == 1) {
        """
            wsclean -v -log-time -name ${imname} -data-column ${datacol} -pol ${pol} -weight ${weight} -scale ${scale} -size ${size} ${size} -make-psf -niter ${niter} -gridder wgridder -reorder ${mses} > ${params.out.logs}/wsclean_${imname}_image.log

            ls *-I-image.fits > imlist_I.txt
            ls *-V-image.fits > imlist_V.txt
            ls *-Q-image.fits > imlist_Q.txt
            ls *-U-image.fits > imlist_U.txt

            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_I.txt --filename ${imname}_I
            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_V.txt --filename ${imname}_V
            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_Q.txt --filename ${imname}_Q
            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_U.txt --filename ${imname}_U
            """
    }
    else {
        """
            wsclean -v -log-time -name ${imname} -data-column ${datacol} -pol ${pol} -weight ${weight} -scale ${scale} -size ${size} ${size} -niter ${niter} -apply-primary-beam -make-psf -join-channels -channels-out ${chansout} -gridder wgridder -no-update-model-required -no-dirty -no-mf-weighting ${mses} > ${params.out.logs}/wsclean_${imname}_image.log

            ls *-I-image.fits > imlist_I.txt
            ls *-V-image.fits > imlist_V.txt
            ls *-Q-image.fits > imlist_Q.txt
            ls *-U-image.fits > imlist_U.txt

            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_I.txt --filename ${imname}_I
            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_V.txt --filename ${imname}_V
            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_Q.txt --filename ${imname}_Q
            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist_U.txt --filename ${imname}_U

            """
    }
}

//   -minuv-l ${minuvl} -maxuv-l ${maxuvl}
//   -minuv-l ${minuvl} -maxuv-l ${maxuvl}
// -save-source-list -fit-spectral-pol !{spectral_pol_fit} -multiscale -no-update-model-required -auto-mask 3 -auto-threshold 1 -mgain 0.6 -local-rms


// Collect data quality statistics
process AOqualityCollect {
    label 'sing'
    maxForks 5

    input:
    val ready
    path full_ms_path
    val data_column

    output:
    val true

    script:
    time = getTime()
    """
        aoquality collect -d ${data_column} ${full_ms_path}  > ${params.out.logs}/aoq_collect_${full_ms_path.getName()}_${time}.log 2>&1
        """
}


process GetData {
    debug true

    input:
    val ready
    val nodes
    val glob
    val txtname

    output:
    path "${txtname}"

    script:
    """
        ls -d ${glob} >> ${txtname}
        cp ${txtname} ${params.data.path}
        """
}


// process AOqualityCombine {
//     label 'singDPPP'

//     input:
//         val ready
//         val mses
//         val output_name

//     output:
//         val true, emit: qstats

//     shell:
//         '''
//         mkdir -p !{params.data.path}/!{params.out.results}/aoquality
//         aoquality combine !{params.data.path}/!{params.out.results}/aoquality/!{output_name}.qs !{mses} > !{params.out.logs}/aoquality_combine_!{output_name}.log
//         # python3 !{projectDir}/templates/plot_aoqstats.py -q !{params.data.path}/!{params.out.results}/aoquality/!{output_name}.qs -o !{params.data.path}/!{params.out.results}/aoquality/!{output_name}.png >> !{params.out.logs}/aoquality_combine.log
//         '''
// }

process AOqualityCombine {
    label 'sing'
    publishDir "${params.data.path}/${params.out.results}/aoquality/${output_name}", pattern: "*.qs", mode: "copy", overwrite: true
    publishDir "${params.data.path}/${params.out.results}/aoquality/${output_name}", pattern: "*.pkl", mode: "copy", overwrite: true
    publishDir "${params.data.path}/${params.out.results}/aoquality/${output_name}/plots", pattern: "*.pdf", mode: "copy", overwrite: true

    input:
    val ready
    val file_list
    val output_name

    output:
    path "${output_name}.qs"
    path "*.pdf"
    val true, emit: done

    script:
    time = getTime()
    def tlist: List = file(file_list).readLines()
    def mses: String = tlist.collect { "${it}" }.join(" ")
    """
        aoquality combine ${output_name}.qs ${mses} > ${params.out.logs}/combine_${output_name}_${time}.log 2>&1
        python3 ${projectDir}/templates/plot_flags.py plot_occ ${file_list} --filename ${output_name} >> ${params.out.logs}/combine_${output_name}_${time}.log 2>&1
        python3 ${projectDir}/templates/plot_aoqstats.py plot_aoq "${output_name}.qs" --name ${output_name} >> ${params.out.logs}/combine_${output_name}_${time}.log 2>&1
        """
}


process H5ParmCollect {
    debug true
    label 'sing'

    publishDir "${params.data.path}/${params.out.results}/solutions/${output_name}", pattern: "*.h5", mode: "move", overwrite: true
    publishDir "${params.data.path}/${params.out.results}/solutions/${output_name}/plots", pattern: "*.png", mode: "move", overwrite: true

    input:
    val ready
    val solution_files
    val output_name

    output:
    path "${output_name}.h5", emit: combined_sols
    path "*.png", emit: plots
    val true, emit: done

    script:
    """
        H5parm_collector.py ${solution_files} -o ${output_name}.h5 > "${params.out.logs}/h5parm_collect.log" 2>&1
        #soltool plot --plot_dir \$(pwd) ${output_name}.h5 >> "${params.out.logs}/h5parm_collect.log" 2>&1
        python3 ${projectDir}/templates/plot_sols.py plot ${output_name}.h5 >> "${params.out.logs}/h5parm_collect.log" 2>&1
        python3 ${projectDir}/templates/plot_gains.py plot ${output_name}.h5 >> "${params.out.logs}/h5parm_collect.log" 2>&1
        """
}


process MergeChansSplitTime {
    label 'pspipe'
    maxForks 1

    input:
    val ready
    val msfiles
    val nodes
    val ntimes
    val column
    val msout
    val mses_per_node
    val outtxt

    output:
    val true

    shell:
    nd = nodes.join(' ')
    """
        python3 !{projectDir}/templates/concat_split.py --mslist !{msfiles} --msout !{msout} --ntimes !{ntimes} --nodes !{nd} --datapath !{params.data.path} --datacolumn !{column} --output_ms_list_file !{outtxt} --nmses_per_node !{mses_per_node} -t 0.1 > ${params.out.logs}/concat_freqs_split_time.log 2>&1
        """
}


process SplitMSToSubbands {
    label 'sing'

    input:
    val ready
    val msfile
    val nchans_per_msout
    val datacolumn

    output:
    val true

    shell:
    """
        python3 ${projectDir}/templates/split_mschans.py -m !{msfile} -n !{nchans_per_msout} -d !{datacolumn} > ${msfile}/split_mschans_to_subbands.log 2>&1
        """
}


process GetTimeChunksPerSubband {
    label 'sing'

    input:
    val ready
    val msin
    val nmses_per_node
    val from_nodes
    val to_nodes

    output:
    val true

    shell:
    fnds = from_nodes.join(' ')
    tnds = to_nodes.join(' ')
    """
        python3 ${projectDir}/templates/write_subband_time_chunks.py -m !{msin} -d !{params.data.path} -n !{nmses_per_node} -f !{fnds} -t !{tnds} > ${params.out.logs}/write_subband_time_chunks.log 2>&1
        """
}


process ConcatMSesinTime {
    label 'sing'
    publishDir "${params.data.path}", mode: "move", overwrite: true

    input:
    val ready
    tuple val(msfiles), val(msout)

    output:
    path "${msout}"
    val true, emit: done

    shell:
    """
        python3 ${projectDir}/templates/concatenate_msfiles.py !{msfiles} --msout !{msout} --concat_property time > ${params.out.logs}/!{msout}_concat_mses_in_time.log 2>&1
        """
}


process WriteMSlist {
    publishDir params.out.logs, mode: "copy", overwrite: true

    input:
    val ready
    val nodes
    val glob_pattern
    val txtname

    output:
    path "${txtname}", emit: per_line_mslist
    path "${txtname}.ps", emit: single_line_mslist

    script:
    """
#!/usr/bin/python3
from glob import glob
mses=[]
for node in ${nodes}:
    mslist=glob(f"/net/node{node}/${glob_pattern}")
    mses+=mslist
with open("${txtname}", "w") as out:
    for ms in sorted(mses):
        out.write(f"{ms}\\n")
with open("${txtname}.ps", "w") as out:
        out.write(f"{' '.join(sorted(mses))}")
    """
}


// process WriteDDMSlist {
//     input:
//         val ready
//         val nodes

//     output:
//         val true

//     script:
//     """
// #!/usr/bin/python3
// from glob import glob
// mses=[]
// for node in ${nodes}:
//     mslist=glob(f"/net/node{node}/${params.data.ms_files.dd}")
//     mses+=mslist
// with open("${params.out.logs}/${params.data.dd_mslist}", "w") as out:
//     for ms in mses:
//         out.write(f"{ms}\\n")
// with open("${params.out.logs}/${params.data.dd_mslist}.ps", "w") as out:
//         out.write(f"{' '.join(mses)}")
//     """
// }


process ApplyBEAM {
    label 'singDPPP'

    input:
    val ready
    path ms
    path parset
    val incol
    val outcol

    output:
    // path "${ms}"
    val true

    script:
    time = getTime()

    """
        DP3 ${parset} msin=${ms} msin.datacolumn=${incol} msout.datacolumn=${outcol} > "${ms}/apply_beam_${time}.log" 2>&1
        """
}

process ReadTxtLinesandAppend {
    input:
    val ready
    val dirname
    val txtname
    val postfix

    output:
    val ms_string, emit: list_str
    val ms_postfix_string, emit: list_postfix_str

    exec:
    tlist = file(dirname).resolve(txtname).readLines()

    ms_string = tlist.collect { "${it}" }.join(" ")

    ms_postfix_string = tlist.collect { "${it}" + postfix }.join(" ")
}


process AOFlag {

    debug true
    label 'singDPPP'
    publishDir "${params.data.path}", mode: 'move'
    maxForks 6

    input:
    val ready
    path ms
    val column
    path aoflagger_strategy
    val aoflag
    val interpolate

    output:
    val true, emit: done

    script:
    // time=getTime()

    // if ( interpolate == 1 )
    //     """
    //     DP3 steps=[aoflag,interpolate] msin=${ms} msin.datacolumn=${column} aoflag.type=aoflagger aoflag.strategy=${aoflagger_strategy} msout=. msout.overwrite=True > "${params.out.logs}/flag_${ms}_${column}_${time}.log" 2>&1
    //     """
    // else
    //     """
    //     DP3 steps=[aoflag] msin=${ms} msin.datacolumn=${column} aoflag.type=aoflagger aoflag.strategy=${aoflagger_strategy} msout=. msout.overwrite=True > "${params.out.logs}/flag_${ms}_${column}_${time}.log" 2>&1
    //     """

    time = getTime()

    // Build steps list based on selections
    steps = []
    if (aoflag == 1) {
        steps << "aoflag"
    }
    if (interpolate == 1) {
        steps << "interpolate"
    }

    // Convert to string representation for DP3
    steps_str = "[" + steps.join(",") + "]"

    // Build step-specific parameters
    step_params = []
    if (aoflag == 1) {
        step_params << "aoflag.type=aoflagger"
        step_params << "aoflag.strategy=${aoflagger_strategy}"
    }
    if (interpolate == 1) {
        step_params << "interpolate.type=interpolate"
    }

    step_params_str = step_params.join(" ")

    """
        DP3 steps=${steps_str} msin=${ms} msin.datacolumn=${column} ${step_params_str} msout=. msout.overwrite=True > "${params.out.logs}/flag_interp_${ms}_${column}_${time}.log" 2>&1
        """
}

process UVWFlag {

    debug true
    label 'sing'
    publishDir "${params.data.path}", mode: 'move'

    input:
    val ready
    path msin
    val data_column
    val uvlambdamin
    val uvlambdamax

    output:
    path "${msin}.l${uvlambdamin}to${uvlambdamax}", emit: msout
    val true, emit: done

    script:
    time = getTime()
    """
        DP3 steps=[uvwflag] msin=${msin} msin.datacolumn=${data_column} uvwflag.uvlambdamin=${uvlambdamin} uvwflag.uvlambdamax=${uvlambdamax} msout=${msin}.l${uvlambdamin}to${uvlambdamax} msout.overwrite=True > "${params.out.logs}/uvwflag_${msin}_${data_column}_${time}.log" 2>&1
        """
}



process Compress {
    debug true
    label 'singDPPP'
    maxForks 3
    publishDir "${params.data.path}", mode: 'move'

    input:
    path ms
    val nbits
    val normalization
    val distribution
    val disttruncation

    output:
    // path "${ms.getSimpleName()}_DC.MS"
    path "${ms.getName().replace(".MS", ".DCMS")}"

    script:
    time = getTime()
    """
        DP3 steps=[aoflag,interpolate] msin=${ms} msin.datacolumn=DATA aoflag.type=aoflagger aoflag.memoryperc=20 msout="${ms.getName().replace(".MS", ".DCMS")}" msout.storagemanager=dysco msout.storagemanager.databitrate=${nbits} msout.storagemanager.distribution=${distribution} msout.storagemanager.normalization=${normalization} msout.storagemanager.disttruncation=${disttruncation} msout.overwrite=True > "${ms}/compress_after_flagging_col_${time}.log" 2>&1
        """
}

process Demix {
    debug true
    label 'sing_demix'
    maxForks 3
    publishDir "${params.data.path}", mode: 'move'

    input:
    val enabled
    path msin
    path parset
    path sourcedb

    output:
    path "${msin.getName().replace(".FMS", ".FDMS")}"

    script:
    time = getTime()

    if (enabled == 1) {
        """
            DP3 ${parset} msin=${msin} demix.skymodel=${sourcedb} msout=${msin.getName().replace(".FMS", ".FDMS")} msout.overwrite=true > "${msin}/demix_${time}.log" 2>&1
            """
    }
    else {
        """
            DP3 msin=${msin} steps=[] msout=${msin.getName().replace(".FMS", ".FDMS")} msout.overwrite=true > "${msin}/demix_disabled_${time}.log" 2>&1
            """
    }
}



// process AverageDItoDDMS {
//     label 'singDPPP'
//     publishDir "${params.data.path}", mode: 'move'

//     input:
//         val ready
//         path ms002
//         val data_column
//         val timestep
//         val freqstep


//     output:
//         path "*_SAP*_SB*_uv_003*.MS"
//         val true , emit: done_averaging

//     shell:
//         '''
//         ms003=$(echo "!{ms002}" | sed "s/002/003/")
//         DP3 steps=[avg] msin=!{ms002} msin.datacolumn=!{data_column} msout=${ms003} avg.type=average avg.timestep=!{timestep} avg.freqstep=!{freqstep}
//         '''
// }


process Average {
    label 'sing'
    publishDir params.data.path, mode: 'move'
    maxForks 6

    input:
    val ready
    tuple path(msin), val(msout)
    val data_column
    val timestep
    val freqstep

    output:
    path "${msout}", emit: averaged_ms
    val true, emit: done_averaging

    script:
    time = getTime()
    """
        python3 ${projectDir}/templates/fill_flagged_weights.py ${msin} --mode nn > "${params.out.logs}/fix_flagged_weights_${msin}_${time}.log" 2>&1
        DP3 steps=[avg] msin=${msin} msin.datacolumn=${data_column} msout=${msout} avg.type=average avg.timestep=${timestep} avg.freqstep=${freqstep} msout.overwrite=True > "${params.out.logs}/average_${msin}_${timestep}tstep_${freqstep}freqstep_${time}.log"
        """
}

//taql update ${msin} set WEIGHT_SPECTRUM=1 
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REMOVE TAQL BEFORE AVERAGING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//Also flags intrastation baselines
process FilterBaselinesAndAverage {
    label 'singDPPP'
    publishDir params.data.path, mode: 'move'

    input:
    val ready
    tuple path(msin), val(msout)
    val data_column
    val timestep
    val freqstep

    output:
    path "${msout}", emit: filtered_averaged_ms
    val true, emit: done

    script:
    time = getTime()
    """
        DP3 steps=[filter,avg] msin=${msin} msin.datacolumn=${data_column} filter.remove=true filter.baseline="[CR]S*&&" avg.type=average avg.timestep=${timestep} avg.freqstep=${freqstep} msout=${msout} msout.overwrite=True > "${params.out.logs}/average_${msin}_${timestep}tstep_${freqstep}freqstep_${time}.log"
        """
}


process MakeClusters {
    label 'sing'

    input:
    path input_model
    val number_of_clusters
    val output_model

    output:
    path "${output_model}"

    shell:
    """
        cluster !{input_model} !{output_model} !{number_of_clusters}
        """
}


// process ClipGains {

//     input:
//         path solsfile
//         val nsigma
//         val mode

//     output:
//         val true

//         script:
//     """
// #!/usr/bin/python3
// from losoto.h5parm import h5parm
// H = h5parm('file.h5', readonly=False)
// soltab = H.getSolset('solset000').getSoltab('soltab000')
//     """

// }


def readTxtIntoString(txt) {
    def tlist: List = file(txt).readLines()
    def tstring: String = tlist.collect { "${it}" }.join(" ")

    return tstring
}


def readTxtAndAppendString(txt, str) {
    def tlist: List = file(txt).readLines()
    def tstring: String = tlist.collect { "${it}" + str }.join(" ")

    return tstring
}



/*
pssh needs a file listing the nodes to run commands on
This command makes sucha a file given a nodeslist and a file path.
input:      list of strings e.g. ["node100", "node101"]
            file name string e.g "hosts_list.txt"
output:     The written file
*/
def writeHosts(nodes_list, hosts) {
    hosts = new File(hosts)
    if (hosts.exists()) {
        hosts.delete()
    }
    hosts.createNewFile()
    hosts.withWriter { out ->
        nodes_list.each {
            out.println("${it}")
        }
    }
}


/*
return a list of the nodes with the 'node' prefix given an input string e.g. [node129]
also assign the master node and number of processes if not provided.
*/
def parseNodes(nodes) {
    // log.info """[nextleap*] init > verifying nodes list"""

    if (nodes instanceof String) {
        nodes_list = nodes.split(',').collect { "node${it}" } as List
    }
    else if (nodes instanceof Integer) {
        nodes_list = nodes.collect { "node${it}" } as List
    }
    else {
        log.error("Error: The `data.nodes` parameter is not valid. Got `--data.nodes=${nodes}`")
        exit(0)
    }

    return nodes_list
}


def makeDirectory(fileName) {
    def file = new File(fileName)
    if (!file.exists()) {
        file.mkdir()
    }
}

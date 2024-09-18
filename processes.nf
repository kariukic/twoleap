#!/usr/bin/env nextflow

include {
    getTime;
} from './modules/utils.nf'


process ScaleData {
    debug true

    input:
    path ms

    output: 
    path "${ms}"

    script:
        time = getTime()
        """
        python3 /home/codex/chege/software/pipelines/nextleap/scaledata.py -i ${ms} -f 0 > "${ms}/scale_${time}.log" 2>&1
        """
}


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
        python3 /home/codex/chege/software/pipelines/nextleap/templates/clip_data.py -i ${ms} --flag_intrastations --flag_badbaselines -c DATA  > "${ms}/clip_${time}.log" 2>&1
        """
}


process DP3Calibrate {
    debug true
    label 'sing'
    maxForks 1
    publishDir "${ms}" , mode: 'copy'

    input:
        val ready
        path ms
        path parset
        path sourcedb
        val solsfile
        val incol
        val solint
        // val maxforks

    output:
        path "${solsfile}"

    script:

        time = getTime()

        """
        DP3 ${parset} msin=${ms} msin.datacolumn=${incol} ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} > "${ms}/cal_${solsfile}_${time}.log"
        """
}


process ApplyGains {
    debug true
    label 'sing'

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



//Subtract a sky direction(s)
process SubtractSources {
    label 'sing'
    input:
        val ready
        tuple path(full_ms_path), path(sourcedb_name), path(calibration_solutions_file), path(sources_to_subtract_file)
        path subtraction_parset
        val input_datacolumn
        val output_datacolumn

    output:
        path "${full_ms_path}"
    
    shell:
        '''
        directions_to_subtract=$(<!{sources_to_subtract_file})
        DP3 !{subtraction_parset} msin=!{full_ms_path} sub.applycal.parmdb=!{calibration_solutions_file} sub.sourcedb=!{sourcedb_name} sub.directions=${directions_to_subtract} msin.datacolumn=!{input_datacolumn} msout.datacolumn=!{output_datacolumn}> di_sub.log
        '''
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
    publishDir "${params.data.path}/${params.out.results}/images", pattern: "*.fits", mode: "move", overwrite: true
    // publishDir "${params.data.path}/${params.out.results}/images", pattern: "*.txt", mode: "copy", overwrite: true

    input:
        val ready
        tuple val(mses), val(image_name)
        val size
        val scale
        val spectral_pol_fit
        val data_column
        

    output:
        path "*.fits"
        val true , emit: done
        // path "${image_name}-sources.txt", emit: model

    shell:
        '''
        wsclean -name !{image_name} -data-column !{data_column} -pol IV -weight briggs -0.1  -minuv-l 50 -maxuv-l 300 -scale !{scale} -size !{size} !{size} -make-psf -niter 0 -join-channels -channels-out 69 -gridder wgridder -wgridder-accuracy 1e-5 -reorder -no-mf-weighting !{mses} > !{params.out.logs}/!{image_name}_wsclean_image.log
        '''
}

// -save-source-list -fit-spectral-pol !{spectral_pol_fit} -multiscale -no-update-model-required -auto-mask 3 -auto-threshold 1 -mgain 0.6 -local-rms
// Collect data quality statistics
process AOqualityCollect {
    label 'sing'

    input:
        val ready
        path full_ms_path
        val data_column

    output:
        path "${full_ms_path}"

    shell:
        """
        aoquality collect -d !{data_column} !{full_ms_path}
        """
}


process AOqualityCombine {
    label 'sing'
    
    input:
        val ready
        val mses
        val output_name

    output:
        val true, emit: qstats

    shell:
        '''
        mkdir -p !{params.data.path}/!{params.out.results}/aoquality
        aoquality combine !{params.data.path}/!{params.out.results}/aoquality/!{output_name}.qs !{mses} > !{params.out.logs}/aoquality_combine.log
        # python3 !{projectDir}/templates/plot_aoqstats.py -q !{params.data.path}/!{params.out.results}/aoquality/!{output_name}.qs -o !{params.data.path}/!{params.out.results}/aoquality/!{output_name}.png >> !{params.out.logs}/aoquality_combine.log
        '''
}


process H5ParmCollect {

    publishDir "${params.data.path}/${params.out.results}/solutions/${output_name}", pattern: "*.h5", mode: "move", overwrite: true
    publishDir "${params.data.path}/${params.out.results}/solutions/${output_name}", pattern: "*.png", mode: "move", overwrite: true

    input:
        val ready
        val solution_files
        val output_name

    output:
        path "${output_name}.h5", emit: combined_sols
        path "*.png"

    shell:
        '''
        H5parm_collector.py !{solution_files} -o !{output_name}.h5 > !{params.out.logs}/h5parm_collect.log
        soltool plot --plot_dir $(pwd) !{output_name}.h5 >> !{params.out.logs}/h5parm_collect.log
        '''
}


process ConcatFrequencySplitTime {
    label 'sing'
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
        val true //"${outtxt}"

    shell:
        nd  = nodes.join(' ')
        """
        python3 !{projectDir}/templates/concat_split.py --mslist !{msfiles} --msout !{msout} --ntimes !{ntimes} --nodes !{nd} --datapath !{params.data.path} --datacolumn !{column} --output_ms_list_file !{outtxt} --nmses_per_node !{mses_per_node} > ${params.out.logs}/freq_concat.log 2>&1
        """

}

process WriteDDMSlist {
    input:
        val ready
        val nodes

    output:
        val true
    
    script:
    """
#!/usr/bin/python3
from glob import glob
mses=[]
for node in ${nodes}:
    mslist=glob(f"/net/node{node}/${params.data.path}/${params.average.ditodd.msout}_T*.MS")
    mses+=mslist
with open("${params.out.logs}/dd_mses.txt", "w") as out:
    for ms in mses:
        out.write(f"{ms}\\n")
    """
}


process ApplyBEAM {
    label 'sing'

    input:
        val ready
        path ms
        path parset
        val incol
        val outcol

    output:
        path "${ms}"

    script:
        time = getTime()

        """
        DP3 ${parset} msin=${ms} msin.datacolumn=${incol} msout.datacolumn=${outcol} > "${params.out.logs}/${ms}_apply_beam_${time}.log" 2>&1
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
        tlist = file( dirname ).resolve( txtname ).readLines()

        ms_string = tlist.collect {"${it}"}.join(" ")

        ms_postfix_string = tlist.collect {"${it}" + postfix}.join(" ")

}


process Flag {

    debug true
    label 'sing'
    publishDir "${params.data.path}", mode: 'move'

    input:
        path ms
        val column

    output:
        path "${ms.getSimpleName()}_flagged.MS"
    

    script:
        time=getTime()
        """
        DP3 steps=[aoflag,interpolate] msin=${ms} msin.datacolumn=${column} aoflag.type=aoflagger aoflag.memoryperc=20 "msout=${ms.getSimpleName()}_flagged.MS" > "${ms}/flag_${column}_${time}.log" 2>&1
        """
} 




process AverageDItoDDMS {
    label 'sing'
    publishDir "${params.data.path}", mode: 'move'

    input:
        val ready
        path ms002
        val data_column
        val timestep
        val freqstep


    output:
        path "*_SAP*_SB*_uv_003*.MS"
        val true , emit: done_averaging

    shell:
    '''
    ms003=$(echo "!{ms002}" | sed "s/002/003/")
    DP3 steps=[avg] msin=!{ms002} msin.datacolumn="!{data_column}" msout=${ms003} avg.type=average avg.timestep=!{timestep} avg.freqstep=!{freqstep}
    '''
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

ClipGains {

    input:
        path solsfile
        val nsigma
        val mode
    
    output:
        val true

        script:
    """
#!/usr/bin/python3
from losoto.h5parm import h5parm
H = h5parm('file.h5', readonly=False)
soltab = H.getSolset('solset000').getSoltab('soltab000')
    """

}

def readTxtIntoString (txt) {
    List tlist = file(txt).readLines()
    String tstring = tlist.collect {"${it}"}.join(" ")

    return tstring
}

def readTxtAndAppendString (txt, str) {
    List tlist = file(txt).readLines()
    String tstring = tlist.collect {"${it}" + str}.join(" ")

    return tstring
}

/*
pssh needs a file listing the nodes to run commands on
This command makes sucha a file given a nodeslist and a file path.
input:      list of strings e.g. ["node100", "node101"]
            file name string e.g "hosts_list.txt"
output:     The written file
*/
def writeHosts ( nodes_list, hosts ) {
    hosts = new File( hosts)
    if (hosts.exists()){
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
def parseNodes ( nodes ) {
    // log.info """[nextleap*] init > verifying nodes list"""

    if (nodes instanceof String){
        nodes_list = nodes.split(',').collect{"node${it}"} as List
    }
    // when a single node is given..
    else if (nodes instanceof Integer) {
        nodes_list = nodes.collect{"node${it}"} as List
    }
    else {
        log.error("Error: The `data.nodes` parameter is not valid. Got `--data.nodes=${nodes}`")
        exit 0
    }

    return nodes_list
}


def makeDirectory ( fileName ) {
    def file = new File(fileName)
    if(!file.exists()) {
        file.mkdir()
    }
}
#!/usr/bin/env nextflow

include {
    BP;
} from "./stages.nf"

workflow {
   BP ( true )
}



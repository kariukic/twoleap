#!/usr/bin/env nextflow

def get_input_ms_list(msfiles) -> void {
    def ms: String = 'ms'
    def txt: String = 'txt'

    def f = file(msfiles)
    def ext: String = f.getExtension().toLowerCase()

    if (ext == txt) {
        if (f.exists()) {
            log.info("Reading MS files from txt file: ${msfiles}")
            def msetsList = new File(msfiles).collect { it }
            return msetsList
        }
        else {
            log.error("${msfiles} does not exist")
        }
    }
    else if (ext == ms) {

        if (f.isDirectory()) {
            log.info("Got MS file: ${msfiles}")
            def msetsList = [msfiles]
            return msetsList
        }
        else {
            log.info("Globbing MS files with pattern: ${msfiles}")
            def msetsList = file(msfiles, glob: true, checkIfExists: true, type: 'dir').collect { it }
            return msetsList
        }
    }
    else {
        log.error("${msfiles} not understood")
    }
}


def get_input_h5_list(solsfiles) -> void {

    def h5: String = 'h5'
    def txt: String = 'txt'

    def f = file(solsfiles)
    def ext: String = f.getExtension().toLowerCase()

    if (ext == txt) {
        if (f.exists()) {
            log.info("Reading solsfiles from txt file: ${solsfiles}")
            def msetsList = new File(solsfiles).collect { it }
            return msetsList
        }
        else {
            log.error("${solsfiles} does not exist")
        }
    }
    else if (ext == h5) {

        if (f.isFile()) {
            log.info("Got h5 file: ${solsfiles}")
            def msetsList = [solsfiles]
            return msetsList
        }
        else {
            log.info("Globbing solsfiles with pattern: ${solsfiles}")
            def msetsList = file(solsfiles, glob: true, checkIfExists: true, type: 'dir').collect { it }
            return msetsList
        }
    }
    else {
        log.error("${solsfiles} not understood")
    }
}


// Get string-formatted current time
def getTime() {
    def date = new Date()
    def sdf = new java.text.SimpleDateFormat("dd_MM_yyyy_HH_mm_ss")
    def start_time = sdf.format(date)
    return start_time
}

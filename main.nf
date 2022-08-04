nextflow.enable.dsl=2

process copySamples {
input: 
	path(samplesCsv)
output: 
	path("samples.out.csv")
publishDir "${params.outDir}", mode: 'copy'

"""
	cat samplesCsv > samples.out.csv
"""
}

workflow {
	// Inputs
	samplesCsv = file(params.samples)
	samples = Channel.fromPath(samplesCsv).splitCsv(header:true, strip:true)
	copySamples(samplesCsv)
}

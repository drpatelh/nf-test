nextflow.enable.dsl=2

// Run bcParser to extract and correct cell-barcodes
// Optionally demuxes fastq files based on some barcodes
process barcodeDemux {
	publishDir "${params.outDir}/demux/", pattern: "$outDir/*gz"
	publishDir "${params.outDir}/demux/", mode: 'copy', pattern: "$outDir/*{txt,tsv,json}" //, saveAs: {"${fqName}.${it.getName()}"}

	input:
	path sheet // Samples CSV
	path libDir // Directory containing the library definition file (barcode sequence lists are loaded from here)
	val libName // Filename of the library definition file
	tuple val(fqName), path(fqFiles) // Input fastq file
	
	output:
	tuple val(fqName), path("$outDir/*_S[1-9]*.fastq.gz"), emit: fastq
	tuple val(fqName), path("$outDir/metrics.json")      , emit: metrics
	path "$outDir/*_S0_*.fastq.gz"                       , emit: unknown, optional: true
	path "$outDir/*.tsv"                                 , emit: tsv

	script:
	outDir = "${fqName}.demux"
	def lib = "$libDir/$libName"
	"""
	bc_parser \\
		--library $lib \\
		--demux $sheet \\
		--fastq-name $fqName \\
		-v \\
		--reads ${fqFiles.join(" ")} \\
		--write-fastq \\
		--out $outDir
	"""
}

process sampleReport {
	publishDir "$params.outDir", mode: 'copy'

	input: 
	val sample
	path samplesheet
	path libJson
	
	output:
	path "${sample}.csv"

	script:
	"""
	cat $samplesheet $libJson > ${sample}.csv
	"""
}

workflow inputReads {
	take:
	samplesCsv // samples.csv file
	libJson // library definition json
	fqDir // Path to directory with input fastqs; empty if running from BCL 

	main:
	// Get list of FastQ files
	fqDir
		.filter { it.name =~ /.fastq.gz$/ }
		.set { fqs }

	// Group fastq files by sample
	fqs
		.map { [ it.getName().toString().tokenize('_')[0], it ] }
		.groupTuple()
		.set { fqFiles }

	Channel
		.fromPath(samplesCsv)
		.splitCsv(header:true, strip:true)
		.map({ it.fastqName })
		.unique()
		.join(fqFiles)
		.set { fqSamples }

	// Process cell-barcodes and (optionally) split fastqs into samples based on tagmentation barcode
	barcodeDemux(samplesCsv, libJson.getParent(), libJson.getName(), fqSamples)

	barcodeDemux
		.out
		.fastq
		.flatMap({it[1]})
		.map { [ it.getName().toString().tokenize('_')[0], it ] }
		.groupTuple(size:2)
		.set { demuxFqs }

	emit:
	fqs = demuxFqs
	metrics = barcodeDemux.out.metrics
}

workflow {
	// Initialise params
	samplesCsv = file(params.samples)
	libJson = file(params.libStructure)
	fqDir = Channel.fromPath("${params.fastqDir}/*", checkIfExists:true) | take (params.fqcount)

	// Run inputReads subworkflow
	inputReads(samplesCsv, libJson, fqDir)
	
	// Run sampleReport subworkflow
	inputReads
		.out
		.fqs
		.map { it[0] }
		.set { foo }

	foo.view { "Foo: $it" }

	sampleReport(foo, samplesCsv, libJson)
}

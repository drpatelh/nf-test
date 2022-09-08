nextflow.enable.dsl=2
import groovy.json.JsonSlurper

def loadJson(json) {
	def jsonSlurper  = new JsonSlurper()
	return jsonSlurper.parse( json )
}

// Create a 'file()' from string 'path'
// 'path' could be null, a s3 URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	if (path =~ /^s3/) { return file(path)}
	return baseDir.resolve(path)
}

// Return 0 for null or non-integer strings
def toIntOr0(str) {
	if ((str != null) && str.isInteger())
		return str as int;
	else
		return 0
}

// Reference genome files and parameters
// Paths can be absolute or relative to location of json
def loadGenome(json) {
	def baseDir = json.getParent()
	genome = loadJson(json)
	genome.bowtie_index = expandPath(genome.bowtie_index, baseDir)
	genome.gtf = expandPath(genome.gtf, baseDir)
	genome.tss = expandPath(genome.tss, baseDir)
	return genome
}

// Load per-sample read-counts from bcParser output
def loadDemuxReadCounts(demuxMetrics) {
	def jsonSlurper  = new JsonSlurper()
	def counts = []
		json = jsonSlurper.parse(demuxMetrics)
		for (sample in json["samples"]) {
			counts.add(tuple(sample.value["name"], sample.value["reads"][0]))
		}
	return counts
}

// Create a bcl-convert samplessheet for 'fastq_samples' in samples.json
process makeBclConvertSamplesheet {
input: 
	path(samplesCsv)
	path(lib)
	path(runinfo)
	path(fqIndex)
output: 
	path("samplesheet.csv")
publishDir "${params.outDir}/fastq", mode: 'copy'
tag 'small'

"""
	bclConvertSheet.py $samplesCsv $lib $runinfo --fastqIndex $fqIndex > samplesheet.csv
"""
}

// Make TSS Regions BED from gene annoation (GTF) if not already provided as input
process makeTssRegions {
input: path(gtf)
output: path("*.bed")
tag 'small'

"""
	tss_regions.py $gtf
"""
}

/*Run bcl-convert, used when starting from a sequencing run-folder
  Requires a separate bcl-convert samplesheet*/
process bclconvert {
input: 
	path(run)
	path(samplesheet)
output: 
	path("fastq/*fastq.gz"), emit: fastq
	path("fastq/Reports"), emit: stats
publishDir "${params.outDir}/", pattern: 'fastq/Reports/*', mode: 'copy'
publishDir "${params.outDir}/", pattern: 'fastq/*.fastq.gz'

script:
	pthreads = ((task.cpus-4)/3).round()
"""
	bcl-convert --sample-sheet $samplesheet --bcl-num-conversion-threads $pthreads --bcl-num-compression-threads $pthreads --bcl-num-decompression-threads $pthreads --bcl-input-directory $run  --output-directory fastq
"""
}

/*Remove adapter sequences from ends of reads
 Run for fastq input only; otherwise adapters should be removed during bcl-convert*/
process trimFq {
input:
	path(fastq)
output: 
	path("trimmed/${basename}.fastq.gz"), emit: fastq
	path("trimmed/${basename}.trim_stats"), emit: stats

script:
	basename = fastq.getSimpleName()
"""
	mkdir trimmed
	cutadapt -j${task.cpus} -a ${params.adapter} -o trimmed/${basename}.fastq.gz $fastq > trimmed/${basename}.trim_stats
"""
}

// Run bcParser to extract and correct cell-barcodes
// Optionally demuxes fastq files based on some barcodes
process barcodeDemux {
input:
	path(sheet) // Samplesheet json
	path(libDir) // Directory containing the library definition file (barcode sequence lists are loaded from here)
	val(libName) // Filename of the library definition file
	tuple(val(fqName), path(fqFiles)) // Input fastq file
output:
	tuple(val(fqName), path("$outDir/*_S[1-9]*.fastq.gz"), emit: fastq)
	path("$outDir/*_S0_*.fastq.gz"), emit: unknown optional true
	path("$outDir/*.tsv")
	tuple(val(fqName), path("$outDir/metrics.json"), emit: metrics)
publishDir "${params.outDir}/demux/", pattern: "$outDir/*gz"
publishDir "${params.outDir}/demux/", mode: 'copy', pattern: "$outDir/*{txt,tsv,json}" //, saveAs: {"${fqName}.${it.getName()}"}

script:
	outDir = "${fqName}.demux"
	lib = "$libDir/$libName"
"""
	bc_parser --library $lib --demux $sheet --fastq-name $fqName -v --reads ${fqFiles.join(" ")} --write-fastq --out $outDir
"""
}

process align {
input: 
	path(indexDir) // Bowtie2 index directory
	val(indexName) // Base-name of the bowtie2 index
	tuple(val(sample), path(reads)) // Input reads (fastq)
output: 
	tuple(val(sample), path("${sample}.bam"), emit: bam)
	path("*.bam.csi")
	path("*.log")
publishDir "$params.outDir/align", pattern: '*bam*'
publishDir "$params.outDir/align", pattern: '*.log', mode:'copy'

script:
	index = "$indexDir/$indexName"
	athreads = task.cpus - 1
	sthreads = 4
"""
	bowtie2 -p $athreads -x $index -1 ${reads[0]} -2 ${reads[1]} 2> ${sample}.log | samtools view -b | samtools sort --threads ${sthreads} --write-index -o ${sample}.bam
"""
}

process sampleReport {
input: 
	tuple(val(sample), path("metrics.json"))
	path(samplesheet)
	path(libJson)
output:
	path("${sample}.csv")
publishDir "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
"""
	#generateReport.py --sample ${sample} --samplesheet ${samplesheet} --libStruct ${libJson}
	cat $samplesheet > ${sample}.csv
"""
}

//// Sub-Workflows
// Fastq generation, trimming, QC, barcode extraction and sample demux
workflow inputReads {
take:
	samples // samples.csv parsed into a channel
	samplesCsv // samples.csv file
	libJson // library definition json
	runDir // Path to sequencing run-folder (BCLs); empty if running from fastq input
	fqDir // Path to directory with input fastqs; empty if running from BCL 
main:
	runInfo = runDir.map{it.resolve("RunInfo.xml")}

	if (params.fastqSamplesheet == null) {
		fqIndex = expandPath(params.fastqIndex, file("${projectDir}/references/"))
		makeBclConvertSamplesheet(samplesCsv, libJson, runInfo, fqIndex)
		fqSheet = makeBclConvertSamplesheet.out
	} else {
		fqSheet = file(params.fastqSamplesheet)
	}
	bclconvert(runDir, fqSheet)
	fqs = fqDir.flatMap{file(it).listFiles()}.filter{it.name =~ /.fastq.gz$/}
	fqs.dump(tag:'fqs')
	fqs = bclconvert.out.fastq.flatten().mix(fqs)
	// Organize fastq files by sample
	fqs.dump(tag:'fqs2')
	fqFiles = fqs.map { file ->
		def ns = file.getName().toString().tokenize('_')
		return tuple(ns.get(0), file)
	}.groupTuple()
	fqSamples = samples.map({it.fastqName}).unique().join(fqFiles)
	fqSamples.dump(tag:'fqSamples')

	// Process cell-barcodes and (optionally) split fastqs into samples based on tagmentation barcode
	barcodeDemux(samplesCsv, libJson.getParent(), libJson.getName(), fqSamples)
	demuxFqs = barcodeDemux.out.fastq.flatMap({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('_')
		return tuple(ns.get(0), file)
	}.groupTuple(size:2)
	demuxFqs.dump(tag:'demuxFqs')
emit:
	fqs = demuxFqs
	metrics = barcodeDemux.out.metrics
}

// Alignments, fragments, peaks and quantification
workflow scATAC {
take:
	samples
	sampleFqs // Fastq files for each sample for all reads (including index reads)
	demuxMetrics
main:
	align(genome.bowtie_index.getParent(), genome.bowtie_index.getName(), sampleFqs)
	bams = align.out.bam
	sampleDemuxMetrics = demuxMetrics.cross(samples.map{[it.fastqName,it.sample]}).map({[it[1][1], it[0][0], it[0][1]]})
emit:
	sampleDemuxMetrics = sampleDemuxMetrics
}

// QC, reporting and preliminary analysis
workflow atacReport {
take:
	samples
	demuxMetrics
	samplesheet
	libJson
main:
	demux = demuxMetrics.map({[it[0], it[2]]})
	sampleReport(demux, samplesheet, libJson)
}


//// Main entry point
// Run the workflow for one or multiple samples
// either from one runFolder or one / multiple sets of fastq files
workflow {
	Helper.initialise(workflow, params, log)
	//// Inputs
	samplesCsv = file(params.samples)
	samples = Channel.fromPath(samplesCsv).splitCsv(header:true, strip:true)
	samples.dump(tag:'samples')
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	// Input reads from runfolder xor fastq directory
	// runDir or fqDir can either be set as parameter or as a column in samples.csv
	if (params.runFolder) {
		runDir = Channel.fromPath(params.runFolder, checkIfExists:true)
		fqDir = Channel.empty()
	} else if (params.fastqDir) {
		fqDir = Channel.fromPath(params.fastqDir, checkIfExists:true)
		runDir = Channel.empty()
	} else {
		//todo should ensure only one of the two is set
		runDir = samples.map{it.runFolder}.filter{it}.first().map{file(it,checkIfExists:true)}
		fqDir = samples.map{it.runFolder}.filter{it}.first().map{file(it,checkIfExists:true)}
	}
	inputReads(samples, samplesCsv, libJson, runDir, fqDir)
	genome = loadGenome(file(params.genome))
	scATAC(samples, inputReads.out.fqs, inputReads.out.metrics)
	atacReport(samples, scATAC.out.sampleDemuxMetrics, samplesCsv, libJson)
}

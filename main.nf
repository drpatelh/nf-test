nextflow.enable.dsl=2
import groovy.json.JsonSlurper

// Load .json file into object
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

def throwError(errMessage) {
	log.error errMessage
	sleep(200)
	System.exit(1)
}

/*
Ensures that FASTQS for each lane of each sample in the @fqFilesChannel
have four labelled R1, R2, I1 and I2 
*/
def validateFastqs(fqFilesChannel) {
	def fastqRegexes = [/_R1/, /_R2/, /_I1/, /_I2/]
	filesGroupedBySampleAndName = fqFilesChannel.map { sampleName, sampleFileList ->
		readNumberRegex = /_[RI][1-2]/
		[sampleName, sampleFileList.groupBy({x ->
			def y = x.getName().toString().tokenize('_')
			def name = x.getName().toString()
			def match = (name =~ readNumberRegex)
			if (match.find()) {
				return name.replace(match.group(0), "<readNumber>")
			} else {
				throwError("FASTQ: '${name}' doesn't contain a readNumber specifying it as a read or index:\n Must contain one of: ${fastqRegexes}")
			}
		})]
	}
	filesGroupedBySampleAndName.map({ sampleName, fileNameHashMap -> 
		def keys = fileNameHashMap.keySet()
		for (name in keys) {
			fastqs = fileNameHashMap.get(name)
			for (pattern in fastqRegexes) {
				requiredFile = fastqs.findAll { it.getBaseName() =~ pattern }
				readType = pattern.replace("_","")
				if ( requiredFile.size == 0 ) {
					throwError("Sample '${sampleName}' is missing FASTQ with readNumber ${readType} for FASTQs following pattern:\n ${name}")
				} else if (requiredFile.size > 1) {
					throwError("Sample '${sampleName}' has more than one FASTQ specified of type ${readType} for group following pattern ${name}")
				}
			}
		}
	})
}

/*
Ensures each row in params.samples(@sampleCsvRows) contains valid 
values for both required and optional columns
*/
def validateSamples(sampleCsvRows) {
	def requiredColumns = ['sample', 'fastqName']
	def optionalNaturalNumberColumns = ['expectedCells', 'subsample']
	sampleCsvRows.map{ csvRow ->
		sampleName = csvRow.get('sample', "unknown")
		for (column in requiredColumns) {
			def val = csvRow.get(column)
			def regexPattern = /[0-9a-zA-Z][0-9a-zA-Z\-\.]*/
			if (val == null) {
				throwError("Required column '${column}' missing from specified samples table: ${params.samples}")
			} else if (!val.matches(regexPattern)) {
				throwError("Invalid value specified for sample '${sampleName}' in column '${column}': `${val}`\nIn samples csv: ${params.samples}\nValues in the `${column}` column must:\n- Begin with an alphanumeric character\n- Consist of only alphanumeric and the following special characters: '-' '.'")
			}
		}
		for (column in optionalNaturalNumberColumns) {
			def val = csvRow.get(column)
			// Allow 0 or '' to indicate 'default'
			def regexPattern = /[0-9]*/
			if (val != null && !val.matches(regexPattern)) {
				throwError("Invalid value specified for sample `${sampleName}` in column `${column}`: `${val}`\n- Must be an integer >= 0")
			}
		}
	}
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

// Run fastQC on input files inputs or from bcl-convert outputs)
// Note that these are the input fastq files, not outputs of bcParser
process fastqc {
input:
	path(fq)
output:
	path("fastqc/*.html"), emit: html
	path("fastqc/*.zip"), emit: zip
tag 'small'

"""
	mkdir fastqc
	fastqc -o fastqc $fq
"""
}

// Use multiQC to report a single QC report for all fastQC and bcl-convert reports
process multiqc {
input:
	path(reports)
output:
	path("multiqc_report.html")
publishDir "${params.outDir}/fastq", mode: 'copy'
tag 'small'

"""
	multiqc .
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

// bowtie2 alignment to the genome
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

// Subsample per-sample BAM files to a fixed depth
process subsampleBam {
input: 
	tuple(val(sample), path(bam), val(sampleReadCount), val(targetReadCount))
output: 
	tuple(val(sample), path(outFn))
publishDir "$params.outDir/align/subsample", pattern: '*bam*'

script:
	fraction = targetReadCount / sampleReadCount
	if ((fraction) > 0 && (fraction < 1)) {
		outFn = "${sample}.subsampled.bam"
"""
	samtools view -b -s $fraction -o $outFn $bam
"""
	} else { // Pass entire input bam through
	outFn = bam
"""
	echo "Taking all reads for $bam"
"""
	}
}

// Run scDedup to remove duplicate reads based on cell-barcode and mapping position
// Also generates cell and fragment statistics
process dedup {
input:
	tuple(val(sample), file(bam))
output: 
	tuple(val(sample), path("${sample}.dedup.bam"), emit: bam)
	tuple(val(sample), path("${sample}.fragments.tsv.gz"), path("${sample}.fragments.tsv.gz.tbi"), emit: fragments)
	tuple(val(sample), path("${sample}.cell_stats.tsv"), path("${sample}.dedup_stats.tsv"), path("${sample}.fragment_hist.tsv"), emit: stats)
publishDir "$params.outDir/align/dedup/", pattern: '*bam*'
publishDir "$params.outDir/align/", pattern: '*.fragments.tsv.*'
publishDir "$params.outDir/align/dedup/", pattern: '*.tsv', mode:'copy'

"""
	sc_dedup $bam --barcode-input Qname --write-fragments --out-prefix $sample
	samtools index ${sample}.dedup.bam &
	tabix -p bed ${sample}.fragments.tsv.gz &
	wait
"""
}

// Convert BAM to BED with one entry for each tagmentation event (read end)
// Includes unpaired reads
process tagSites {
input: 
	tuple(val(sample), file(bam))
output: 
	tuple(val(sample), path("${sample}.bed.gz"))
publishDir "$params.outDir/peaks"

"""
	bedtools bamtobed -i $bam -tag XC | gzip -c >${sample}.bed.gz
"""
}

// Per-sample de-novo peak calling from tagmentation events
process callPeaks {
input: 
	tuple(val(sample), file(tagSites))
output: 
	tuple(val(sample), path("*_peaks.narrowPeak"))
publishDir "$params.outDir/peaks", mode: 'copy'

script:
"""
	macs3 callpeak -t ${tagSites} -f BED --nomodel --shift -100 --extsize 200 --keep-dup all -n $sample
"""
}

// Cell X peak count matrix
// Either using sample-specific peak calls or pre-defined peak set
process countPeaks {
input: 
	tuple(val(sample), path(tagSites), path(peaks))
output: 
	tuple(val(sample), path("${sample}.counts"))
publishDir "$params.outDir/peaks", mode: 'copy'

script:
"""
	bedtools intersect -a <(gunzip -c $tagSites) -b $peaks -loj | sc_counter --out-dir ${sample}.counts
"""
}

// Cell X gene count matrix
// Based on window around all transcript 5' ends
process countTss {
input:
	path(tss)
	tuple(val(sample), path(tagSites))
output: 
	tuple(val(sample), path("${sample}.tss_counts"))
publishDir "$params.outDir/peaks", mode: 'copy'

"""
	bedtools intersect -a <(gunzip -c $tagSites) -b $tss -loj | sc_counter --out-dir ${sample}.tss_counts
"""
}

/* 
Generates: 
- table for each sample denoting which cells have passed set qcFilters (qc.tsv)
- calculated thresholds for set filters (qc.json)
- figures displaying the number of cells filtered at each step (Figs)
*/
process cellFilter {
	// dedup_stats & fragment_hist.tsv are not currently in use by atacQcFilter but are named as they are 
	// passed in the same channel cell_stats.tsv is received from (dedup.out.stats) may be used in future
input:
	tuple(val(sample), val(expectedCells), path("${sample}.cell_stats.tsv"), path("dedup_stats.tsv"), path("fragment_hist.tsv"), path("*"))
output:
	tuple(val(sample), path("QC/${sample}/${sample}_thresholds.json"), path("QC/${sample}/${sample}_QC.tsv"), emit: qcStats) 
	tuple(val(sample), path("QC/${sample}/fripKneePlot.png"), path("QC/${sample}/uniqueReadsKneePlot.png"), emit: qcFigs)
publishDir "$params.outDir", mode: 'copy'

script:
	qcAndArchRParams = Helper.getqcAndArchRParams(params.qcAndArchRParams)
	expectedCells ?= 0
"""
	atacQcFilter.py --sample $sample --minUniqueReads $qcAndArchRParams.minUniqueReads --expectedCells $expectedCells --background $qcAndArchRParams.backgroundRatio --topPercentCells $qcAndArchRParams.topPercentCells
""" 
} 

/*
Runs automated ArchR analysis for each sample 
- path("*") captures thresholds.json which is emitted along with needed QC.tsv by the cellFilter process. 
  It is not needed for this process  
*/
process archrAnalysis { 
input:
	tuple(val(sample), path("${sample}.fragments.tsv.gz"), path("${sample}.fragments.tsv.gz.tbi"), path("*"), path("QC.tsv"))
output:
	path("ArchR/${sample}")
publishDir "$params.outDir" 
errorStrategy 'ignore'
	if (workflow.profile == 'conda'){
		if (params.archrCondaEnv != null) {
			conda params.archrCondaEnv
		} else {
			throwError("In order to run archrAnalysis process with conda profile\na custom conda environment needs to be created and specified using the\n'archrCondaEnv' field in your params file. See the Dependencies section of README for more info")
		} 
	}
script:
	additionalArgs=""
	// Use gtf and specified BSGenome 
	if (genome.BSGenome != null){
		additionalArgs="$additionalArgs -g ${genome.gtf} -b ${genome.BSGenome}"
	}
	// Use Built in ArchRGenome && hasBuiltInGenome
	else if (genome.archrAlias != null) {
		additionalArgs="$additionalArgs -a $genome.archrAlias"
	} 
"""
	ArchR_analysis.R -s '${genome.speciesName}' -f ${sample}.fragments.tsv.gz -sN $sample -o $sample -r ${genome.annotationConvention} -th ${task.cpus} $additionalArgs -qc QC.tsv -q $params.qcAndArchRParams
"""
}

process sampleReport {
input: 
	tuple(val(sample), path("metrics.json"), path("${sample}.cell_stats.tsv"), path("${sample}.dedup_stats.tsv"), path("${sample}.fragment_hist.tsv"), path("*"), path("${sample}_QC.json"), path("${sample}_QC.tsv"))
	path(samplesheet)
	path(libJson)
output:
	path("reports/${sample}.report.html")
	path("reports/${sample}.reportStatistics.tsv")
publishDir "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
"""
	generateReport.py --sample ${sample} --samplesheet ${samplesheet} --libStruct ${libJson}
"""
}

process fastqReport {
input: 
	tuple(val(fastqName), path(files))
	path(samplesheet)
	path(libJson)
output:
	path("reports/${fastqName}.report.html")
publishDir "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
"""
	generateReport.py --fastqName ${fastqName} --samplesheet ${samplesheet} --libStruct ${libJson}
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
	if (params.trimFastq) {
		readFqs = fqs.filter { it.getBaseName() =~ /_R\d_/ }
		indexFqs = fqs.filter { it.getBaseName() =~ /_I\d_/ }
		trimFq(readFqs)
		fqs = trimFq.out.fastq.mix(indexFqs)
	}
	fqs = bclconvert.out.fastq.flatten().mix(fqs)
	// Organize fastq files by sample
	fqs.dump(tag:'fqs2')
	fqFiles = fqs.map { file ->
		def ns = file.getName().toString().tokenize('_')
		return tuple(ns.get(0), file)
	}.groupTuple()
	fqSamples = samples.map({it.fastqName}).unique().join(fqFiles)
	validateFastqs(fqSamples)
	fqSamples.dump(tag:'fqSamples')

	// Process cell-barcodes and (optionally) split fastqs into samples based on tagmentation barcode
	barcodeDemux(samplesCsv, libJson.getParent(), libJson.getName(), fqSamples)
	demuxFqs = barcodeDemux.out.fastq.flatMap({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('_')
		return tuple(ns.get(0), file)
	}.groupTuple(size:2)
	demuxFqs.dump(tag:'demuxFqs')
	if (params.fastqc) {
		fastqc(demuxFqs.flatMap({it.get(1)}))
		reports = fastqc.out.zip
		if (runDir != null) {
			reports = reports.mix(bclconvert.out.stats)
		}
		if (params.trimFastq) {
			reports = reports.mix(trimFq.out.stats)
		}
		multiqc(reports.collect())
	}
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
	tssBed
main:
	align(genome.bowtie_index.getParent(), genome.bowtie_index.getName(), sampleFqs)
	if (params.subsample) {
		readCounts = demuxMetrics
		.flatMap{loadDemuxReadCounts(it[1])}
		.join(samples.map{[it.sample, toIntOr0(it.subsample)]})
		readCounts.dump(tag:'readCounts')
		subsampleBam(align.out.bam.combine( readCounts, by:0))
		bams = subsampleBam.out
	} else {
		bams = align.out.bam
	}
	dedup(bams)
	tagSites(dedup.out.bam)
	if (params.peaks == null) {
		callPeaks(tagSites.out)
		peaks = callPeaks.out
	} else {
		peaks = samples.map{[it.sample, file(params.peaks)]}
	}
	peakInput = tagSites.out.join(peaks)
	peakInput.dump(tag:'peakInput')
	countPeaks(peakInput)
	if (tssBed != null) { 
		tssBed.dump(tag:'tssBed')
		countTss(tssBed, tagSites.out)
		tssCounts = countTss.out
	} else {
		tssCounts = Channel.empty()
	}
	counts = countPeaks.out.join(tssCounts, remainder: true).map { [it[0], it[1..-1].findAll {it}] }
	// Using cross instead of join to multimatch samples to their fastqName
	// Only getting sampleName, fastqName and demuxPath using map 
	sampleDemuxMetrics = demuxMetrics.cross(samples.map{[it.fastqName,it.sample]}).map({[it[1][1], it[0][0], it[0][1]]})
	sampleDemuxMetrics.dump(tag:'sampleDemuxMetrics')
emit:
	sampleDemuxMetrics = sampleDemuxMetrics
	sampleDedupMetrics = dedup.out.stats
	counts = counts
	fragments = dedup.out.fragments
}

// QC, reporting and preliminary analysis
workflow atacReport {
take:
	samples
	demuxMetrics
	dedupMetrics
	counts
	fragments
	samplesheet
	libJson
main:
	expectedCells = samples.map{[it.sample, toIntOr0(it.expectedCells)]}

	dedupStatsAndCounts = dedupMetrics.join(counts)
	cellFilter(expectedCells.join(dedupStatsAndCounts))
	canUseBuiltInGenome = (genome.archrAlias != null)
	canUseCustomGenome = (genome.BSGenome != null)
	if ( params.runArchR && (canUseBuiltInGenome || canUseCustomGenome)) {
		archrAnalysis(fragments.join(cellFilter.out.qcStats))
	}
	if (params.generateReport) {
		// Removing fastqName from demuxMetrics as not needed for cellFilter or sampleReport 
		demuxAndDedup = demuxMetrics.map({[it[0], it[2]]}).join(dedupMetrics)
		allStats = demuxAndDedup.join(counts)
		allStats.dump(tag:'allStats')
		allStatsAndQC = allStats.join(cellFilter.out.qcStats)
		sampleReport(allStatsAndQC, samplesheet, libJson)

		// Creating a seperate channel from allStats since fastqReport only needs QC and demux metrics 
		multiReportData = demuxMetrics.join(cellFilter.out.qcStats)
		// by:[1,2] where 1 and 2 are fastqName and demux metrics.json path respectively
		// .map grouped into [fastqName, [all needed metrics files for all samples originating from fastqName]]
		multiReportDataFormatted = multiReportData.groupTuple(by:[1,2]).map({[it[1], it[2..-1].flatten()]})
		fastqReport(multiReportDataFormatted, samplesheet, libJson)
	}
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
	validateSamples(samples)
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
	// Reference Genome
	genome = loadGenome(file(params.genome))
	if ((genome.tss == null) && (genome.gtf != null)) {
		makeTssRegions(genome.gtf)
		tssBed = makeTssRegions.out
	} else {
		tssBed = genome.tss
	}

	//// scATAC Analysis and reporting
	scATAC(samples, inputReads.out.fqs, inputReads.out.metrics, tssBed)
	atacReport(samples, scATAC.out.sampleDemuxMetrics,  scATAC.out.sampleDedupMetrics, scATAC.out.counts, scATAC.out.fragments, samplesCsv, libJson)
}

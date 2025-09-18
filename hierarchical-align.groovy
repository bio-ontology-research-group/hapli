#!/usr/bin/env groovy

import groovy.transform.Field

@Field def features = [:]
@Field def sequences = [:]
@Field def children = [:].withDefault { [] }
@Field def roots = []
@Field def threads
@Field def pafFile = null
@Field def refSeq
@Field def getFeatureSequence

def cli = new CliBuilder(usage: 'hierarchical_align.groovy [options]')
cli.with {
    g(longOpt: 'gff', args: 1, argName: 'file', required: true, 'GFF3 annotation file')
    r(longOpt: 'reference-fasta', args: 1, argName: 'file', required: true, 'Reference genome FASTA file')
    p(longOpt: 'paths', args: 1, argName: 'file', required: true, 'Paths FASTA file')
    o(longOpt: 'output', args: 1, argName: 'file', required: false, 'Output file (default: hierarchical_alignments.txt)')
    s(longOpt: 'start', args: 1, argName: 'int', required: true, 'Region start coordinate')
    e(longOpt: 'end', args: 1, argName: 'int', required: true, 'Region end coordinate')
    c(longOpt: 'chr', args: 1, argName: 'string', required: false, 'Chromosome (default: chr1)')
    t(longOpt: 'threads', args: 1, argName: 'int', required: false, 'Number of threads (default: 8)')
    gn(longOpt: 'genes', args: 1, argName: 'list', required: false, 'Comma-separated list of gene names to align')
    tp(longOpt: 'target-paths', args: 1, argName: 'list', required: false, 'Comma-separated list of paths to align against')
    a(longOpt: 'paf-output', args: 1, argName: 'file', required: false, 'Save raw minimap2 PAF alignments to a file')
    h(longOpt: 'help', 'Show usage')
}

def options = cli.parse(args)
if (!options || options.h) {
    cli.usage()
    return
}

def gffFile = new File(options.g)
def referenceFastaFile = new File(options.r)
def pathsFile = new File(options.p)
def outputFile = new File(options.o ?: 'hierarchical_alignments.txt')
def regionStart = options.s as Integer
def regionEnd = options.e as Integer
def chromosome = options.c ?: 'chr1'
threads = (options.t ?: '8') as Integer
def targetGenes = options.gn ? options.gn.split(',').collect { it.trim() } : []
def targetPaths = options.tp ? options.tp.split(',').collect { it.trim() } as Set : null
pafFile = options.a ? new File(options.a) : null

if (pafFile) {
    pafFile.text = '' // Clear PAF file at start
    println "Will save raw PAF alignments to ${pafFile.name}"
}

if (targetPaths) {
    println "Filtering for ${targetPaths.size()} target paths."
    def filteredPathsFile = File.createTempFile("filtered_paths", ".fa")
    filteredPathsFile.deleteOnExit()
    def writer = filteredPathsFile.newWriter()
    def currentPathId = null
    def writing = false
    pathsFile.eachLine { line ->
        if (line.startsWith('>')) {
            currentPathId = line.substring(1).split()[0]
            writing = targetPaths.contains(currentPathId)
        }
        if (writing) {
            writer.println line
        }
    }
    writer.close()
    pathsFile = filteredPathsFile // Use filtered paths for all subsequent operations
}

// Parse GFF3 hierarchy
features = [:]
children = [:].withDefault { [] }
roots = []

println "Parsing GFF3 hierarchy..."
gffFile.eachLine { line ->
    if (line.startsWith('#')) return
    def parts = line.split('\t')
    if (parts.size() < 9) return
    
    def chr = parts[0]
    def type = parts[2]
    def start = parts[3] as Integer
    def end = parts[4] as Integer
    def strand = parts[6]
    def attributes = parts[8]
    
    if (chr != chromosome || end < regionStart || start > regionEnd) return
    
    def id = (attributes =~ /ID=([^;]+)/).findResult { it[1] }
    def parent = (attributes =~ /Parent=([^;]+)/).findResult { it[1] }
    def geneName = (attributes =~ /gene_name=([^;]+)/).findResult { it[1] } ?: ""
    
    if (id) {
        def adjStart = Math.max(0, start - regionStart)
        def adjEnd = Math.min(regionEnd - regionStart, end - regionStart)
        
        features[id] = [
            id: id,
            type: type,
            start: adjStart,
            end: adjEnd,
            origStart: start,
            origEnd: end,
            strand: strand,
            geneName: geneName,
            parent: parent,
            children: []
        ]
        
        if (parent) {
            children[parent] << id
        } else {
            roots << id
        }
    }
}

children.each { parentId, childIds ->
    if (features[parentId]) {
        features[parentId].children = childIds
    }
}

println "Found ${roots.size()} root features and ${features.size()} total features"

// Load reference sequence for the target chromosome
println "Loading reference sequence for chromosome ${chromosome}..."
def refSeqBuilder = new StringBuilder()
def currentChrId = null
def readingChr = false
referenceFastaFile.eachLine { line ->
    if (line.startsWith('>')) {
        if (readingChr) {
            readingChr = false // Stop reading on next chromosome
        }
        currentChrId = line.substring(1).split()[0]
        if (currentChrId == chromosome) {
            readingChr = true
        }
    } else if (readingChr) {
        refSeqBuilder.append(line.trim())
    }
}
refSeq = refSeqBuilder.toString()

if (refSeq.length() == 0) {
    System.err.println "Error: Chromosome '${chromosome}' not found or has no sequence in reference FASTA file ${referenceFastaFile.name}"
    return
}
println "Loaded reference sequence for ${chromosome} (${refSeq.length()} bp)"

// Sequence extraction function
getFeatureSequence = { feature, referenceSequence ->
    if (!feature) return null
    // GFF is 1-based, substring is 0-based exclusive end
    def start = feature.origStart - 1
    def end = feature.origEnd
    
    if (start < 0 || end > referenceSequence.length() || start >= end) {
        System.err.println "Warning: Feature ${feature.id} coordinates (${feature.origStart}-${feature.origEnd}) are invalid or out of bounds for reference."
        return null
    }
    
    def seq = referenceSequence.substring(start, end)
    
    if (feature.strand == '-') {
        seq = seq.reverse().tr('ACGTNacgtn', 'TGCANtgcan')
    }
    return seq
}

// Alignment function
def alignFeature(featureId, pathFile, numThreads, parentRegion = null) {
    def feature = features[featureId]
    if (!feature) return null
    
    def seq = getFeatureSequence(feature, refSeq)
    if (!seq) {
        System.err.println "Warning: Could not extract sequence for feature ${featureId}, skipping alignment."
        return []
    }
    def seqLen = seq.length()
    def results = []
    
    def queryFile = File.createTempFile("query", ".fa")
    queryFile.text = ">${featureId}\n${seq}\n"
    
    def cmd
    if (seqLen < 30) {
        cmd = "minimap2 -c -cx sr -k5 -w2 --min-dp-score 10 --min-chain-score 10 -B2 -O2,12 -E1,0 --score-N 0 --secondary=yes -N 500 -t ${numThreads}"
    } else if (seqLen < 100) {
        cmd = "minimap2 -c -cx sr -k7 -w3 --min-dp-score 20 --min-chain-score 20 -B3 -O3,18 --secondary=yes -N 500 -t ${numThreads}"
    } else if (seqLen < 500) {
        cmd = "minimap2 -c -cx asm20 -k10 -w5 --min-dp-score 50 --secondary=yes -N 500 -t ${numThreads}"
    } else {
        cmd = "minimap2 -c -cx asm5 --secondary=yes -N 500 -t ${numThreads}"
    }
    
    def targetFile = pathFile
    if (parentRegion && parentRegion.tempFile) {
        targetFile = parentRegion.tempFile
    }
    
    cmd += " ${targetFile} ${queryFile}"
    
    def proc = cmd.execute()
    def output = proc.in.text 
    def errors = proc.err.text
    
    proc.waitFor()

    if (pafFile) {
        pafFile.append(output)
    }
    
    if (proc.exitValue() != 0) {
        System.err.println "Error running minimap2 for feature ${featureId}. Exit code: ${proc.exitValue()}"
        System.err.println "Stderr: ${errors}"
        return []
    }

    output.eachLine { line ->
        def parts = line.split('\t')
        if (parts.size() >= 11) {
            def targetPath = parts[5]
            def targetStart = parts[7] as Integer
            def targetEnd = parts[8] as Integer
            def matches = parts[9] as Integer
            def alignLen = parts[10] as Integer
            def identity = alignLen > 0 ? (matches / alignLen * 100) : 0
            
            def inRegion = true
            if (parentRegion) {
                def tolerance = Math.max(50, seqLen)
                inRegion = (targetStart >= parentRegion.start - tolerance && 
                           targetEnd <= parentRegion.end + tolerance)
            }
            
            if (inRegion) {
                results << [
                    feature: featureId,
                    type: feature.type,
                    path: targetPath,
                    start: targetStart,
                    end: targetEnd,
                    identity: identity,
                    length: seqLen
                ]
            }
        }
    }
    
    queryFile.delete()
    return results
}

println "\nStarting hierarchical alignment..."
outputFile.withWriter { writer ->
    writer.println "HIERARCHICAL FEATURE ALIGNMENT RESULTS"
    writer.println "Region: ${chromosome}:${regionStart}-${regionEnd}"
    writer.println "=" * 80
    
    def allRoots = roots.sort()
    def filteredRoots = allRoots
    if (targetGenes.size() > 0) {
        filteredRoots = allRoots.findAll { geneId ->
            def gene = features[geneId]
            gene && targetGenes.contains(gene.geneName)
        }
        println "\nFiltered to ${filteredRoots.size()} of ${allRoots.size()} total root features based on specified gene names."
    }

    def totalGenes = filteredRoots.size()
    filteredRoots.eachWithIndex { geneId, i ->
        def gene = features[geneId]
        
        // **NEW**: Print progress for each gene
        println "[${i + 1}/${totalGenes}] Aligning gene: ${gene.geneName ?: geneId}"
        System.out.flush()

        writer.println "\nGene: ${gene.geneName ?: geneId} (${gene.type})"
        writer.println "  Original coordinates: ${gene.origStart}-${gene.origEnd}"
        writer.println "  Region coordinates: ${gene.start}-${gene.end}"
        
        def geneAlignments = alignFeature(geneId, pathsFile, threads)
        
        if (!geneAlignments) {
            writer.println "  No alignments found"
            return
        }
        
        def uniquePaths = geneAlignments.collect { it.path }.unique()
        writer.println "  Aligned to ${uniquePaths.size()} paths"
        
        geneAlignments.groupBy { it.path }.each { path, alns ->
            def aln = alns[0]
            writer.println "\n  Path: ${path}"
            writer.println "    Gene: ${aln.start}-${aln.end} (${String.format('%.1f%%', aln.identity)}, ${aln.length}bp)"
            
            gene.children.each { transcriptId ->
                def transcript = features[transcriptId]
                if (!transcript) return

                println "  -> Aligning transcript: ${transcriptId}"
                System.out.flush()
		
                def region = [start: aln.start, end: aln.end]
                def txAlignments = alignFeature(transcriptId, pathsFile, threads, region)
                
                if (txAlignments) {
                    def txAln = txAlignments.find { it.path == path }
                    if (txAln) {
                        writer.println "    ${transcript.type}: ${transcriptId}"
                        writer.println "      ${txAln.start}-${txAln.end} (${String.format('%.1f%%', txAln.identity)}, ${txAln.length}bp)"
                        
                        transcript.children.each { childId ->
                            def child = features[childId]
                            if (!child) return

                            println "    -> Aligning sub-feature: ${childId} (${child.type})"
                            System.out.flush()
			    
                            def childRegion = [start: txAln.start, end: txAln.end]
                            def childAlns = alignFeature(childId, pathsFile, threads, childRegion)
                            if (childAlns) {
                                def childAln = childAlns.find { it.path == path }
                                if (childAln) {
                                    writer.println "        ${child.type}: ${childAln.start}-${childAln.end} (${String.format('%.1f%%', childAln.identity)}, ${childAln.length}bp)"
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    writer.println "\n" + "=" * 80
    writer.println "Analysis complete"
}

println "\nHierarchical alignment complete. Results in ${outputFile.name}"

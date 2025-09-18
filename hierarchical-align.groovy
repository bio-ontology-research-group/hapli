#!/usr/bin/env groovy

import groovy.transform.Field

@Field def features = [:]
@Field def sequences = [:]
@Field def children = [:].withDefault { [] }
@Field def roots = []
@Field def threads

def cli = new CliBuilder(usage: 'hierarchical_align.groovy [options]')
cli.with {
    g(longOpt: 'gff', args: 1, argName: 'file', required: true, 'GFF3 annotation file')
    f(longOpt: 'features', args: 1, argName: 'file', required: true, 'Features FASTA file')
    p(longOpt: 'paths', args: 1, argName: 'file', required: true, 'Paths FASTA file')
    o(longOpt: 'output', args: 1, argName: 'file', required: false, 'Output file (default: hierarchical_alignments.txt)')
    s(longOpt: 'start', args: 1, argName: 'int', required: true, 'Region start coordinate')
    e(longOpt: 'end', args: 1, argName: 'int', required: true, 'Region end coordinate')
    c(longOpt: 'chr', args: 1, argName: 'string', required: false, 'Chromosome (default: chr1)')
    t(longOpt: 'threads', args: 1, argName: 'int', required: false, 'Number of threads (default: 8)')
    h(longOpt: 'help', 'Show usage')
}

def options = cli.parse(args)
if (!options || options.h) {
    cli.usage()
    return
}

def gffFile = new File(options.g)
def fastaFile = new File(options.f)
def pathsFile = new File(options.p)
def outputFile = new File(options.o ?: 'hierarchical_alignments.txt')
def regionStart = options.s as Integer
def regionEnd = options.e as Integer
def chromosome = options.c ?: 'chr1'
threads = (options.t ?: '8') as Integer

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

// Load feature sequences
sequences = [:]
def currentId = null
def currentSeq = new StringBuilder()

fastaFile.eachLine { line ->
    if (line.startsWith('>')) {
        if (currentId) {
            sequences[currentId] = currentSeq.toString()
        }
        currentId = line.substring(1)
        currentSeq = new StringBuilder()
    } else {
        currentSeq.append(line)
    }
}
if (currentId) {
    sequences[currentId] = currentSeq.toString()
}

println "Loaded ${sequences.size()} feature sequences"

// Alignment function
def alignFeature(featureId, pathFile, numThreads, parentRegion = null) {
    def feature = features[featureId]
    if (!feature) return null
    
    def seqKey = sequences.keySet().find { it.contains(featureId) }
    if (!seqKey) return null
    
    def seq = sequences[seqKey]
    def seqLen = seq.length()
    def results = []
    
    def queryFile = File.createTempFile("query", ".fa")
    queryFile.text = ">${featureId}\n${seq}\n"
    
    def cmd
    if (seqLen < 30) {
        cmd = "minimap2 -cx sr -k5 -w2 --min-dp-score 10 --min-chain-score 10 -B2 -O2,12 -E1,0 --score-N 0 --secondary=yes -N 500 -t ${numThreads}"
    } else if (seqLen < 100) {
        cmd = "minimap2 -cx sr -k7 -w3 --min-dp-score 20 --min-chain-score 20 -B3 -O3,18 --secondary=yes -N 500 -t ${numThreads}"
    } else if (seqLen < 500) {
        cmd = "minimap2 -cx asm20 -k10 -w5 --min-dp-score 50 --secondary=yes -N 500 -t ${numThreads}"
    } else {
        cmd = "minimap2 -cx asm5 --secondary=yes -N 500 -t ${numThreads}"
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
    
    def totalGenes = roots.size()
    roots.sort().eachWithIndex { geneId, i ->
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
                        writer.println "    ${transcript.type}: ${transcriptId.take(30)}"
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

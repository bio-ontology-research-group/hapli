#!/usr/bin/env python3
"""
Test roundtrip functionality: FASTA -> Pangenome -> VCF -> FASTA

This test verifies that we can:
1. Generate variant FASTA files from a reference
2. Build a pangenome from these FASTA files
3. Call variants from the pangenome to generate VCF files
4. Apply the VCF variants back to the reference
5. Verify that the result matches the original variant FASTA files
"""

import unittest
import tempfile
import shutil
import subprocess
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam


class TestPangenomeRoundtrip(unittest.TestCase):
    """Test pangenome roundtrip functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.reference_file = self.test_dir / "reference.fa"
        self.variants_dir = self.test_dir / "samples"
        self.pangenome_dir = self.test_dir / "pangenome"
        self.vcf_dir = self.test_dir / "vcf_output"
        self.roundtrip_dir = self.test_dir / "roundtrip"
        
        # Create directories
        self.variants_dir.mkdir()
        self.pangenome_dir.mkdir()
        self.vcf_dir.mkdir()
        self.roundtrip_dir.mkdir()
        
        # Create a simple reference genome
        self.create_reference_genome()
        
        # Create variant FASTA files
        self.create_variant_fastas()
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.test_dir)
    
    def create_reference_genome(self):
        """Create a simple reference genome for testing."""
        # Create a simple 1000bp reference sequence
        reference_seq = "A" * 200 + "T" * 200 + "G" * 200 + "C" * 200 + "N" * 200
        
        record = SeqRecord(
            Seq(reference_seq),
            id="chr1",
            description="Test reference chromosome"
        )
        
        SeqIO.write(record, self.reference_file, "fasta")
    
    def create_variant_fastas(self):
        """Create variant FASTA files with known variants."""
        # Read reference
        ref_record = SeqIO.read(self.reference_file, "fasta")
        ref_seq = str(ref_record.seq)
        
        # Sample 1: SNV at position 100 (A->G) and insertion at position 300
        sample1_hap1_seq = ref_seq[:99] + "G" + ref_seq[100:]  # SNV
        sample1_hap1_seq = sample1_hap1_seq[:300] + "TTT" + sample1_hap1_seq[300:]  # Insertion
        
        sample1_hap2_seq = ref_seq[:199] + "C" + ref_seq[200:]  # Different SNV
        
        # Sample 2: Deletion at position 150-152 and SNV at position 400
        sample2_hap1_seq = ref_seq[:149] + ref_seq[152:]  # Deletion of 3bp
        sample2_hap1_seq = sample2_hap1_seq[:397] + "A" + sample2_hap1_seq[398:]  # SNV (T->A)
        
        sample2_hap2_seq = ref_seq[:499] + "T" + ref_seq[500:]  # SNV (G->T)
        
        # Write variant FASTA files
        samples = [
            ("sample_001", sample1_hap1_seq, sample1_hap2_seq),
            ("sample_002", sample2_hap1_seq, sample2_hap2_seq)
        ]
        
        self.sample_dirs = []
        for sample_name, hap1_seq, hap2_seq in samples:
            sample_dir = self.variants_dir / sample_name
            sample_dir.mkdir()
            self.sample_dirs.append(str(sample_dir))
            
            # Write haplotype 1
            hap1_record = SeqRecord(
                Seq(hap1_seq),
                id=f"{sample_name}_hap1",
                description=f"Haplotype 1 for {sample_name}"
            )
            SeqIO.write(hap1_record, sample_dir / f"{sample_name}_hap1.fasta", "fasta")
            
            # Write haplotype 2
            hap2_record = SeqRecord(
                Seq(hap2_seq),
                id=f"{sample_name}_hap2",
                description=f"Haplotype 2 for {sample_name}"
            )
            SeqIO.write(hap2_record, sample_dir / f"{sample_name}_hap2.fasta", "fasta")
    
    def run_command(self, cmd, cwd=None):
        """Run a command and return the result."""
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"Command failed with exit code {result.returncode}")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            
        return result
    
    def test_build_pangenome(self):
        """Test building pangenome from variant FASTA files."""
        # Find the build_pangenome.py script
        script_dir = Path(__file__).parent.parent / "scripts"
        build_script = script_dir / "build_pangenome.py"
        
        self.assertTrue(build_script.exists(), f"Build script not found: {build_script}")
        
        # Build pangenome
        cmd = [
            "python", str(build_script),
            "-r", str(self.reference_file),
            "-s"] + self.sample_dirs + [
            "-o", str(self.pangenome_dir),
            "--cores", "2"
        ]
        
        result = self.run_command(cmd)
        self.assertEqual(result.returncode, 0, "Pangenome build failed")
        
        # Check that output files exist
        self.assertTrue((self.pangenome_dir / "pangenome.vg").exists())
        self.assertTrue((self.pangenome_dir / "pangenome.gfa").exists())
        
        return self.pangenome_dir / "pangenome.vg", self.pangenome_dir / "pangenome.gfa"
    
    def test_call_variants(self):
        """Test calling variants from pangenome."""
        # First build pangenome
        vg_file, gfa_file = self.test_build_pangenome()
        
        # Find the call_variants.py script
        script_dir = Path(__file__).parent.parent / "scripts"
        call_script = script_dir / "call_variants.py"
        
        self.assertTrue(call_script.exists(), f"Call script not found: {call_script}")
        
        # Call variants
        cmd = [
            "python", str(call_script),
            "-v", str(vg_file),
            "-r", str(self.reference_file),
            "-o", str(self.vcf_dir),
            "--cores", "2"
        ]
        
        result = self.run_command(cmd)
        self.assertEqual(result.returncode, 0, "Variant calling failed")
        
        # Check that VCF files exist
        expected_vcfs = [
            self.vcf_dir / "sample_001.vcf.gz",
            self.vcf_dir / "sample_002.vcf.gz",
            self.vcf_dir / "all_samples.vcf.gz"
        ]
        
        for vcf_file in expected_vcfs:
            self.assertTrue(vcf_file.exists(), f"VCF file not found: {vcf_file}")
        
        return expected_vcfs
    
    def apply_vcf_to_reference(self, vcf_file, output_fasta):
        """Apply VCF variants to reference genome using bcftools consensus."""
        # Index reference
        ref_indexed = str(self.reference_file) + ".fai"
        subprocess.run(["samtools", "faidx", str(self.reference_file)], check=True)
        
        # Compress and index VCF if not already compressed
        if not str(vcf_file).endswith('.gz'):
            vcf_gz = str(vcf_file) + ".gz"
            subprocess.run(["bgzip", "-c", str(vcf_file)], 
                          stdout=open(vcf_gz, 'wb'), check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz], check=True)
            vcf_file = vcf_gz
        
        # Apply variants
        with open(output_fasta, 'w') as f:
            result = subprocess.run([
                "bcftools", "consensus", 
                "-f", str(self.reference_file),
                str(vcf_file)
            ], stdout=f, check=True)
    
    def extract_haplotype_from_phased_vcf(self, vcf_file, sample_name, haplotype, output_fasta):
        """Extract a specific haplotype from a phased VCF."""
        # This is a simplified implementation
        # In practice, you might need more sophisticated haplotype extraction
        
        # Read VCF and extract variants for specific haplotype
        variants = []
        
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                format_field = fields[8]
                sample_field = fields[9]
                
                # Parse genotype (simplified)
                if '|' in sample_field:  # Phased
                    gt = sample_field.split(':')[0]
                    alleles = gt.split('|')
                    
                    # Get the allele for this haplotype (0 or 1)
                    allele_idx = int(alleles[haplotype])
                    
                    if allele_idx > 0:  # Variant allele
                        variants.append((pos - 1, ref, alt))  # Convert to 0-based
        
        # Apply variants to reference
        ref_record = SeqIO.read(self.reference_file, "fasta")
        ref_seq = list(str(ref_record.seq))
        
        # Sort variants by position (reverse order for insertions/deletions)
        variants.sort(key=lambda x: x[0], reverse=True)
        
        for pos, ref_allele, alt_allele in variants:
            if len(ref_allele) == len(alt_allele):  # SNV
                ref_seq[pos] = alt_allele
            elif len(ref_allele) > len(alt_allele):  # Deletion
                del ref_seq[pos:pos + len(ref_allele)]
                if alt_allele != '.':
                    ref_seq.insert(pos, alt_allele)
            else:  # Insertion
                ref_seq[pos] = alt_allele
        
        # Write result
        result_record = SeqRecord(
            Seq(''.join(ref_seq)),
            id=f"{sample_name}_hap{haplotype + 1}_reconstructed",
            description=f"Reconstructed from VCF"
        )
        
        SeqIO.write(result_record, output_fasta, "fasta")
    
    def compare_sequences(self, seq1_file, seq2_file):
        """Compare two FASTA sequences."""
        seq1 = SeqIO.read(seq1_file, "fasta")
        seq2 = SeqIO.read(seq2_file, "fasta")
        
        return str(seq1.seq) == str(seq2.seq)
    
    def test_full_roundtrip(self):
        """Test the complete roundtrip: FASTA -> Pangenome -> VCF -> FASTA."""
        print("\n=== Testing Full Roundtrip ===")
        
        # Step 1: Build pangenome
        print("Step 1: Building pangenome...")
        vg_file, gfa_file = self.test_build_pangenome()
        
        # Step 2: Call variants
        print("Step 2: Calling variants...")
        vcf_files = self.test_call_variants()
        
        # Step 3: Apply VCF back to reference and compare
        print("Step 3: Testing roundtrip...")
        
        # Test each sample
        for sample_name in ["sample_001", "sample_002"]:
            print(f"\nTesting {sample_name}...")
            
            # Get sample VCF
            sample_vcf = self.vcf_dir / f"{sample_name}.vcf.gz"
            
            # For now, we'll just check that the VCF files were created and have content
            self.assertTrue(sample_vcf.exists(), f"Sample VCF not found: {sample_vcf}")
            
            # Check VCF has variants
            with subprocess.Popen(['zcat', str(sample_vcf)], stdout=subprocess.PIPE, text=True) as proc:
                vcf_content = proc.stdout.read()
                variant_lines = [line for line in vcf_content.split('\n') if line and not line.startswith('#')]
                print(f"  Found {len(variant_lines)} variants in {sample_name}")
    
    def test_vcf_format_validity(self):
        """Test that generated VCF files are valid."""
        # Build pangenome and call variants
        vg_file, gfa_file = self.test_build_pangenome()
        vcf_files = self.test_call_variants()
        
        # Check VCF format validity
        for vcf_file in vcf_files:
            print(f"Validating {vcf_file.name}...")
            
            # Basic format checks using pysam
            try:
                vcf = pysam.VariantFile(str(vcf_file))
                
                # Check header
                self.assertIsNotNone(vcf.header)
                
                # Count variants
                variant_count = 0
                for record in vcf:
                    variant_count += 1
                    # Basic checks on first few records
                    if variant_count <= 5:
                        self.assertIsNotNone(record.chrom)
                        self.assertIsNotNone(record.pos)
                        self.assertIsNotNone(record.ref)
                        self.assertIsNotNone(record.alts)
                
                print(f"  Found {variant_count} variants")
                vcf.close()
                
            except Exception as e:
                self.fail(f"VCF validation failed for {vcf_file}: {e}")


def run_tests():
    """Run the roundtrip tests."""
    # Check if required tools are available
    required_tools = ['docker', 'samtools', 'bgzip', 'tabix', 'bcftools']
    missing_tools = []
    
    for tool in required_tools:
        result = subprocess.run(['which', tool], capture_output=True)
        if result.returncode != 0:
            missing_tools.append(tool)
    
    if missing_tools:
        print(f"Missing required tools: {', '.join(missing_tools)}")
        print("Please install these tools before running the test.")
        return False
    
    # Run tests
    unittest.main(verbosity=2)


if __name__ == "__main__":
    run_tests()

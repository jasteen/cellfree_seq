'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

import os
import math
from runner import run_stage


def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 2
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)


def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)


class Stages(object):
    def __init__(self, state):
        self.state = state
        # Reference genome and interval files
        self.reference = self.get_options('ref_grch37')
        self.interval_file = self.get_options('interval_file')
        self.interval_file_QC = self.get_options('interval_file_QC')
        self.AF_THR = self.get_options('AF_THR')
        self.vardict_bed = self.get_options('vardict_bed')
        self.vardict_bed_QC = self.get_options('vardict_bed_QC')
        # Programs and program settings
        self.other_vep = self.get_options('other_vep')
        # Annotation resources
        self.dbsnp_b37 = self.get_options('dbsnp_b37')
        self.brcaex = self.get_options('vep_brcaex')
        self.gnomad = self.get_options('vep_gnomad')
        self.revel = self.get_options('vep_revel')
        self.maxentscan = self.get_options('vep_maxentscan')
        self.exac = self.get_options('vep_exac')
        self.dbnsfp = self.get_options('vep_dbnsfp')
        self.dbscsnv = self.get_options('vep_dbscsnv')
        self.cadd = self.get_options('vep_cadd')

    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)

    def original_fastqs(self, output):
        '''grab the fastq files to map'''
        pass

    def align_bwa(self, inputs, bam_out, sample_id, lib):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        read_group = '"@RG\\tID:{sample}\\tSM:{sample}\\tPU:lib1\\tPL:Illumina"' \
            .format(sample=sample_id)
        
        command = 'bwa mem -M -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -u -h -q 1 -f 2 -F 4 -F 8 -F 256 - ' \
                  '| samtools sort -@ {cores} -o {bam}; samtools index {bam}'.format(
                          cores=cores,
                          read_group=read_group,
                          fastq_read1=fastq_read1_in,
                          fastq_read2=fastq_read2_in,
                          reference=self.reference,
                          bam=bam_out)
        run_stage(self.state, 'align_bwa', command)

    def run_vardict(self, bam_in, vcf_out, sample_name):
        '''Run vardict on each BAM file'''
        
        if "QC" in bam_in: 
            vardict_bed = self.vardict_bed_QC 
        else:  
            vardict_bed = self.vardict_bed 

        command = 'export PATH=/home/jste0021/scripts/git_controlled/VarDict:$PATH; ' \
                  'vardict -G {reference} -f {AF_THR} -N {sample_name} -b {bam_in} -c 1 -S 2 -E 3 -g 4 {vardict_bed} | ' \
                  'teststrandbias.R | ' \
                  'var2vcf_valid_b37_chrnames.pl -N {sample_name} -E -f {AF_THR} > {vcf_out}'.format(
                             reference=self.reference,
                             AF_THR=self.AF_THR,
                             sample_name=sample_name,
                             bam_in=bam_in,
                             vardict_bed=vardict_bed,
                             vcf_out=vcf_out)
        run_stage(self.state, 'run_vardict', command)
                               
    def sort_vcfs(self, vcf_in, vcf_out):
        '''sort undr_rover vcf files'''
        command = 'bcftools sort -o {vcf_out} -O z {vcf_in}'.format(vcf_out=vcf_out, vcf_in=vcf_in)
        run_stage(self.state, 'sort_vcfs', command)

    def index_vcfs(self, vcf_in, vcf_out):
        '''Index undr_rover vcf files'''
        command = 'bcftools index -f --tbi {vcf_in}'.format(vcf_in=vcf_in)
        run_stage(self.state, 'index_vcfs', command)

    # Write as a collate
    def concatenate_vcfs(self, vcf_files_in, vcf_out):
        merge_commands = []
        temp_merge_outputs = []
        for n in range(0, int(math.ceil(float(len(vcf_files_in)) / 200.0))):
            start = n * 200
            filelist = vcf_files_in[start:start + 200]
            filelist_command = ' '.join([vcf for vcf in filelist])
            temp_merge_filename = vcf_out.rstrip('.vcf') + ".temp_{start}.vcf".format(start=str(start))
            command1 = 'bcftools merge -O z -o {vcf_out} {join_vcf_files} && bcftools index -t -f {vcf_out}; '.format(vcf_out=temp_merge_filename, join_vcf_files=filelist_command)
            merge_commands.append(command1)
            temp_merge_outputs.append(temp_merge_filename)

        final_merge_vcfs = ' '.join([vcf for vcf in temp_merge_outputs])
        command2 = 'bcftools merge -O z -o {vcf_out} {join_vcf_files} '.format(vcf_out=vcf_out, join_vcf_files=final_merge_vcfs)

        merge_commands.append(command2)
        final_command = ''.join(merge_commands)
        run_stage(self.state, 'concatenate_vcfs', final_command)

    def vt_decompose_normalise(self, vcf_in, vcf_out):
        '''Decompose multiallelic sites and normalise representations'''
        command = "vt decompose -s {vcf_in} | vt normalize -r {reference} -o " \
                  "{vcf_out} -".format(reference=self.reference,
                                       vcf_in=vcf_in,
                                       vcf_out=vcf_out)
        run_stage(self.state, 'vt_decompose_normalise', command)

    def apply_vep(self, vcf_in, vcf_out):
        '''Apply VEP'''
        cores = self.get_stage_options('apply_vep', 'cores')
        vep_command = "vep --cache --dir_cache {other_vep} " \
                      "--assembly GRCh37 --refseq --offline " \
                      "--fasta {reference} " \
                      "--sift b --polyphen b --symbol --numbers --biotype " \
                      "--total_length --hgvs --format vcf " \
                      "--vcf --force_overwrite --flag_pick --no_stats " \
                      "--custom {brcaexpath},brcaex,vcf,exact,0,Clinical_significance_ENIGMA," \
                      "Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA," \
                      "Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature " \
                      "--custom {gnomadpath},gnomAD,vcf,exact,0,AF_NFE,AN_NFE " \
                      "--custom {revelpath},RVL,vcf,exact,0,REVEL_SCORE " \
                      "--plugin MaxEntScan,{maxentscanpath} " \
                      "--plugin ExAC,{exacpath},AC,AN " \
                      "--plugin dbNSFP,{dbnsfppath},REVEL_score,REVEL_rankscore " \
                      "--plugin dbscSNV,{dbscsnvpath} " \
                      "--plugin CADD,{caddpath} " \
                      "--fork {cores} " \
                      "-i {vcf_in} " \
                      "-o {vcf_out}".format(other_vep=self.other_vep,
                                            cores=cores,
                                            vcf_out=vcf_out,
                                            vcf_in=vcf_in,
                                            reference=self.reference,
                                            brcaexpath=self.brcaex,
                                            gnomadpath=self.gnomad,
                                            revelpath=self.revel,
                                            maxentscanpath=self.maxentscan,
                                            exacpath=self.exac,
                                            dbnsfppath=self.dbnsfp,
                                            dbscsnvpath=self.dbscsnv,
                                            caddpath=self.cadd)
        run_stage(self.state, 'apply_vep', vep_command)


    #### generate stats for mapping ####
    def intersect_bed(self, bam_in, bam_out):
        '''intersect the bed file with the interval file '''
        
        if "QC" in bam_in: 
            interval_file = self.interval_file_QC 
        else:  
            interval_file = self.interval_file

        command = "intersectBed -abam {bam_in} -b {interval_file} > {bam_out} " \
                .format(bam_in=bam_in,
                        interval_file=interval_file,
                        bam_out=bam_out)
        run_stage(self.state, 'intersect_bed', command)

    def coverage_bed(self, bam_in, txt_out):
        ''' make coverage files '''
        
        if "QC" in bam_in:
            interval_file = self.interval_file_QC
        else:
            interval_file = self.interval_file

        command = "coverageBed -b {bam_in} -a {interval_file} -hist | grep all > {txt_out}" \
                .format(bam_in=bam_in,
                        interval_file=interval_file,
                        txt_out=txt_out)
        run_stage(self.state, 'coverage_bed', command)

    def genome_reads(self, bam_in, txt_out):
        '''count reads that map to the genome'''
        command = 'samtools view -c -F4 {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'genome_reads', command)

    def target_reads(self, bam_in, txt_out):
        '''count reads that map to target panel'''
        command = 'samtools view -c -F4 {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'target_reads', command)

    def total_reads(self, bam_in, txt_out):
        '''count the total number of reads that we started with'''
        command = 'samtools view -c {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'total_reads', command)
    
    # Try and get rid of the R script
    def generate_stats(self, inputs, txt_out, samplename, joint_output):
        '''run R stats script'''
        # Assigning inputfiles to correct variables based on suffix
        for inputfile in inputs:
            if inputfile.endswith('.bedtools_hist_all.txt'):
                a = inputfile
            elif inputfile.endswith('.mapped_to_genome.txt'):
                b = inputfile
            elif inputfile.endswith('.mapped_to_target.txt'):
                c = inputfile
            elif inputfile.endswith('.total_raw_reads.txt'):
                d = inputfile
        e = samplename
        command = 'touch {txt_out};  Rscript --vanilla /projects/vh83/pipelines/code/modified_summary_stat.R ' \
                  '{hist_in} {map_genome_in} {map_target_in} {raw_reads_in} {sample_name} ' \
                  '{summary_out}'.format(txt_out=txt_out,
                                      hist_in=a,
                                      map_genome_in=b,
                                      map_target_in=c,
                                      raw_reads_in=d ,
                                      sample_name=e ,
                                      summary_out=joint_output)
        run_stage(self.state, 'generate_stats', command)


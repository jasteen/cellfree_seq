'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from, regex
from stages import Stages
import glob
from utils import safe_make_dir

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='hiplexpipe_vardict')
    # Stages are dependent on the state
    stages = Stages(state)

    safe_make_dir('alignments')
    safe_make_dir('metrics')
    safe_make_dir('metrics/summary')

    # The original FASTQ files
    fastq_files = glob.glob('fastqs/*')

    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_R2_{lib[0]}.fastq.gz'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}', '{lib[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}.sort.hq.bam')

    pipeline.transform(
        task_func=stages.intersect_bed,
        name='intersect_bed',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='metrics/summary/{sample[0]}.intersectbed.bam')

    pipeline.transform(
        task_func=stages.coverage_bed,
        name='coverage_bed',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.bedtools_hist_all.txt')

    pipeline.transform(
        task_func=stages.genome_reads,
        name='genome_reads',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='metrics/summary/{sample[0]}.mapped_to_genome.txt')

    pipeline.transform(
        task_func=stages.target_reads,
        name='target_reads',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.mapped_to_target.txt')

    pipeline.transform(
        task_func=stages.total_reads,
        name='total_reads',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='metrics/summary/{sample[0]}.total_raw_reads.txt')

    pipeline.collate(
        task_func=stages.generate_stats,
        name='generate_stats',
        input=output_from('coverage_bed', 'genome_reads', 'target_reads', 'total_reads'),
        filter=regex(r'.+/(.+)\.(bedtools_hist_all|mapped_to_genome|mapped_to_target|total_raw_reads)\.txt'),
        output=r'metrics/summary/all_sample.summary.\1.txt',
        extras=[r'\1', 'all_sample.summary.txt'])
    
    safe_make_dir('variants')
    safe_make_dir('variants/vardict')
   
    pipeline.transform(
        task_func=stages.run_vardict,
        name='run_vardict',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='variants/vardict/{sample[0]}.vcf',
        extras=['{sample[0]}'])   
 
    pipeline.transform(
        task_func=stages.sort_vcfs,
        name='sort_vcfs',
        input=output_from('run_vardict'),
        filter=formatter('variants/vardict/(?P<sample>[a-zA-Z0-9_-]+).vcf'),
        output='variants/vardict/{sample[0]}.sorted.vcf.gz')

    pipeline.transform(
        task_func=stages.index_vcfs,
        name='index_vcfs',
        input=output_from('sort_vcfs'),
        filter=suffix('.sorted.vcf.gz'),
        output='.sorted.vcf.gz.tbi')

    (pipeline.merge(
        task_func=stages.concatenate_vcfs,
        name='concatenate_vcfs',
        input=output_from('sort_vcfs'),
        output='variants/vardict/combined.vcf.gz')
        .follows('index_vcfs'))
 
    pipeline.transform(
        task_func=stages.vt_decompose_normalise,
        name='vt_decompose_normalise',
        input=output_from('concatenate_vcfs'),
        filter=suffix('.vcf.gz'),
        output='.decomp.norm.vcf.gz')
      
    pipeline.transform(
        task_func=stages.index_vcfs,
        name='index_final_vcf',
        input=output_from('vt_decompose_normalise'),
        filter=suffix('.decomp.norm.vcf.gz'),
        output='.decomp.norm.vcf.gz.tbi')
   
    (pipeline.transform(
        task_func=stages.apply_vep,
        name='apply_vep',
        input=output_from('vt_decompose_normalise'),
        filter=suffix('.decomp.norm.vcf.gz'),
        output='.decomp.norm.vep.vcf')
        .follows('index_final_vcf'))

    return pipeline

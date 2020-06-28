
from gwf import Workflow 
from pathlib import Path
import os, string, re

#################################################################################
# Templates
#################################################################################


def intervals_to_sites(intervals_file, sites_file): 

    options = {'memory': '4g',
               'walltime': '11:00:00'
               }

    spec = """
    python ./scripts/liftover_sites.py {intervals_file} {sites_file}

    """

    return [str(intervals_file)], [str(sites_file)], options, spec


def liftover(bed_file, chain_file, mapped_file, unmapped_file):

    options = {'memory': '4g',
               'walltime': '36:00:00'
               }

    spec = """
    liftOver -bedPlus=3 {bed_file} {chain_file} {mapped_file} {unmapped_file}

    """.format(bed_file=bed_file, chain_file=chain_file, mapped_file=mapped_file, unmapped_file=unmapped_file)

    return [str(bed_file), str(chain_file)], [str(mapped_file), str(unmapped_file)], options, spec


def bed_difference(bed_file1, bed_file2, output_file): 

    options = {'memory': '4g',
               'walltime': '11:00:00'
               }
    shell_spec = """
    python ./scripts/liftover_diff.py {bed_file1} {bed_file2} | sort -k1,1 -k2,2n -k3,3n > {output_file}

    """.format(bed_file1=bed_file1, bed_file2=bed_file2, output_file=output_file)

    return [str(bed_file1), str(bed_file2)], [str(output_file)], options, shell_spec


def bed_merge_and_split(input_files, output_files):

    output_base_names = [Path(x).name for x in output_files]

    # get output dir and make sure it is unambiguous:
    output_dirs = [Path(x).parent for x in output_files]
    assert len(set(output_dirs)) == 1
    output_dir = output_dirs[0]

    input_files = [str(x) for x in input_files]
    output_files = [str(x) for x in output_files]

    options = {'memory': '10g',
               'walltime': '11:00:00'}

    shell_spec = """
    sort -k1,1 -k2,2n -k3,3n --merge {input_files} -T /scratch/$GWF_JOBID | python ./scripts/bed_split.py {output_dir} {output_base_names}

    """.format(input_files=" ".join(input_files),
               output_dir=output_dir,
               output_base_names=" ".join(output_base_names))

    return input_files, output_files, options, shell_spec


#################################################################################
# Workflow components
#################################################################################

gwf = Workflow(defaults={'account': 'simons'})


def split_file(input_file, split_files_dir, n_files=50):
    """
    Split a file into n_files chunks.
    """
    if not split_files_dir.exists():
        os.makedirs(str(split_files_dir))
    prefix = str(split_files_dir / input_file.with_suffix('').name)
    split_files = [Path('{}.{}{}'.format(prefix, x, input_file.suffix)) for x in range(n_files)]
    gwf.target('split_files', inputs=[str(input_file)], outputs=[str(p) for p in split_files]) << """

        python ./scripts/file_chunks.py $((1 + `wc -l < {input_file}`/{n_files})) {prefix} {suffix} {input_file}

        """.format(input_file=str(input_file), prefix=prefix, suffix=input_file.suffix, n_files=n_files)
    return split_files



def reciprocal_liftover(intervals_files, forwards_chain_file, backwards_chain_file, 
                        slurm_tag, steps_dir, target_chromosomes):
    """
    Does reciprocal lift over of a set of intervals.
    """

    if not steps_dir.exists():
        os.makedirs(str(steps_dir))

    # output files
    mapped_files= [steps_dir / x.with_suffix('.mapped').name for x in intervals_files]
    unmapped_files = [x.with_suffix('.unmapped') for x in mapped_files]
    backmapped_files = [x.with_suffix('.backmapped') for x in mapped_files]
    unbackmapped_files = [x.with_suffix('.nobackmapped') for x in mapped_files]
    filtered_files = [x.with_suffix('.filtered') for x in mapped_files]

    lifted_files = [steps_dir / 'sorted' / "{}.bed".format(x) for x in target_chromosomes]

    for i, chrom in enumerate(intervals_files):

        # lift over intervals
        gwf.target_from_template('{}_lift_{}'.format(slurm_tag, i),
            liftover(bed_file=intervals_files[i], chain_file=forwards_chain_file,
            mapped_file=mapped_files[i], unmapped_file=unmapped_files[i]))

        # lift back to orginal coordinates to ensure one to one correspondence
        gwf.target_from_template('{}_liftback_{}'.format(slurm_tag, i),
            liftover(bed_file=mapped_files[i], chain_file=backwards_chain_file,
            mapped_file=backmapped_files[i], unmapped_file=unbackmapped_files[i]))

        # filter out intervals that does not map both ways
        gwf.target_from_template('{}_filter_{}'.format(slurm_tag, i),
            bed_difference(bed_file1=mapped_files[i], bed_file2=unbackmapped_files[i],
            output_file=filtered_files[i]))

    # filter out intervals that does not map both ways
    gwf.target_from_template('{}_merge_and_split'.format(slurm_tag),
        bed_merge_and_split(input_files=filtered_files, output_files=lifted_files))

    return lifted_files


#################################################################################
# Workflow
#################################################################################

human_chromosomes = ['chr{}'.format(x) for x in list(range(1,23)) + ['X']]
#ape_chromosomes = ['chr1', 'chr2a', 'chr2b'] + ['chr{}'.format(x) for x in list(range(3,23)) + ['X']]

split_bed_files = split_file(Path('decode_hg38_sexavg_per_gen.tsv'), Path('steps/split_files'), n_files=50)

sites_files_lifted_to_species = reciprocal_liftover(split_bed_files, 
    forwards_chain_file=Path('chain_files/hg38ToHg19.over.chain'), 
    backwards_chain_file=Path('chain_files/hg19ToHg38.over.chain'),
    slurm_tag='siteslift',
    steps_dir=Path(os.getcwd(), 'steps', 'liftover'),
    target_chromosomes=human_chromosomes)


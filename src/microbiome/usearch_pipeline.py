from loguru import logger 
from pathlib import Path
import os 
from subprocess import run
from joblib import Parallel, delayed

from typing import Literal


class UsearchPipeline:
    '''Usearch pipeline for 16S rRNA amplicon data analysis.
    USEARCH11 (https://www.drive5.com/usearch/manual/),
    USEARCH12 (https://rcedgar.github.io/usearch12_documentation/cmds.html),
    Including quality control, merging paired-end reads, filtering low-quality reads, dereplication, removing singletons, 
    preorder, clustering to OTUs, denoising to zOTUs, building OTU and zOTU feature tables.
    '''

    def __init__(self, home_dir: str, threads: int | None = None, db_file: str | None = None, usearch_version: str = 'usearch12'):
        '''Initialize the pipeline with home directory and number of threads.
        Args:
            home_dir: str, the home directory for the pipeline, which should contain a 'datas' directory with input files and will be used to store results in a 'results' directory.
                |__ datas (Requires*)
                    |__ sample_1_1.fastq/fq
                    |__ sample_1_2.fastq/fq
                    |__ sample_2_1.fastq/fq
                    |__ sample_2_2.fastq/fq
                |__ results (Optional, will be created if not exist)
                    |__ 001_quality_control
                    |__ 002_merge_paired_end_reads
                    |__ 003_all_samples_merged
                        |__ all_samples.fastq
                        |__ all_samples_filtered.fastq
                        |__ all_samples_filtered_dereplicated.fasta
                        |__ all_samples_filtered_dereplicated_no_singleton.fasta
                        |__ all_samples_filtered_dereplicated_no_singleton_preorder.fasta
                        |__ all_samples_filtered_dereplicated_no_singleton_preorder_otus.fasta
                        |__ all_samples_filtered_dereplicated_no_singleton_preorder_zotus.fasta
                        |__ all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table.txt
                        |__ all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table_map.txt
                        |__ all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table.txt
                        |__ all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table_map.txt
            threads: int, the number of threads to use for parallel processing, default is None which means using all available CPU cores.
            db_file: str, the path to the database file for taxonomic classification.
        '''

        self.usearch_version: str = usearch_version
        self.home_dir: Path = Path(home_dir)
        self.data_dir: Path = self.home_dir.joinpath('datas')
        self.results_dir: Path = self.home_dir.joinpath('results')
        self.threads: int = threads if threads is not None else os.cpu_count()
        self.db_file: str | None = db_file

        logger.info(f'Home directory set to {home_dir}')
        if not self.results_dir.exists():
            self.results_dir.mkdir(parents=True)
            logger.info(f'Created results directory at {self.results_dir}')

    def check_directories(self, dir: Path, notices: str):
        '''Check directory then create one if not exist.'''
        if not dir.exists():
            dir.mkdir(parents=True)
            logger.info(notices)

    def run_cmd(self, cmd: str, show: bool = False):
        '''Run bash commands.'''
        try:
            if isinstance(cmd, list):
                cmd = ' '.join(cmd) 
            run(cmd, shell=True, check=True, capture_output=not show, text=True)
        except Exception as e:
            logger.error(f"Command failed: {cmd}")
            logger.error(f"Error: {e}")

    def show_files(self):
        logger.info(f'Data directory is set to: {self.data_dir}')
        for read_file in self.data_dir.iterdir():
            # print(read_file)
            pass

    def quality_control(self):
        logger.info('Quality control step starting...')
        qc_dir: Path = self.results_dir.joinpath('001_quality_control')
        self.check_directories(qc_dir, f'Created quality control directory at {qc_dir}')

    def merge_paired_end_reads(self, connect: Literal['_', '.'] = '.', fastq_format: Literal['fastq', 'fq'] = 'fastq'):
        '''Merge all paired reads.
        Examples:
        >>> usearch11 -fastq_mergepairs sample_1.fastq -reverse sample_2.fastq -fastqout sample_merged.fastq -relabel @
        '''
        
        def get_cmds(file_dir: Path, file_name:str, output_dir: Path, connect: Literal['_', '.'] = '_', fastq_format: Literal['fastq', 'fq'] = 'fastq') -> str | None:
            '''Generate command for merging paired-end reads.'''
            file_1, file_2 = tuple(f"{file_dir}/{file_name}{connect}{i}.{fastq_format}" for i in [1, 2])
            merge_file = output_dir.joinpath(f"{file_name}.fastq")

            # If the forward and reverse files do not exist, skip merging for this sample.
            if Path(file_1).exists() and Path(file_2).exists():
                return f"{self.usearch_version} -fastq_mergepairs {file_1} -reverse {file_2} -fastqout {merge_file} -relabel @" ## @ 很重要
            else:
                return None 
        
        logger.info('Merging paired-end reads...')
        merge_dir: Path = self.results_dir.joinpath('002_merge_paired_end_reads')
        self.check_directories(merge_dir, f'Created merge directory at {merge_dir}')
        file_name_list = {file.stem.split(connect)[0] for file in self.data_dir.iterdir() if file.suffix in {'.fastq', '.fq'}} # Get the unique sample names.
        cmds, paired_nums, single_nums = [], 0, 0
        for file_name in file_name_list:
            cmd = get_cmds(self.data_dir, file_name, merge_dir, connect, fastq_format)
            if cmd is not None:
                cmds.append(cmd)                                                                                                # Paired-end files.
                paired_nums += 1
            else:
                cp_cmd: str = f'cp {self.data_dir.joinpath(f"{file_name}*.{fastq_format}")} {merge_dir}'
                cmds.append(cp_cmd)                                                                                             # Single-end files, just copy to the merge directory.
                single_nums += 1
        tasks = [delayed(self.run_cmd)(cmd) for cmd in cmds]
        parallel = Parallel(n_jobs=self.threads, verbose=50)
        parallel(tasks)

        if single_nums > 0:
            logger.warning(f'Finished merging paired-end reads. Total paired samples: {paired_nums}, total single samples: {single_nums}.')

    def merge_all_samples(self):
        '''Merge all samples reads into one file.
        Examples:
        >>> cat 002_merge_paired_end_reads/*.fastq > 003_all_samples_merged/all_samples.fastq
        >>> echo "Total reads in all samples:" && grep @ 003_all_samples_merged/all_samples.fastq | wc -l
        '''

        logger.info('Merging all samples reads into one file...')
        merge_dir: Path = self.results_dir.joinpath('002_merge_paired_end_reads')
        all_samples_dir: Path = self.results_dir.joinpath('003_all_samples_merged')
        self.check_directories(all_samples_dir, f'Created directory for merged samples at {all_samples_dir}')
        merge_cmd = f'cat {merge_dir}/*q > {all_samples_dir.joinpath("all_samples.fastq")}'                             # fastq or fq.
        self.run_cmd(merge_cmd)
        self.run_cmd(f'echo "Total reads in all samples:" && grep @ {all_samples_dir.joinpath("all_samples.fastq")} | wc -l', show=True)

    def filter_low_quality_reads(self, fastq_maxee: float = 1.0):
        '''Filter low-quality merged reads (according to the maximum expected errors).
        Args:
            fastq_maxee: float, the maximum expected errors allowed for a read to pass the filter, default is 1.0.
        Return:
            None
        Examples:
        >>> usearch11 -fastq_filter all_samples.fastq -fastaout all_samples_filtered.fastq -fastq_maxee 1.0 -relabel Filt -threads 8
        '''

        logger.info(f'Filtering low-quality reads step starting with maximum expected errors of {fastq_maxee} ...')
        all_samples_file: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples.fastq')
        all_samples_filtered: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered.fastq')
        filter_cmd = f'{self.usearch_version} -fastq_filter {all_samples_file} -fastaout {all_samples_filtered} -fastq_maxee {fastq_maxee} -relabel Filt -threads {self.threads}'
        self.run_cmd(filter_cmd)        
        self.run_cmd(f'echo "Total reads after filtering:" && grep ">" {all_samples_filtered} | wc -l', show=True)

    def dereplication(self):
        '''Dereplicate merged and filtered reads.
        Examples:
        >>> usearch11 -fastx_uniques all_samples_filtered.fastq -fastaout all_samples_filtered_dereplicated.fasta -relabel Uniq -sizeout -threads 8
        '''

        logger.info('Dereplication step starting...')
        all_samples_filtered: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered.fastq')
        derep_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated.fasta')
        derep_cmd = f'{self.usearch_version} -fastx_uniques {all_samples_filtered} -fastaout {derep_fasta} -relabel Uniq -sizeout -threads {self.threads}'
        self.run_cmd(derep_cmd)        
        self.run_cmd(f'echo "Total unique reads after dereplication:" && grep ">" {derep_fasta} | wc -l', show=True)

    def remove_singletons(self, minsize: int = 2):
        '''Remove singletons from the dereplicated fasta file.
        Args:
            minsize: int=2, the minimum size of reads to keep, default is 2, which means removing singletons.
        Return:
            None
        Examples:
        >>> usearch11 -sortbysize all_samples_filtered_dereplicated.fasta -fastaout all_samples_filtered_dereplicated_no_singleton.fasta -minsize 2
        '''

        logger.info('Remove singletons step starting...')
        derep_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated.fasta')
        no_singleton_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton.fasta')
        remove_cmd = f'{self.usearch_version} -sortbysize {derep_fasta} -fastaout {no_singleton_fasta} -minsize {minsize} '
        self.run_cmd(remove_cmd)        
        self.run_cmd(f'echo "Total unique reads after removing singletons:" && grep ">" {no_singleton_fasta} | wc -l', show=True)

    def preorder(self, minsize: int = 8):
        '''Preorder sequences by abundance, and remove low-abundance reads accroding to the abundance.
        Args:
            minsize: int=8, the minimum size of reads to keep, default is 8, which means removing low-abundance reads.
        Return:
            None
        Examples:
        >>> usearch11 -sortbysize all_samples_filtered_dereplicated_no_singleton.fasta -fastaout all_samples_filtered_dereplicated_no_singleton_preorder.fasta -minsize 8
        '''

        logger.info('Preorder step starting...')
        no_singleton_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton.fasta')
        preorder_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder.fasta')
        preorder_cmd = f'{self.usearch_version} -sortbysize {no_singleton_fasta} -fastaout {preorder_fasta} -minsize {minsize}'
        self.run_cmd(preorder_cmd)
        self.run_cmd(f'echo "Total OTUs after preorder:" && grep ">" {preorder_fasta} | wc -l', show=True)

    def cluster_to_outs(self, method: Literal['usearch', 'vsearch'] = 'usearch'):
        '''Cluster sequences into OTUs.
        Args:
            method: str, the method to use for clustering, default is 'usearch', can also be 'vsearch'.
        Return:
            None
        Examples:
        >>> usearch11 -cluster_otus all_samples_filtered_dereplicated_no_singleton_preorder.fasta -otus all_samples_filtered_dereplicated_no_singleton_preorder_otus.fasta -relabel Otu -threads 8 > all_samples_filtered_dereplicated_no_singleton_preorder_otus.log 2>&1
        '''

        logger.info('Clustering sequences into OTUs step starting...')
        preorder_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder.fasta')
        otu_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus.fasta')
        otu_fasta_log: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus.log')
        if method == 'usearch':
            cluster_cmd = f'{self.usearch_version} -cluster_otus {preorder_fasta} -otus {otu_fasta} -relabel Otu -threads {self.threads} > {otu_fasta_log}.{method} 2>&1'
        elif method == 'vsearch':
            cluster_cmd = f'vsearch --cluster_size {preorder_fasta} --centroids {otu_fasta} --log {otu_fasta_log}.{method}'
        else:
            raise ValueError(f"Invalid method: {method}. Must be 'usearch' or 'vsearch'.")
        self.run_cmd(cluster_cmd)
        self.run_cmd(f'echo "Total OTUs after clustering:" && grep ">" {otu_fasta} | wc -l', show=True)

    def denoise_to_zotus(self, method: Literal['usearch', 'vsearch'] = 'usearch'):
        '''Denoise sequences to zOTUs.
        Args:
            method: str, the method to use for denoising, default is 'usearch', can also be 'vsearch'.
        Return:
            None
        Examples:
        >>> usearch11 -unoise3 all_samples_filtered_dereplicated_no_singleton_preorder.fasta -zotus all_samples_filtered_dereplicated_no_singleton_preorder_zotus.fasta -threads 8 > all_samples_filtered_dereplicated_no_singleton_preorder_zotus.log 2>&1
        '''

        logger.info('Denoising sequences to obtain zOTUs step starting...')
        no_singleton_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder.fasta')
        denoised_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus.fasta')
        denoised_fasta_log: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus.log')
        if method == 'usearch':
            denoise_cmd = f'{self.usearch_version} -unoise3 {no_singleton_fasta} -zotus {denoised_fasta} -threads {self.threads} > {denoised_fasta_log}.{method} 2>&1'
        elif method == 'vsearch':
            denoise_cmd = f'vsearch --cluster_size {no_singleton_fasta} --centroids {denoised_fasta} --log {denoised_fasta_log}.{method}'
        else:
            raise ValueError(f"Invalid method: {method}. Must be 'usearch' or 'vsearch'.")
        self.run_cmd(denoise_cmd)
        self.run_cmd(f'echo "Total ZOTUs after denoising:" && grep ">" {denoised_fasta} | wc -l', show=True)

    def build_otu_feature_table(self):
        '''Building OTU feature table.
        Examples:
        >>> usearch11 -otutab all_samples.fastq -otus all_samples_filtered_dereplicated_no_singleton_preorder_otus.fasta -otutabout all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table.txt -mapout all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table_map.txt -threads 8 > all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table.log 2>&1
        '''

        logger.info('Building OTU feature table step starting...')
        all_samples: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples.fastq')
        otu_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus.fasta')
        feature_table: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table.txt')
        feature_table_log: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table.log')
        feature_table_map: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus_feature_table_map.txt')
        build_cmd = f'{self.usearch_version} -otutab {all_samples} -otus {otu_fasta} -otutabout {feature_table} -mapout {feature_table_map} -threads {self.threads} > {feature_table_log} 2>&1'
        self.run_cmd(build_cmd)
        self.run_cmd(f'echo "First 10 lines of the feature table:" && head -n 10 {feature_table}', show=True)

    def build_zotu_feature_table(self):
        '''Building zOTU feature table.
        Examples: 
        >>> usearch11 -otutab all_samples.fastq -zotus all_samples_filtered_dereplicated_no_singleton_preorder_zotus.fasta -otutabout all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table.txt -mapout all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table_map.txt -threads 8 > all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table.log 2>&1
        '''

        logger.info('Building the zOTU feature table step starting...')
        all_samples: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples.fastq')
        zotu_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus.fasta')
        zotu_feature_table: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table.txt')
        zotu_feature_table_log: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table.log')
        zotu_feature_table_map: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus_feature_table_map.txt')
        build_cmd = f'{self.usearch_version} -otutab {all_samples} -zotus {zotu_fasta} -otutabout {zotu_feature_table} -mapout {zotu_feature_table_map} -threads {self.threads} > {zotu_feature_table_log} 2>&1'
        self.run_cmd(build_cmd)
        self.run_cmd(f'echo "First 10 lines of the zOTU feature table:" && head -n 10 {zotu_feature_table}', show=True)

    def classification(self, db_file: str | None = None):
        '''Classify OTUs or zOTUs using usearch.
        Args:
            db_file: str, the path to the database file for taxonomic classification, default is None which means using the default database file.
        Return:
            None
        Examples:
        >>> usearch11 -sintax all_samples_filtered_dereplicated_no_singleton_preorder_otus.fasta -db . ./dbs/rdp_16s_v18_usearch.udb -tabbedout all_samples_filtered_dereplicated_no_singleton_preorder_otus_sintax.txt -strand both 
        '''

        if db_file is None or not Path(db_file).exists():
            db_file = self.db_file if self.db_file is not None and Path(self.db_file).exists() else '../dbs/rdp_16s_v18_usearch.udb'
        db_file = Path(db_file)
        logger.info(f'Using database file for classification: {db_file}')

        logger.info('Classification step starting...')
        otu_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus.fasta')
        zotu_fasta: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus.fasta')
        otu_sintax: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_otus_sintax.txt')
        zotu_sintax: Path = self.results_dir.joinpath('003_all_samples_merged/all_samples_filtered_dereplicated_no_singleton_preorder_zotus_sintax.txt')
        for file in [(otu_fasta, otu_sintax), (zotu_fasta, zotu_sintax)]:
            if file[0].exists():
                logger.info(f'Classifying {file[0].name} ...')
                cmd = f'{self.usearch_version} -sintax {file[0]} -db {db_file} -tabbedout {file[1]} -strand both -threads {self.threads} > {file[1]}.log 2>&1'
                self.run_cmd(cmd)


    def run(self, connect: Literal['_', '.'] = '.', fastq_format: Literal['fastq', 'fq'] = 'fq'):
        '''One step to run all steps in the pipeline.
        Args:
            connect: str, the character connecting paired-end read files, default is '_', can also be '.'
            fastq_format: str, the format of the input fastq files, default is 'fq', can also be 'fastq'
            method: str, the version of usearch to use, default is 'usearch12', can also be 'usearch11'
        Return:
            None

        Examples:
        >>> up = UsearchPipeline(home_dir='/path/to/home/dir')
        >>> up.run()
        '''

        logger.info(f'Running the Usearch pipeline with {self.usearch_version} ...')
        
        # Quality control.
        self.show_files()
        self.quality_control()
        self.merge_paired_end_reads(connect=connect, fastq_format=fastq_format)
        self.merge_all_samples()
        self.filter_low_quality_reads()
        self.dereplication()

        # Remove singletons and Preorder
        temp_usearch_version = self.usearch_version # Temporarily switch to usearch11 for removing singletons and preorder, which is faster than usearch12, then switch back to usearch12 for clustering and denoising, which is faster than usearch11
        self.usearch_version = 'usearch11'          # Use usearch11 for removing singletons and preorder, which is faster than usearch12
        self.remove_singletons(minsize=2)           # Deprecate in Usearch12 version
        self.preorder(minsize=8)                    # Deprecate in Usearch12 version
        self.usearch_version = temp_usearch_version # Restore usearch version for clustering and denoising, which is faster in usearch12 than usearch11

        # Cluster or Denoise, then build feature table.
        # self.cluster_to_outs()
        self.denoise_to_zotus()
        # self.build_otu_feature_table()
        self.build_zotu_feature_table()

        # Classify OTUs or zOTUs
        self.classification()

        logger.info('Pipeline execution completed.')

        return None 


if __name__ == "__main__":
    up = UsearchPipeline(home_dir='/bmp/backup/zhaosy/ws/china_16s_pipeline/results_demo', db_file='../dbs/rdp_16s_v18_usearch.udb')
    up.run()
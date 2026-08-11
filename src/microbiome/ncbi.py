from loguru import logger
from pathlib import Path
from subprocess import run
import time

from typing import Literal

import pandas as pd 


def run_command(command: str, retries: int = 3, show: bool = True) -> None:
    '''Run command in shell and check for errors.
    Args:
        command: Command to run in shell.
        retries: Number of times to retry the command.
        show: Whether to display the command output.
    Returns:
        None
    '''

    if isinstance(command, list):
        shell = False
    elif isinstance(command, str):
        shell = True

    for i in range(retries):
        try:
            result = run(command, shell=shell, check=True, capture_output=not show, text=True)
            return 0
        
        except Exception as e:
            # logger.error(f'Error running command: {command}')
            # logger.error(e)
            logger.warning(f'Retrying {i+1}/{retries} after error: {e}')
            time.sleep(1)  # Wait before retrying

    logger.error(f'Failed to run command after {retries} attempts: {command}')

    return None


class SRATools:
    '''Download SRA data using Aspera connect according to Project ID and metadata file.
    Examples:
    >>> sratoools = SRATools()
    >>> meta_file = sratools.query(project_id='PRJNA123456', output_meta_file='./metadata.csv')
    >>> srr_ids, sequence_types = sratools.parse_meta(meta_file='./metadata.csv')
    >>> sratoools.download(srr_ids, sequence_types, output_dir='./data', method='asp', backend='threading', pre_dispath=5, verbose=10)
    '''

    def __init__(
        self,
        aspera_path: str | None = None,
        aspera_key_path: str | None = None,
    ):

        self.aspera_path = aspera_path if aspera_path is not None else '~/.aspera/connect/bin/ascp'
        self.aspera_key_path = aspera_key_path if aspera_key_path is not None else '~/.aspera/connect/etc/asperaweb_id_dsa.openssh'

    def query(self, project_id: str, output_meta_file: str): ...

    def parse_meta(self, meta_file: str) -> tuple[list, list]: 
        '''Parse metadata file to extract SRR IDs and Sequence types.
        Args:
            meta_file: Path to metadata file (.csv).
        Returns:
            Tuple of list of SRR IDs and list of sequence types.

        Examples:
        >>> meta_file = './metadata.csv'
        >>> srr_ids, sequence_types = sratools.parse_meta(meta_file)
        '''

        meta_df = pd.read_csv(meta_file)

        return meta_df['Run'].tolist(), meta_df['LibraryLayout'].tolist()

    def aspera(self, srr_id: str, suffix: str, output_dir: str, method: str = 'ascp'):
        '''Download SRA data using Aspera connect according to SRR ID.
        Args:
            srr_id: SRA Run ID to download.
            suffix: Suffix for output file name. Options are '_1', '_2', or ''.
            output_dir: Directory to save downloaded files.
        Returns:
            None

        Examples:
        >>> srr_id = 'SRR12345678'
        >>> suffix = '_1'
        >>> output_dir = './data'
        >>> sratoools.aspera(srr_id, suffix, output_dir)
        '''

        def ena_path(srr):
            if len(srr) == 10:  # SRR + 7位
                return f"/vol1/fastq/{srr[:6]}/{srr}/{srr}{suffix}.fastq.gz"
            elif len(srr) == 11:  # SRR + 8位
                return f"/vol1/fastq/{srr[:6]}/0{srr[-2:]}/{srr}/{srr}{suffix}.fastq.gz"
            elif len(srr) == 12:  # SRR + 9位
                return f"/vol1/fastq/{srr[:6]}/{srr[-6:-3]}/{srr[-3:]}/{srr}/{srr}{suffix}.fastq.gz"

        ascp_cmd = [
            f'{self.aspera_path}',
            # ' --help',
            f' -QT -k 1 -l 300m -P 33001',
            f' -i {self.aspera_key_path}',
            f' era-fasp@fasp.sra.ebi.ac.uk:{ena_path(srr_id)}',
            f' {output_dir}',
        ]

        wget_cmd = [
            'wget',
            f' -O {Path(output_dir) / f"{srr_id}{suffix}.fastq.gz"}',
            f' ftp://ftp.sra.ebi.ac.uk/{ena_path(srr_id)}'
        ]

        if method == 'ascp':
            run_command(" ".join(ascp_cmd))
        elif method == 'wget':
            run_command(" ".join(wget_cmd))
        else:
            raise ValueError('Invalid method.')

        return None

    def download(
        self,
        srr_ids: list,
        sequence_types: list[Literal['single', 'paired']],
        output_dir: str = '.',
        method: str = 'ascp',
        backend: Literal['loky', 'multiprocessing', 'threading'] = 'threading',
        pre_dispatch: int = 5,
        verbose: int = 10,
    ):
        '''Download SRA data according to SRR IDs using specified method.
        Args:
            srr_ids: List of SRA IDs to download.
            sequence_types: List of sequence types to download. Options are 'paired' and 'single'.
            method: Method to use for downloading. Options are 'ascp' and 'wget'.
            output_dir: Directory to save downloaded files.
        Returns:
            None

        Examples:
        >>> srr_ids = ['SRR12345678', 'SRR12345679']
        >>> sequence_types = ['paired', 'single']
        >>> output_dir = './data'
        >>> sratoools = SRATools()
        >>> sratoools.download(srr_ids, sequence_types, output_dir, method='ascp', backend='threading', pre_dispath=5, verbose=10)
        '''

        from joblib import Parallel, delayed

        # Generate tasks for parallel execution
        tasks: list = []

        for srr_id, sequence_type in zip(srr_ids, sequence_types):
            if 'paired' in sequence_type.lower():
                paired_left_tasks = [delayed(self.aspera)(srr_id=srr_id, suffix='_1', output_dir=output_dir, method=method)]
                tasks.extend(paired_left_tasks)
                paired_right_tasks = [delayed(self.aspera)(srr_id=srr_id, suffix='_2', output_dir=output_dir, method=method)]
                tasks.extend(paired_right_tasks)

            elif 'single' in sequence_type.lower():
                single_tasks = [delayed(self.aspera)(srr_id=srr_id, suffix='', output_dir=output_dir, method=method)]
                tasks.extend(single_tasks)

            else:
                logger.warning(f'Unknown sequence type: {sequence_type} for SRR ID: {srr_id}. Skipping.')

        # Run tasks in parallel
        parallel = Parallel(n_jobs=-1, backend=backend, pre_dispatch=pre_dispatch, verbose=verbose)
        parallel(tasks)


        return None
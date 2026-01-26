
from typing import List, Tuple, Dict, Optional, Union, Any

from loguru import logger

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 


from skbio.stats import subsample_counts
from skbio.diversity import alpha_diversity


class AutoRarefaction:
    '''Rarefaction analysis for OTU/ASV tables. Includes subsampling to multiple depths, and visualization of rarefaction curves.
    Examples:
    >>> ar = AutoRarefaction(otu_table=otu_table_df)
    >>> rarefaction_data_df = ar.rarefaction_curves(max_depth=100000, steps=10, repeats=30, alpha_metric_name='Shannon')
    >>> otu_rarefied_table: pd.DataFrame = ar.rarefy(depths=[3947], repeats=1, alpha_metric_name='chao1', is_final_rarefaction=True)
    '''

    def __init__(self, otu_table: pd.DataFrame):
        self.otu_table = otu_table.astype(int) # OTU_table，row：OTU_ID, column：SampleID
        self.depth_per_sample = self.otu_table.sum(axis=0)  # 样本测序深度

        self.min_depth = self.depth_per_sample.min()
        logger.info(f'Minimum sample depth: {self.min_depth}')


    def _subsample_counts(self, counts: np.ndarray, depth: int) -> np.ndarray:
        '''Subsample counts to the given depth.'''

        total = counts.sum()
        if total < depth:
            return None # ignore samples with insufficient depth
        else:
            return subsample_counts(counts=counts, n=depth)
    
    def _rarefy_all(self, depths: List[int], repeats: int = 30, alpha_metric_name: Optional[str] = None) -> pd.DataFrame:
        '''Rarefy the OTU table to multiple depths.
        Args:
            depths : List[int], List of depths to rarefy to.
            repeats : int, Number of repeated subsamplings per depth.
            alpha_metric_name : str, Name of the alpha diversity metric to compute.
        Returns:
            Dict[int, pd.DataFrame], Mapping from depth to rarefied OTU table.
        Examples:
        >>> all_rarefaction_table_df = ar.rarefaction_curves(max_depth=100000, steps=10, repeats=30, alpha_metric_name='chao1')
        >>> otu_rarefied_table: pd.DataFrame = ar.rarefy(depths=[3947], repeats=1, alpha_metric_name='chao1', is_final_rarefaction=True)
        '''

        rarefied_tables = {
            'SampleID': [],
            'Depth': [],
            alpha_metric_name: [],
        }

        # per sample
        for sample in self.otu_table.columns:
            counts = self.otu_table[sample].values

            # per depth
            for depth in depths:

                alpha_metrics: List[float] = []
                # per repeat
                for _ in range(repeats):
                    subsampled = self._subsample_counts(counts, depth)

                    if subsampled is not None:
                        alpha_metric = alpha_diversity(
                            metric = alpha_metric_name,
                            counts = subsampled,
                            ids= [sample],
                        )
                        alpha_metrics.append(alpha_metric[sample])
                    else:
                        continue # ignore samples with insufficient depth (sum of counts < depth)
                if len(alpha_metrics) > 0:
                    avg_alpha_metric = np.mean(alpha_metrics)
                else:
                    avg_alpha_metric = np.nan

                rarefied_tables['SampleID'].append(sample)
                rarefied_tables['Depth'].append(depth)
                rarefied_tables[alpha_metric_name].append(avg_alpha_metric)

        return pd.DataFrame(rarefied_tables)
    
    def rarefaction_curves(self, max_depth: Optional[int] = None, steps: int = 20, repeats: int = 30, alpha_metric_name: str = "Shannon"):
        '''Plot rarefaction curves for each sample.
        Args:
            max_depth : Optional[int], Maximum subsampling depth. Default = min(sample_depth).
            steps : int, Number of depths.
            repeats : int, Number of repeated subsamplings per depth.
            alpha_metric_name: str, Name of the alpha diversity metric to compute.
        Returns:
            pd.DataFrame : DataFrame containing rarefaction data.
        '''

        # Determine max_depth
        if max_depth is None:
            max_depth = self.min_depth
            logger.info(f'Using min sample depth as max_depth: {max_depth}')

        depths = np.linspace(100, max_depth, steps)

        # Perform rarefaction
        rarefaction_data_df: pd.DataFrame = self._rarefy_all(
            depths= depths.tolist(),
            repeats= repeats,
            alpha_metric_name= alpha_metric_name,
        )

        # plot rarefaction curves
        sns.lineplot(
            data=rarefaction_data_df,
            x='Depth',
            y=alpha_metric_name,
            hue='SampleID',
            marker='o',
            alpha=0.6,
            palette='husl',
            legend=False, # disable legend for clarity, cause many samples.
        )

        return rarefaction_data_df
    
    def choose_best_depth(self, all_rarefied_table_df: pd.DataFrame, alpha_metric_name: str, coverage_threshold: float = 0.9) -> int:
        '''Choose the best rarefaction depth based on coverage threshold.
        Args:
            all_rarefied_table_df : pd.DataFrame, DataFrame containing rarefaction data.
            coverage_threshold : float, Coverage threshold to choose depth.
        Returns:
            int, Chosen depth.
        '''

        depth_counts_df = all_rarefied_table_df.groupby(by=['Depth'])[alpha_metric_name].mean().reset_index()
        depths = depth_counts_df['Depth'].values
        best_depth = None
        for d in depths:
            coverage = (self.depth_per_sample >= d).mean()
            if coverage >= coverage_threshold:
                best_depth = d
                break

        if best_depth is not None:
            logger.info(f'Chosen depth: {best_depth} with coverage >= {coverage_threshold}')
            return best_depth
        else:
            logger.warning(f'No depth found with coverage >= {coverage_threshold}, using max depth: {depths[-1]}')
            return depths[-1]
     
    def rarefy(self, depth: int) -> pd.DataFrame:
        '''Rarefy the OTU table to multiple depths.
        Args:
            depth : int, Depth to rarefy to.
        Returns:
            pd.DataFrame, Rarefied OTU table.
        Examples:
        >>> all_rarefaction_table_df = ar.rarefaction_curves(max_depth=100000, steps=10, repeats=30, alpha_metric_name='chao1')
        >>> otu_rarefied_table: pd.DataFrame = ar.rarefy(depth=3947)
        '''

        rarefied_tables_final: pd.DataFrame = pd.DataFrame(index=self.otu_table.index)

        # per sample
        for sample in self.otu_table.columns:
            counts = self.otu_table[sample].values

            subsampled = self._subsample_counts(counts, depth)

            if subsampled is not None:
                subsampled_df: pd.DataFrame = pd.DataFrame(data=subsampled, index=self.otu_table.index, columns=[sample])  # for the final rarefied table, when the depth is got.
                rarefied_tables_final = pd.concat([rarefied_tables_final, subsampled_df], axis=1)
            else:
                continue # ignore samples with insufficient depth (sum of counts < depth)
            
        return rarefied_tables_final
    


from loguru import logger
from typing import List, Tuple, Dict, Optional, Union

from pathlib import Path

import re
from tqdm import tqdm

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns 

from sklearn.preprocessing import StandardScaler, MinMaxScaler, LabelEncoder, OneHotEncoder


def sintax_to_taxonomy(sintax_file_path: str, fill_unclassified: bool = True):
    '''Transformer sintax format to taxonomy format.
    Args:
        sintax_file_path: str, path to sintax file.
        fill_unclassified: bool, whether to fill unclassified taxonomy with "Unclassified_<Rank>".
    
    Return taxonomy format:
        OTU_ID	Kingdom	Phylum	Class	Order	Family	Genus	Species
        OTU_1	Bacteria	Firmicutes	Bacilli	Lactobacillales	Lactobacillaceae	Lactobacillus	NA
        OTU_2	Bacteria	Proteobacteria	Gammaproteobacteria	...	...	...	...

    Examples:
    >>>taxonomy_df = sintax_to_taxonomy(sintax_file_path='./asv.sintax')
    '''

    otu_sintax_df = pd.read_csv(sintax_file_path, sep='\t', header=None, names=['OTU_ID', 'Taxonomy_with_confidence', 'Strand', 'Taxonomy'])

    RANK_MAP = {
        "d": "Domain",
        "k": "Kingdom",
        "p": "Phylum",
        "c": "Class",
        "o": "Order",
        "f": "Family",
        "g": "Genus",
        "s": "Species"
    }
    
    taxonomy_records: list = []

    for idx, row in otu_sintax_df.iterrows():
        otu_id = row['OTU_ID']
        tax_str = row["Taxonomy_with_confidence"]

        # init firstly taxonomy information with NA.
        tax_dict: dict = {rank: "NA" for rank in RANK_MAP.values()} #         Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species
        tax_dict['OTU_ID'] = otu_id                                 # OTU_ID, Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species

        if isinstance(tax_str, str) and tax_str.strip():
            items = tax_str.split(",") # e.g. d:Bacteria(1.0000),p:Actinomycetota(1.0000),
            for item in items:
                m = re.match(r"([dkpcfogs]):([^()]+)\([\d\.]+\)", item.strip())
                if m:
                    rank_short, tax_name = m.group(1), m.group(2)
                    full_rank = RANK_MAP.get(rank_short)
                    if full_rank:
                        tax_dict[full_rank] = tax_name

        taxonomy_records.append(tax_dict)

    taxonomy_df = pd.DataFrame(taxonomy_records)
    taxonomy_df = taxonomy_df[["OTU_ID"] + list(RANK_MAP.values())] # according rank order to order.

    if fill_unclassified:
        for rank in RANK_MAP.values():
            taxonomy_df[rank] = taxonomy_df[rank].fillna("NA")
            taxonomy_df[rank] = taxonomy_df[rank].replace("NA", f"Unclassified_{rank}")

    return taxonomy_df


class Amplicon:
    '''Operator of feature tables, e.g. OTU, Taxonomy, Metadata and so on.'''

    def __init__(self, 
        otu_file_path: str,
        sintax_file_path: str, 
        metadata_file_path: str,   
        otu_id_file_path: str = None
    ):
        self.otu_file_path = otu_file_path
        self.sintax_file_path = sintax_file_path
        self.metadata_file_path = metadata_file_path
        self.otu_id_file_path = otu_id_file_path

    def features_parser(self):
        '''Load features (e.g. otu_table, taxonomy_table, metadata_table)
        Return:
            if otu_id_file_path is not None:
                (index2id, id2index, otu_table, taxonomy_table, metadata_table)
            else:
                (otu_table, taxonomy_table, metadata_table)
        '''

        # load otu table
        otu_table = pd.read_csv(self.otu_file_path, sep='\t', index_col=0)
        otu_table.columns.name = "SampleID"

        # load taxonomy
        taxonomy_table = sintax_to_taxonomy(sintax_file_path=self.sintax_file_path)
        taxonomy_table.set_index('OTU_ID', inplace=True)

        # load metadata
        metadata_table = pd.read_excel(self.metadata_file_path, sheet_name='clean')
        metadata_table.set_index('SampleID', inplace=True)

        if self.otu_id_file_path is not None:
            # load otu id mapping
            otu_id_df = pd.read_csv(self.otu_id_file_path, sep='\t', header=None, names=['Index', 'OTU_ID'])
            index2id = dict(zip(otu_id_df['Index'], otu_id_df['OTU_ID']))
            id2index = dict(zip(otu_id_df['OTU_ID'], otu_id_df['Index']))
            return (index2id, id2index, otu_table, taxonomy_table, metadata_table)

        return (otu_table, taxonomy_table, metadata_table)

    @staticmethod
    def merge_otu_metadata_taxonomy_table(otu_table_df: pd.DataFrame, metadata_table_df: pd.DataFrame = None, taxonomy_table_df: pd.DataFrame = None):
        '''Merge otu_table and metadata_table or taxonomy_table.
        Args:
            otu_table_df: otu table with 'OUT_ID' as the index.
            metadata_table_df: metadata table with "SampleID' as the index.
            taxonomy_table_df: taxonomy table with "OTU_ID'' as th index.
        Returns:
            otu_metadata_table_df: otu table merged with metadata table.
            otu_taxonomy_table_df: otu table merged with taxonomy table.

        Examples:
        >>>otu_meta, otu_taxa = merge_otu_metadata_taxonomy_table(otu_table, metadata_table, taxonomy_table)
        '''
        
        # wide to long with columns: [OUT_ID, SampleID, Abundance]
        otu_long_df = otu_table_df.reset_index().melt(id_vars=['OTU_ID'], var_name='SampleID', value_name='Abundance')

        # merge otu_long_df with metadata_table_df with columns: [OUT_ID, SampleID, Abundance, Group, SoilType, ...]
        if metadata_table_df is not None:
            otu_metadata = pd.merge(left=otu_long_df, right=metadata_table_df, left_on='SampleID', right_index=True, how='left')

        # merge otu_long_df with taxonomy_table_df with columns: [OUT_ID, SampleID, Abundance, Domain, Kingdom, ...]
        if taxonomy_table_df is not None:
            otu_taxonomy = pd.merge(left=otu_long_df, right=taxonomy_table_df, left_on='OTU_ID', right_index=True, how='left')

        # return
        if (metadata_table_df is not None) and (taxonomy_table_df is not None):
            return (otu_metadata.set_index('OTU_ID'), otu_taxonomy.set_index('OTU_ID'))
        elif metadata_table_df is not None:
            return otu_metadata.set_index('OTU_ID')
        elif taxonomy_table_df is not None:
            return otu_taxonomy.set_index('OTU_ID')
        else:
            raise ValueError('ERROR: There are no metadata and taxonomy table.')

    @staticmethod
    def group_by_metadata_taxonomy_table(otu_taxonomy_table_df: pd.DataFrame, taxonomy_rank_col: list = ['Species'], metadata_table_df: pd.DataFrame = None):
        '''Group OTU table by taxonomy rank and metadata group.
        Args:
            otu_taxonomy_table_df: OTU table merged with taxonomy table with columns: [OTU_ID, SampleID, Abundance, Domain, Kingdom, ...]
            taxonomy_rank_col: list, taxonomy rank columns to group by, e.g. ['Phylum', 'Genus']
            metadata_table_df: metadata table with "SampleID' as the index, columns: [Group, SoilType, ...]
        Returns:
            otu_grouped_rank_group: DataFrame, grouped OTU table with index: taxonomy_rank_col, columns: [SampleID, Abundance, metadata_group_col...]
        
        Examples:
        >>>grouped_otu_df = group_by_metadata_taxonomy_table(otu_taxonomy_table_df, taxonomy_rank_col=["Phylum"], metadata_table_df=metadata_table)
        '''

        # 1. otu_taxonomy_table: groupby taxonomy rank.
        # index: taxonomy_rank_col; columns: SampleID, Abundance.
        otu_grouped_taxonomy_rank_table = otu_taxonomy_table_df.groupby(by=taxonomy_rank_col + ['SampleID'])['Abundance'].agg('sum').reset_index().set_index(taxonomy_rank_col)

        # 2. merge with metadata_table_df
        # index: SampleID; columns: SampleID, Abundance, metadata_group_col...
        otu_grouped_rank_group = pd.merge(left=otu_grouped_taxonomy_rank_table, right=metadata_table_df, left_on='SampleID', right_index=True, how='left')

        return otu_grouped_rank_group
    
    @staticmethod
    def subset_by_group(taxa_metadata_table_df: pd.DataFrame, group_col: str, group_name: str = None, 
                        save_dir_path: str = None, prefix: str = 'otu'):
        '''Extract DataFrame from taxa_metadata_table according to group_col == group_name.
        Args:
            taxa_metadata_table_df: DataFrame, taxa [OTU, Species, Genus, ...] merged with metadata [SampleID, Abundance, Group, SoilType, ...]
            group_col: str, metadata column to group by, e.g. 'PartName'`
            group_name: str, metadata group name to filter, e.g. 'Bulk'. If None, return all groups as a dict.
        Returns:
            if group_name is None:
                groups_dict: dict, key is group_name, value is DataFrame of the group.
            else:
                df_group: DataFrame, filtered DataFrame of the group.
        '''

        def long_to_wide(df_long: pd.DataFrame):
            df_wide = df_long.reset_index().pivot(
                index=[prefix],
                columns='SampleID',
                values='Abundance',
            )
            df_wide.index.rename('OTU_ID', inplace=True) # rename index, for subsequent FastSpar input requirement.
            return df_wide
        
        if group_name is None:
            groups = taxa_metadata_table_df[group_col].unique().tolist()
            groups_dict = {}
            for g in tqdm(groups, desc=f'Grouping by {group_col}'):
                df_group = taxa_metadata_table_df[taxa_metadata_table_df[group_col] == g]
                
                # save if needed
                if save_dir_path is not None:
                    save_path = Path(save_dir_path) / f"{prefix}_{group_col}_{g}.txt"
                    if not save_path.parent.exists():
                        save_path.parent.mkdir(parents=True, exist_ok=True)
                    long_to_wide(df_group).to_csv(save_path, sep='\t')

                groups_dict[g] = df_group
            return groups_dict
        
        else:
            df_group = taxa_metadata_table_df[taxa_metadata_table_df[group_col] == group_name]
            if save_dir_path is not None:
                save_path = Path(save_dir_path) / f"{prefix}_{group_col}_{group_name}.txt"
                if not save_path.parent.exists():
                    save_path.parent.mkdir(parents=True, exist_ok=True)
                long_to_wide(df_group).to_csv(save_path, sep='\t')
            return df_group


class OTUQC:
    '''OTU table quality control including filter and sparsity curve.

    Examples:
    >>>otu_table_filetered = OTUQC.filter_otu_table(otu_table_df=otu_table, min_prevalence=0.01, min_abundance=0.00001)
    >>>otu_table_filtered_topk = OTUQC.keep_top_otus(otu_table_df=otu_table_filtered, top_k=3000)
    >>>otu_table_for_network = OTUQC.otu_for_network(otu_table_df=otu_table, min_prevalence=0.01, min_abundance=0.00001, top_k=3000)
    >>>sparsities = OTUQC.sparsity_curve(otu_table_df=otu_table_filtered_topk, steps=[1000, 2000, 5000, 10000, 20000])
    '''

    @staticmethod
    def filter_otu_table_by_sample(otu_table_df: pd.DataFrame, method: str = 'iqr', min_max_threshold: Tuple[float, float] = (0.25, 0.75), iqr_scale_factor: float = 1.5, show: bool = True) -> pd.DataFrame:
        '''Filter samples based on total abundance.
        Args:
            otu_table_df: OTU table DataFrame, OTU_ID as index, sample IDs as columns.
            method: str, filtering method, options: 'iqr', 'quantile', 'abundance'; default is 'iqr'.
            min_max_threshold: Tuple[float, float], minimum and maximum quantile to filter samples based on total abundance.
            iqr_scale_factor: float, scale factor for IQR method.
            show: bool, whether to show the filtering plots.
        Returns:
            filtered_otu_table_df: filtered OTU table DataFrame.
        '''

        total_abundance = otu_table_df.sum(axis=0) # for each sample
        min_quantile, max_quantile = min_max_threshold

        if method == 'iqr':
            # IQR method [default]
            q1 = total_abundance.quantile(q=min_quantile)
            q3 = total_abundance.quantile(q=max_quantile)
            iqr = q3 - q1
            min_threshold = q1 - iqr_scale_factor * iqr
            max_threshold = q3 + iqr_scale_factor * iqr
        elif method == 'quantile':
            min_threshold = total_abundance.quantile(q=min_quantile)
            max_threshold = total_abundance.quantile(q=max_quantile)
        elif method == 'abundance':
            min_threshold = min_quantile
            max_threshold = max_quantile
        else:
            raise ValueError(f"Filtering method '{method}' is not supported. Choose from 'iqr' or 'quantile'.")

        logger.info(f"Sample total abundance filtering thresholds: min={min_threshold}, max={max_threshold}")

        keep = (total_abundance >= min_threshold) & (total_abundance <= max_threshold)
        filtered_otu_table_df = otu_table_df.loc[:, keep]

        logger.info(f"Raw/Filtered shape: {otu_table_df.shape} -> {filtered_otu_table_df.shape}, kept samples: {keep.sum() * 100/len(keep):.2f}%")

        if show:
            plt.subplot(221)
            sns.boxplot(data=total_abundance, orient='h')
            plt.title('Total Abundance per Sample')

            plt.subplot(222)
            sns.histplot(total_abundance, kde=True)
            plt.axvline(min_threshold, color='green', linestyle='--', label='Min Threshold', alpha=0.7)
            plt.axvline(max_threshold, color='blue', linestyle='--', label='Max Threshold', alpha=0.7)
            plt.legend()
            plt.title('Total Abundance Distribution')

            plt.subplot(223)
            sns.boxplot(data=filtered_otu_table_df.sum(axis=0), orient='h')
            plt.title('Filtered Total Abundance per Sample')

            plt.subplot(224)
            sns.histplot(filtered_otu_table_df.sum(axis=0), kde=True)
            plt.axvline(min_threshold, color='green', linestyle='--', label='Min Threshold', alpha=0.7)
            plt.axvline(max_threshold, color='blue', linestyle='--', label='Max Threshold', alpha=0.7)
            plt.legend()
            plt.title('Filtered Total Abundance Distribution')

            plt.tight_layout()
            plt.show()
        
        return filtered_otu_table_df

    @staticmethod
    def filter_otu_table_by_otu(otu_table_df: pd.DataFrame, min_prevalence: float = 0.01, min_abundance: float = 0.00001, show: bool = False) -> pd.DataFrame:
        '''Filter OTU table based on prevalence and abundance.
        Meanwhile, [core microbiome] can be extracted via setting min_prevalence=0.5, min_abundance=None.
        Args:
            otu_table_df: OTU table DataFrame, OTU_ID as index, sample IDs as columns.
            min_prevalence: float, minimum prevalence threshold (0-1), default is 0.01.
            min_abundance: float, minimum abundance threshold (0-1), default is 0.00001.
            show: bool, whether to show the filtering plots.
        Returns:
            filtered_otu_table_df: filtered OTU table DataFrame.
        '''

        n_samples = otu_table_df.shape[1]

        # Calculate prevalence threshold
        prevalence = (otu_table_df > 0).sum(axis=1) / n_samples
        keep_prevalence = prevalence >= min_prevalence

        # Calculate total abundance threshold
        if min_abundance is not None:
            rel_abundance = otu_table_df.div(otu_table_df.sum(axis=0), axis=1)
            global_abundance = rel_abundance.mean(axis=1)           # or sum(axis=1), mean operation is more strict.
            keep_abundance = global_abundance >= min_abundance
        else:
            keep_abundance = pd.Series(True, index=otu_table_df.index)

        keep = keep_prevalence & keep_abundance

        filtered_otu_table_df = otu_table_df.loc[keep, :]

        logger.info(f"Raw/Filtered shape: {otu_table_df.shape} -> {filtered_otu_table_df.shape}, kept samples: {filtered_otu_table_df.shape[0] * 100/otu_table_df.shape[0]:.2f}%")

        if show:
            plt.subplot(221)
            sns.histplot(prevalence, kde=True)
            plt.axvline(min_prevalence, color='green', linestyle='--', label='Min Prevalence', alpha=0.7)
            plt.legend(loc='upper right')
            plt.xscale('log')
            plt.title('OTU Prevalence Distribution')

            plt.subplot(222)
            sns.histplot(prevalence.loc[filtered_otu_table_df.index], kde=True)
            plt.xscale('log')
            plt.title('Filtered OTU Prevalence Distribution')

            if min_abundance is not None:
                plt.subplot(223)
                sns.histplot(global_abundance, kde=True)
                plt.axvline(min_abundance, color='green', linestyle='--', label='Min Abundance', alpha=0.7)
                plt.legend(loc='upper right')
                plt.xscale('log')
                plt.title('OTU Global Abundance Distribution')

                plt.subplot(224)
                sns.histplot(rel_abundance.loc[filtered_otu_table_df.index, :].mean(axis=1), kde=True)
                plt.xscale('log')
                plt.title('Filtered OTU Global Abundance Distribution')

            plt.tight_layout()
            plt.show()

        return filtered_otu_table_df 
    
    @staticmethod
    def keep_top_otus(otu_table_df: pd.DataFrame, top_k: int = 3000):
        '''Keep top_k OTUs based on total abundance.
        Args:
            otu_table_df: OTU table DataFrame, OTU_ID as index, sample IDs as columns.
            top_k: int, number of top OTUs to keep based on total abundance.
        Returns:
            filtered_otu_table_df: filtered OTU table DataFrame containing top_k OTUs.
        '''

        total_abundance = otu_table_df.sum(axis=1)
        top = total_abundance.sort_values(ascending=False).head(top_k).index

        return otu_table_df.loc[top, :]
    
    @staticmethod
    def otu_for_network(otu_table_df: pd.DataFrame, min_prevalence: float = 0.01, min_abundance: float = 0.00001, top_k: int = 3000):
        '''Comprehensive filter function including prevalence, abundance and top_k.'''

        filtered_otu_table = OTUQC.filter_otu_table_by_otu(
            otu_table_df=otu_table_df,
            min_prevalence=min_prevalence,
            min_abundance=min_abundance,
            show=False,
        )

        filtered_otu_table = OTUQC.keep_top_otus(
            otu_table_df=filtered_otu_table,
            top_k=top_k
        )

        return filtered_otu_table

    @staticmethod
    def sparsity_curve(otu_table_df: pd.DataFrame, steps=[1000, 2000, 5000, 10000, 20000]):
        '''Plot sparsity curve based on different number of OTUs.'''

        sparsities = []

        for k in steps:
            if k >= otu_table_df.shape[0]: # number of OTUs
                sparsities.append(np.nan)
                continue 

            subset = otu_table_df.iloc[:k, :]
            # Spearman correlation justfor sparsity test
            corr = subset.T.corr(method='spearman')
            sparsity = (corr == 0).mean().mean()
            sparsities.append(sparsity)

            print(f"{k} OTUs -> Sparsity: {sparsity:.4f}")

        plt.plot(steps[:len(sparsities)], sparsities, marker='o')
        plt.xlabel('OTU numbers')
        plt.ylabel('Sparsity')
        plt.title('Sparsity Curve')
        plt.grid(True)

        return sparsities
    
    @staticmethod
    def normalize(otu_table_df: pd.DataFrame, method: str = 'rel') -> pd.DataFrame:
        '''Normalize OTU table with specified method.
        Args:
            otu_table_df: OTU table DataFrame, OTU_ID as index, sample IDs as columns.
            method: str, normalization method, options: 'tss', 'clr', 'rarefy', 'none'.
        Returns:
            normalized_otu_table_df: normalized OTU table DataFrame.
        Examples:
            >>> OTUQC.normalize(otu_table_df, method='tss')
            >>> OTUQC.normalize(otu_table_df, method='clr')
        '''

        if method == 'none':
            normalized_otu_table_df = otu_table_df.copy()

        elif method == 'rel':
            # normalized_otu_table_df = otu_table_df.div(otu_table_df.sum(axis=0), axis=1)
            sums = otu_table_df.sum(axis=0).replace(0, np.nan) # avoid division by zero
            normalized_otu_table_df = otu_table_df.div(sums, axis=1).fillna(0)

        elif method == 'clr':
            rel_abundance = otu_table_df.div(otu_table_df.sum(axis=0), axis=1)
            log_rel_abundance = np.log(rel_abundance.replace(0, np.nan))
            gm = log_rel_abundance.mean(axis=0)
            normalized_otu_table_df = log_rel_abundance.subtract(gm, axis=1).fillna(0)

        elif method == 'rarefy':
            def rarefy_column(col: pd.Series, depth: int):
                if col.sum() < depth:
                    raise ValueError(f"Cannot rarefy sample with total count {col.sum()} to depth {depth}.")
                probabilities = col / col.sum()
                rarefied_counts = np.random.multinomial(depth, probabilities)
                return pd.Series(rarefied_counts, index=col.index)

            min_depth = otu_table_df.sum(axis=0).min()
            normalized_otu_table_df = otu_table_df.apply(lambda col: rarefy_column(col, int(min_depth)), axis=0)

        elif method == 'deseq2':
            # Placeholder for DESeq2 normalization
            # In practice, this would require calling R's DESeq2 package via rpy2 or similar.
            raise NotImplementedError("DESeq2 normalization requires R's DESeq2 package and is not implemented in this function.")
        
        elif method == 'css':
            # Placeholder for CSS normalization
            # In practice, this would require a specific implementation or package.
            raise NotImplementedError("CSS normalization is not implemented in this function.")
        
        elif method == 'tmm':
            # Placeholder for TMM normalization
            # In practice, this would require a specific implementation or package.
            raise NotImplementedError("TMM normalization is not implemented in this function.")
        
        elif method == 'tpm':
            # Placeholder for TPM normalization
            # In practice, this would require gene length information.
            raise NotImplementedError("TPM normalization requires gene length information and is not implemented in this function.")    

        else:
            raise ValueError(f"Normalization method '{method}' is not supported. Choose from 'tss', 'clr', 'rarefy', or 'none', 'deseq2', 'css', 'tmm', 'tpm'.")

        return normalized_otu_table_df
    
    def _choose_rarefaction_depth(self, strategy: str, q: float ): ...
    
    @staticmethod
    def rarefy(otu_table_df: pd.DataFrame, depth_method: Union[str, int] = 'min') -> pd.DataFrame:
        '''Rarefy OTU table to the depth.
        Args:
            otu_table_df: OTU table DataFrame, OTU_ID as index, sample IDs as columns.
        Returns:
            rarefied_otu_table_df: rarefied OTU table DataFrame.
        '''

        if isinstance(depth_method, int):
            depth = depth_method

        elif isinstance(depth_method, str) and depth_method == 'min':
            depth = int(otu_table_df.sum(axis=0).min())

        elif isinstance(depth_method, str) and depth_method.endswith('%'):
            perc = float(depth_method.strip('%')) / 100.0
            depth = int(otu_table_df.sum(axis=0).min() * perc)


class Metadata:
    '''Operator of metadata table.'''

    @staticmethod
    def metadata_continuous_normalize(metadata_table_df: pd.DataFrame, continuous_cols: List[str], method: str = 'standard') -> pd.DataFrame:
        '''Normalize continuous metadata columns.
        Args:
            metadata_table_df: metadata table DataFrame, SampleID as index.
            continuous_cols: list of str, continuous columns to normalize.
            method: str, normalization method, options: 'standard', 'minmax', 'log'.
        Returns:
            normalized_metadata_table_df: normalized metadata table DataFrame.
        '''

        normalized_metadata_table_df = metadata_table_df.copy()[continuous_cols]

        if method == 'standard':
            scaler = StandardScaler()
            normalized_metadata_table_df[continuous_cols] = scaler.fit_transform(normalized_metadata_table_df[continuous_cols])

        elif method == 'minmax':
            scaler = MinMaxScaler()
            normalized_metadata_table_df[continuous_cols] = scaler.fit_transform(normalized_metadata_table_df[continuous_cols])
        
        elif method == 'log':
            normalized_metadata_table_df[continuous_cols] = np.log1p(normalized_metadata_table_df[continuous_cols])

        else:
            raise ValueError(f"Normalization method '{method}' is not supported. Choose from 'standard', 'minmax', or 'log'.")

        return normalized_metadata_table_df

    @staticmethod
    def metadata_categorical_encode(metadata_table_df: pd.DataFrame, categorical_cols: List[str], method: str = 'onehot') -> pd.DataFrame:
        '''Encode categorical metadata columns.
        Args:
            metadata_table_df: metadata table DataFrame, SampleID as index.
            categorical_cols: list of str, categorical columns to encode.
            method: str, encoding method, options: 'label', 'onehot'.
        Returns:
            encoded_metadata_table_df: encoded metadata table DataFrame.
        '''

        encoded_metadata_table_df = metadata_table_df.copy()[categorical_cols]

        if method == 'onehot':
            encoder = OneHotEncoder(sparse_output=False, drop='first')
            onehot_encoded_array = encoder.fit_transform(encoded_metadata_table_df[categorical_cols])
            onehot_encoded_df = pd.DataFrame(
                onehot_encoded_array,
                index=encoded_metadata_table_df.index,
                columns=encoder.get_feature_names_out(categorical_cols)
            )
            return onehot_encoded_df
        
        elif method == 'label':
            for col in categorical_cols:
                le = LabelEncoder()
                encoded_metadata_table_df[col] = le.fit_transform(encoded_metadata_table_df[col])
            return encoded_metadata_table_df            

        else:
            raise ValueError(f"Encoding method '{method}' is not supported. Choose from 'label' or 'onehot'.")



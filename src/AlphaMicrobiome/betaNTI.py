# import numpy as np
# import pandas as pd
# import cupy as cp
# from skbio import TreeNode
# from tqdm import tqdm
# from scipy.spatial.distance import braycurtis


# # ============================================================
# #  BetaNTI (GPU)
# # ============================================================

# class BetaNTI_GPU:
#     """
#     GPU-accelerated betaNTI calculator (Stegen framework compatible)
#     """

#     def __init__(
#         self,
#         otu_table: pd.DataFrame,
#         tree_file: str,
#         metadata: pd.DataFrame,
#         n_null: int = 999,
#         random_state: int = 42,
#     ):
#         # ------------------------------
#         # Basic data
#         # ------------------------------
#         self.otu = otu_table.copy()
#         self.samples = self.otu.columns.tolist()
#         self.otus = self.otu.index.tolist()
#         self.metadata = metadata.loc[self.samples]

#         self.n_null = n_null
#         self.random_state = random_state
#         np.random.seed(random_state)

#         # ------------------------------
#         # Load tree
#         # ------------------------------
#         self.tree = TreeNode.read(tree_file)
#         self.tip_map = {tip.name: tip for tip in self.tree.tips()}
#         self._check_tree_match()

#         # ------------------------------
#         # Phylogenetic distance
#         # ------------------------------
#         self.phylo_dist = self._build_phylo_distance_matrix()
#         self.phylo_dist_gpu = cp.asarray(self.phylo_dist)

#         # OTU matrix to GPU
#         self.otu_gpu = cp.asarray(self.otu.values.astype(np.float32))

#     # --------------------------------------------------
#     # Internal
#     # --------------------------------------------------

#     def _check_tree_match(self):
#         missing = set(self.otus) - set(self.tip_map.keys())
#         if missing:
#             raise ValueError(
#                 f"{len(missing)} OTUs not found in tree tips. "
#                 f"Example: {list(missing)[:5]}"
#             )

#     def _build_phylo_distance_matrix(self):
#         n = len(self.otus)
#         dist = np.zeros((n, n), dtype=np.float32)

#         for i, otu_i in tqdm(
#             enumerate(self.otus),
#             total=n,
#             desc="Phylogenetic distance",
#         ):
#             node_i = self.tip_map[otu_i]
#             for j in range(i + 1, n):
#                 node_j = self.tip_map[self.otus[j]]
#                 d = node_i.distance(node_j)
#                 dist[i, j] = dist[j, i] = d

#         return dist

#     @staticmethod
#     def _beta_mntd_weighted(sample1, sample2, dist_mat):
#         idx1 = cp.where(sample1 > 0)[0]
#         idx2 = cp.where(sample2 > 0)[0]

#         if len(idx1) == 0 or len(idx2) == 0:
#             return cp.nan

#         w1 = sample1[idx1] / cp.sum(sample1[idx1])
#         w2 = sample2[idx2] / cp.sum(sample2[idx2])

#         d12 = dist_mat[idx1][:, idx2]

#         min1 = cp.min(d12, axis=1)
#         min2 = cp.min(d12, axis=0)

#         return cp.sum(w1 * min1) + cp.sum(w2 * min2)

#     # --------------------------------------------------
#     # betaMNTD / betaNTI
#     # --------------------------------------------------

#     def compute_beta_mntd_matrix(self):
#         n = len(self.samples)
#         beta_mntd = np.zeros((n, n), dtype=np.float32)

#         for i in tqdm(range(n), desc="Observed betaMNTD"):
#             for j in range(i + 1, n):
#                 d = self._beta_mntd_weighted(
#                     self.otu_gpu[:, i],
#                     self.otu_gpu[:, j],
#                     self.phylo_dist_gpu,
#                 )
#                 beta_mntd[i, j] = beta_mntd[j, i] = float(d.get())

#         return pd.DataFrame(beta_mntd, index=self.samples, columns=self.samples)

#     def compute_null_beta_mntd(self):
#         n = len(self.samples)
#         null_dist = np.zeros((self.n_null, n, n), dtype=np.float32)

#         for k in tqdm(range(self.n_null), desc="Null betaMNTD"):
#             shuffled = np.random.permutation(len(self.otus))
#             dist_null = self.phylo_dist_gpu[shuffled][:, shuffled]

#             for i in range(n):
#                 for j in range(i + 1, n):
#                     d = self._beta_mntd_weighted(
#                         self.otu_gpu[:, i],
#                         self.otu_gpu[:, j],
#                         dist_null,
#                     )
#                     null_dist[k, i, j] = null_dist[k, j, i] = float(d.get())

#         return null_dist

#     def compute_beta_nti(self):
#         obs = self.compute_beta_mntd_matrix().values
#         null = self.compute_null_beta_mntd()

#         null_mean = null.mean(axis=0)
#         null_sd = null.std(axis=0)
#         null_sd[null_sd == 0] = np.nan

#         beta_nti = (obs - null_mean) / null_sd
#         return pd.DataFrame(beta_nti, index=self.samples, columns=self.samples)


# # ============================================================
# #  Stegen Framework (betaNTI + RCbray)
# # ============================================================

# class Stegen(BetaNTI_GPU):
#     """
#     Full Stegen framework:
#     betaNTI + RCbray → ecological process inference
#     """

#     def __init__(
#         self,
#         otu_table,
#         tree_file,
#         metadata,
#         n_null=999,
#         rc_null=999,
#         random_state=42,
#     ):
#         super().__init__(
#             otu_table=otu_table,
#             tree_file=tree_file,
#             metadata=metadata,
#             n_null=n_null,
#             random_state=random_state,
#         )
#         self.rc_null = rc_null
#         self.random_state = random_state

#     # --------------------------------------------------
#     # RCbray
#     # --------------------------------------------------

#     def compute_rcbray(self):
#         otu = self.otu.loc[:, self.samples]
#         rng = np.random.default_rng(self.random_state)
#         n = len(self.samples)

#         obs = np.zeros((n, n))
#         for i in range(n):
#             for j in range(i + 1, n):
#                 obs[i, j] = obs[j, i] = braycurtis(
#                     otu.iloc[:, i], otu.iloc[:, j]
#                 )

#         null = np.zeros((self.rc_null, n, n))
#         for k in tqdm(range(self.rc_null), desc="RCbray null"):
#             shuffled = otu.apply(lambda x: rng.permutation(x.values), axis=0)
#             for i in range(n):
#                 for j in range(i + 1, n):
#                     null[k, i, j] = null[k, j, i] = braycurtis(
#                         shuffled.iloc[:, i], shuffled.iloc[:, j]
#                     )

#         rc = np.zeros((n, n))
#         for i in range(n):
#             for j in range(i + 1, n):
#                 p = (null[:, i, j] < obs[i, j]).mean()
#                 rc[i, j] = rc[j, i] = 2 * (p - 0.5)

#         return pd.DataFrame(rc, index=self.samples, columns=self.samples)

#     # --------------------------------------------------
#     # ✅ FIXED Stegen process (NO pandas ambiguity)
#     # --------------------------------------------------

#     @staticmethod
#     def stegen_process(beta_nti_df, rcbray_df):
#         bnti = beta_nti_df.values
#         rc = rcbray_df.values
#         n = bnti.shape[0]

#         processes = []

#         for i in range(n):
#             for j in range(i + 1, n):
#                 b = bnti[i, j]
#                 r = rc[i, j]

#                 if np.isnan(b) or np.isnan(r):
#                     continue

#                 if b > 2:
#                     proc = "Variable selection"
#                 elif b < -2:
#                     proc = "Homogeneous selection"
#                 else:
#                     if r > 0.95:
#                         proc = "Dispersal limitation"
#                     elif r < -0.95:
#                         proc = "Homogenizing dispersal"
#                     else:
#                         proc = "Drift"

#                 processes.append(proc)

#         return pd.Series(processes)

#     # --------------------------------------------------
#     # Group-wise summary (Result 3 / Figure 3)
#     # --------------------------------------------------

#     def summarize_by_group(self, group_col="PartName"):
#         beta_nti = self.compute_beta_nti()
#         rcbray = self.compute_rcbray()

#         results = {}

#         for g in self.metadata[group_col].unique():
#             samples = self.metadata.index[
#                 self.metadata[group_col] == g
#             ]

#             # 用 iloc 保证顺序一致
#             idx = [self.samples.index(s) for s in samples]

#             bnti_sub = beta_nti.iloc[idx, idx]
#             rc_sub = rcbray.iloc[idx, idx]

#             proc = self.stegen_process(bnti_sub, rc_sub)
#             results[g] = proc.value_counts(normalize=True)

#         return pd.DataFrame(results).fillna(0)



# import numpy as np
# import pandas as pd
# import cupy as cp
# from skbio import TreeNode
# from scipy.spatial.distance import braycurtis
# from tqdm import tqdm


# class Stegen:
#     """
#     Stegen framework:
#     betaNTI (phylogenetic) + RCbray (taxonomic)
#     """

#     def __init__(
#         self,
#         otu_table: pd.DataFrame,      # OTU x Sample (relative abundance)
#         tree_file: str,
#         metadata: pd.DataFrame,       # Sample x metadata
#         group_col: str = "PartName",
#         n_null: int = 999,
#         rc_null: int = 999,
#         random_state: int = 42,
#     ):
#         # ------------------------
#         # Basic
#         # ------------------------
#         self.otu = otu_table.copy()
#         self.samples = otu_table.columns.tolist()
#         self.otus = otu_table.index.tolist()
#         self.meta = metadata.loc[self.samples]
#         self.group_col = group_col

#         self.n_null = n_null
#         self.rc_null = rc_null
#         self.random_state = random_state
#         np.random.seed(random_state)

#         # ------------------------
#         # Load tree
#         # ------------------------
#         self.tree = TreeNode.read(tree_file)
#         self.tip_map = {tip.name: tip for tip in self.tree.tips()}
#         self._check_tree_match()

#         # ------------------------
#         # Cophenetic distance (CPU → GPU)
#         # ------------------------
#         self.phylo_dist = self._build_cophenetic_matrix()
#         self.phylo_dist_gpu = cp.asarray(self.phylo_dist)

#         # OTU table → GPU
#         self.otu_gpu = cp.asarray(self.otu.values)

#     # ==================================================
#     # Internal
#     # ==================================================
#     def _check_tree_match(self):
#         missing = set(self.otus) - set(self.tip_map)
#         if missing:
#             raise ValueError(
#                 f"{len(missing)} OTUs missing in tree tips, e.g. {list(missing)[:5]}"
#             )

#     def _build_cophenetic_matrix(self):
#         """
#         Fast & stable cophenetic distance matrix
#         Compatible with skbio >=0.5
#         """
#         n = len(self.otus)
#         dist = np.zeros((n, n), dtype=np.float32)

#         # depth from root
#         depth = {}
#         for name, node in self.tip_map.items():
#             d = 0.0
#             cur = node
#             while cur.parent is not None:
#                 d += cur.length or 0.0
#                 cur = cur.parent
#             depth[name] = d

#         for i in tqdm(range(n), desc="Cophenetic distance"):
#             ni = self.tip_map[self.otus[i]]
#             for j in range(i + 1, n):
#                 nj = self.tip_map[self.otus[j]]

#                 # ✅ 正确的 LCA 调用方式
#                 lca = self.tree.lowest_common_ancestor([ni, nj])

#                 d_lca = 0.0
#                 cur = lca
#                 while cur.parent is not None:
#                     d_lca += cur.length or 0.0
#                     cur = cur.parent

#                 dij = depth[self.otus[i]] + depth[self.otus[j]] - 2 * d_lca
#                 dist[i, j] = dist[j, i] = dij

#         return dist


#     @staticmethod
#     def _beta_mntd(sample1, sample2, dist_mat):
#         idx1 = cp.where(sample1 > 0)[0]
#         idx2 = cp.where(sample2 > 0)[0]

#         if idx1.size == 0 or idx2.size == 0:
#             return cp.nan  # 不用 cp.asarray

#         w1 = sample1[idx1] / cp.sum(sample1[idx1])
#         w2 = sample2[idx2] / cp.sum(sample2[idx2])

#         d12 = dist_mat[idx1][:, idx2]
#         min1 = cp.min(d12, axis=1)
#         min2 = cp.min(d12, axis=0)

#         return cp.sum(w1 * min1) + cp.sum(w2 * min2)


#     # ==================================================
#     # betaNTI
#     # ==================================================
#     def compute_beta_nti(self):
#         n = len(self.samples)

#         # observed betaMNTD
#         obs = np.zeros((n, n), dtype=np.float32)
#         for i in tqdm(range(n), desc="Observed betaMNTD"):
#             for j in range(i + 1, n):
#                 d = self._beta_mntd(
#                     self.otu_gpu[:, i],
#                     self.otu_gpu[:, j],
#                     self.phylo_dist_gpu,
#                 )
#                 # obs[i, j] = obs[j, i] = float(d.get())
#                 obs[i, j] = obs[j, i] = float(cp.asnumpy(d))

#         # null models
#         null = np.zeros((self.n_null, n, n), dtype=np.float32)
#         for k in tqdm(range(self.n_null), desc="Null betaMNTD"):
#             shuffled = np.random.permutation(len(self.otus))
#             dist_null = self.phylo_dist_gpu[shuffled][:, shuffled]

#             for i in range(n):
#                 for j in range(i + 1, n):
#                     d = self._beta_mntd(
#                         self.otu_gpu[:, i],
#                         self.otu_gpu[:, j],
#                         dist_null,
#                     )
#                     # null[k, i, j] = null[k, j, i] = float(d.get())
#                     null[k, i, j] = null[k, j, i] = float(cp.asnumpy(d))

#         beta_nti = (obs - null.mean(axis=0)) / null.std(axis=0)

#         return pd.DataFrame(beta_nti, index=self.samples, columns=self.samples)

#     # ==================================================
#     # RCbray
#     # ==================================================
#     def compute_rcbray(self):
#         n = len(self.samples)
#         otu = self.otu
#         rng = np.random.default_rng(self.random_state)

#         obs = np.zeros((n, n))
#         for i in range(n):
#             for j in range(i + 1, n):
#                 obs[i, j] = obs[j, i] = braycurtis(
#                     otu.iloc[:, i], otu.iloc[:, j]
#                 )

#         null = np.zeros((self.rc_null, n, n))
#         for k in tqdm(range(self.rc_null), desc="RCbray null"):
#             shuffled = otu.apply(lambda x: rng.permutation(x), axis=0)
#             for i in range(n):
#                 for j in range(i + 1, n):
#                     null[k, i, j] = null[k, j, i] = braycurtis(
#                         shuffled.iloc[:, i], shuffled.iloc[:, j]
#                     )

#         rc = np.zeros((n, n))
#         for i in range(n):
#             for j in range(i + 1, n):
#                 p = (null[:, i, j] < obs[i, j]).mean()
#                 rc[i, j] = rc[j, i] = 2 * (p - 0.5)

#         return pd.DataFrame(rc, index=self.samples, columns=self.samples)

#     # ==================================================
#     # Stegen process classification
#     # ==================================================
#     @staticmethod
#     def classify_pair(beta_nti, rc):
#         if beta_nti > 2:
#             return "Variable selection"
#         elif beta_nti < -2:
#             return "Homogeneous selection"
#         elif rc > 0.95:
#             return "Dispersal limitation"
#         elif rc < -0.95:
#             return "Homogenizing dispersal"
#         else:
#             return "Drift"

#     def summarize_by_group(self):
#         beta_nti = self.compute_beta_nti()
#         rcbray = self.compute_rcbray()

#         results = {}

#         for g in self.meta[self.group_col].unique():
#             samples = self.meta.index[self.meta[self.group_col] == g]

#             pairs = []
#             for i in range(len(samples)):
#                 for j in range(i + 1, len(samples)):
#                     s1, s2 = samples[i], samples[j]
#                     b = beta_nti.loc[s1, s2]
#                     r = rcbray.loc[s1, s2]
#                     if not np.isnan(b) and not np.isnan(r):
#                         pairs.append(self.classify_pair(b, r))

#             results[g] = pd.Series(pairs).value_counts(normalize=True)

#         return pd.DataFrame(results).fillna(0)


import numpy as np
import pandas as pd
import cupy as cp
# import numpy as cp
from skbio import TreeNode
from scipy.spatial.distance import braycurtis
from tqdm import tqdm
from pathlib import Path


class StegenPipeline:
    """
    Final Stegen framework:
    betaNTI (GPU accelerated) + RCbray + ecological process inference
    """

    def __init__(
        self,
        otu_table: pd.DataFrame,
        tree_file: str,
        metadata: pd.DataFrame,
        group_col: str,
        n_null: int = 999,
        rc_null: int = 999,
        random_state: int = 42,
        temp_dir: str = './tmp',
        device: int = 0,
    ):
        """
        otu_table: OTUs x Samples (RAW COUNTS recommended)
        metadata: index = sample IDs
        group_col: grouping column in metadata
        """

        self.otu = otu_table.copy()
        self.samples = otu_table.columns.tolist()
        self.otus = otu_table.index.tolist()

        self.meta = metadata.loc[self.samples]
        self.group_col = group_col

        self.n_null = n_null
        self.rc_null = rc_null
        self.random_state = random_state
        np.random.seed(random_state)
        self.temp_dir = Path(temp_dir) if temp_dir else None
        print(f"Temporary directory: {self.temp_dir}")
        self.device = device
        cp.cuda.Device(self.device).use()
        # ------------------------
        # Load tree
        # ------------------------
        self.tree = TreeNode.read(tree_file)
        self.tip_map = {tip.name: tip for tip in self.tree.tips()}
        self._check_tree_match()

        # ------------------------
        # Cophenetic distance (CPU → GPU)
        # ------------------------
        # self.phylo_dist = self._build_cophenetic_matrix()
        # self.phylo_dist_gpu = cp.asarray(self.phylo_dist)
        if self.temp_dir: 
            self.temp_dir.mkdir(parents=True, exist_ok=True)
            cache_file = self.temp_dir / "cophenetic_matrix.npy"
            if cache_file.exists():
                print(f'Loading cached cophenetic matrix from {cache_file}')
                self.phylo_dist = np.load(cache_file)
            else:
                self.phylo_dist = self._build_cophenetic_matrix()
                print(f'Saving cophenetic matrix to {cache_file}')
                np.save(cache_file, self.phylo_dist)
        self.phylo_dist_gpu = cp.asarray(self.phylo_dist)

        # ------------------------
        # OTU table → GPU
        # ------------------------
        self.otu_gpu = cp.asarray(self.otu.values)

        # cache nonzero OTUs per sample
        self.sample_nonzero = [
            cp.where(self.otu_gpu[:, i] > 0)[0]
            for i in range(self.otu_gpu.shape[1])
        ]

    # ======================================================
    # Tree utilities
    # ======================================================
    def _check_tree_match(self):
        missing = set(self.otus) - set(self.tip_map.keys())
        if missing:
            raise ValueError(
                f"{len(missing)} OTUs not found in tree tips."
            )

    def _build_cophenetic_matrix(self):
        """
        Fast cophenetic distance matrix
        """
        n = len(self.otus)
        dist = np.zeros((n, n), dtype=np.float32)

        for i in tqdm(range(n), desc="Cophenetic distance"):
            ni = self.tip_map[self.otus[i]]
            for j in range(i + 1, n):
                nj = self.tip_map[self.otus[j]]
                d = ni.distance(nj)
                dist[i, j] = dist[j, i] = d

        return dist

    # ======================================================
    # betaMNTD / betaNTI
    # ======================================================
    @staticmethod
    def _beta_mntd(idx1, idx2, dist):
        sub = dist[idx1][:, idx2]
        return cp.mean(cp.min(sub, axis=1))

    def compute_beta_nti(self):
        n = len(self.samples)

        # ------------------------
        # Observed betaMNTD
        # ------------------------
        obs = np.zeros((n, n), dtype=np.float32)

        for i in tqdm(range(n), desc="Observed betaMNTD"):
            idx_i = self.sample_nonzero[i]
            if idx_i.size == 0:
                continue

            for j in range(i + 1, n):
                idx_j = self.sample_nonzero[j]
                if idx_j.size == 0:
                    continue

                d = self._beta_mntd(
                    idx_i, idx_j, self.phylo_dist_gpu
                )
                val = float(cp.asnumpy(d))
                obs[i, j] = obs[j, i] = val

        # ------------------------
        # Null models (streaming)
        # ------------------------
        null_mean = np.zeros((n, n), dtype=np.float32)
        null_M2 = np.zeros((n, n), dtype=np.float32)

        for k in tqdm(range(self.n_null), desc="Null betaMNTD"):
            perm = np.random.permutation(len(self.otus))
            dist_null = self.phylo_dist_gpu[perm][:, perm]

            cur = np.zeros((n, n), dtype=np.float32)

            for i in range(n):
                idx_i = self.sample_nonzero[i]
                if idx_i.size == 0:
                    continue

                for j in range(i + 1, n):
                    idx_j = self.sample_nonzero[j]
                    if idx_j.size == 0:
                        continue

                    d = self._beta_mntd(
                        idx_i, idx_j, dist_null
                    )
                    val = float(cp.asnumpy(d))
                    cur[i, j] = cur[j, i] = val

            delta = cur - null_mean
            null_mean += delta / (k + 1)
            null_M2 += delta * (cur - null_mean)

        null_sd = np.sqrt(null_M2 / max(self.n_null - 1, 1))
        null_sd[null_sd == 0] = np.nan

        beta_nti = (obs - null_mean) / null_sd
        np.fill_diagonal(beta_nti, np.nan)

        return pd.DataFrame(
            beta_nti,
            index=self.samples,
            columns=self.samples,
        )

    # ======================================================
    # RCbray
    # ======================================================
    def compute_rcbray(self):
        otu = self.otu
        rng = np.random.default_rng(self.random_state)

        n = len(self.samples)
        obs = np.zeros((n, n), dtype=np.float32)

        for i in range(n):
            for j in range(i + 1, n):
                obs[i, j] = obs[j, i] = braycurtis(
                    otu.iloc[:, i], otu.iloc[:, j]
                )

        null = np.zeros((self.rc_null, n, n), dtype=np.float32)

        for k in tqdm(range(self.rc_null), desc="RCbray null"):
            shuffled = otu.apply(
                lambda x: rng.permutation(x.values),
                axis=0,
            )
            for i in range(n):
                for j in range(i + 1, n):
                    null[k, i, j] = null[k, j, i] = braycurtis(
                        shuffled.iloc[:, i],
                        shuffled.iloc[:, j],
                    )

        rc = np.zeros((n, n), dtype=np.float32)
        for i in range(n):
            for j in range(i + 1, n):
                p = (null[:, i, j] < obs[i, j]).mean()
                rc[i, j] = rc[j, i] = 2 * (p - 0.5)

        return pd.DataFrame(
            rc, index=self.samples, columns=self.samples
        )

    # ======================================================
    # Stegen process inference
    # ======================================================
    @staticmethod
    def stegen_process(beta_nti: pd.DataFrame, rcbray: pd.DataFrame):
        """
        Robust Stegen process inference (index-safe)
        """

        bnti = beta_nti.values
        rc = rcbray.values

        n = bnti.shape[0]
        processes = []

        for i in range(n):
            for j in range(i + 1, n):
                b = bnti[i, j]
                r = rc[i, j]

                if np.isnan(b) or np.isnan(r):
                    continue

                if b > 2:
                    processes.append("Variable selection")
                elif b < -2:
                    processes.append("Homogeneous selection")
                else:
                    if r > 0.95:
                        processes.append("Dispersal limitation")
                    elif r < -0.95:
                        processes.append("Homogenizing dispersal")
                    else:
                        processes.append("Drift")

        return pd.Series(processes)


    # ======================================================
    # Group-wise summary (Figure-ready)
    # ======================================================
    def summarize_by_group(self):
        beta_nti = self.compute_beta_nti()
        rcbray = self.compute_rcbray()

        results = {}

        for g in self.meta[self.group_col].unique():
            samples = self.meta.index[
                self.meta[self.group_col] == g
            ]

            b_sub = beta_nti.loc[samples, samples]
            r_sub = rcbray.loc[samples, samples]

            proc = self.stegen_process(b_sub, r_sub)
            results[g] = proc.value_counts(normalize=True)

        return (beta_nti, rcbray, pd.DataFrame(results).fillna(0))




from pathlib import Path 
from loguru import logger
from tqdm import tqdm 

from typing import List, Dict, Tuple, Optional, Union

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns 

import networkx as nx
import random


class NetworkBuilder:
    '''A collection of static methods for building microbial networks from correlation data.
    Containing:
        - Network construction from correlation and p-value matrices
        - Batch processing from directory of files
        - Network saving in edge list format
    Examples:
    >>>G_rhizo = MicrobiomeNetwork.build_network("rhizo_cor.tsv", "rhizo_pvalues.tsv", r_threshold=0.3, p_threshold=0.01)
    >>>G_bulk = MicrobiomeNetwork.build_network("bulk_cor.tsv", "bulk_pvalues.tsv", r_threshold=0.3, p_threshold=0.01)
    '''

    def build_network(self, cor_file: str, p_file: str,  r_threshold: float = 0.3, p_threshold: float = 0.01):
        '''Building network from correlation and p-value matrix caculated by SparCC/FastSpar.
        Args:
            cor_file : str, path to correlation matrix file (tsv)
            p_file : str, path to p-value matrix file (tsv)
            r_threshold : float, correlation coefficient threshold
            p_threshold : float, p-value threshold
        Returns:
            G : networkx.Graph, constructed network

        Example:
        >>>G_rhizo = build_network("rhizo_cor.tsv", "rhizo_pvalues.tsv", "rhizo")
        >>>G_bulk = build_network("bulk_cor.tsv", "bulk_pvalues.tsv", "bulk")
        '''

        cor = pd.read_csv(cor_file, sep='\t', index_col=0)
        pval = pd.read_csv(p_file, sep='\t', index_col=0)

        G = nx.Graph()

        n = len(cor.index)
        assert cor.shape == pval.shape # 确保两个矩阵维度一致

        # 使用位置索引，避免类型比较问题
        for a in range(n):
            for b in range(a + 1, n):  # 避免重复和自我比较
                r = cor.iloc[a, b]
                p = pval.iloc[a, b]

                # if (r >= 0.6 or r <= -0.4) and p < p_threshold:
                if abs(r) >= r_threshold and p < p_threshold:
                    node_i = cor.index[a]
                    node_j = cor.columns[b]
                    G.add_node(node_i, taxa= node_i) # add atribute ['taxa', ...]
                    G.add_node(node_j, taxa= node_j) # add atribute ['taxa', ...]
                    G.add_edge(node_i, node_j, weight=r, sign='positive' if r > 0 else 'negative') # add atribute ['weight', 'sign', ...], sign for gephi.

        return G
    
    def build_network_from_dir(self, dir_path: str, start_with: Tuple[str, ...] = ('Species_',), r_threshold: float = 0.3, p_threshold: float = 0.01):
        '''Batch build networks from a directory containing correlation and p-value files.
        Args:
            dir_path : str, path to the directory containing correlation and p-value files
            start_with : Tuple[str, ...], prefixes to filter subdirectories, e.g. ('Genus_', 'Species_')
            r_threshold : float, correlation coefficient threshold
            p_threshold : float, p-value threshold
        Returns:
            networks : dict, mapping from subdirectory names to constructed networkx.Graph objects      

        Example:
        >>> networks = build_network_from_dir("/path/to/directory", start_with="Genus_", r_threshold=0.3, p_threshold=0.01)
        '''

        dir_path = Path(dir_path)
        networks: Dict[str, nx.Graph] = {}

        for sub in dir_path.iterdir():
            if sub.is_dir() and sub.stem.startswith(start_with):
                cor_file = sub / f"{sub.stem}_cor.tsv"
                pval_file = sub / f"{sub.stem}_pvalues.tsv"

                if cor_file.exists() and pval_file.exists():
                    G = self.build_network(
                        cor_file=cor_file, 
                        p_file=pval_file, 
                        r_threshold=r_threshold,
                        p_threshold=p_threshold
                    )
                    networks[sub.stem] = G
                    print(f"[OK] Finished: {sub.stem} Nodes {G.number_of_nodes()}, Edges {G.number_of_edges()}）")
                else:
                    print(f"[Skipped] {sub.stem} due to missing files")

        return networks
    
    def save_network(self, G: nx.Graph, output_path: str, format: str = 'gexf'):
        '''Save network to a file in edge list format.
        Args:
            G : networkx.Graph, the network to save
            output_path : str, path to save the edge list file
            format : str, file format, 'gexf' or 'edgelist'
        Returns:
            None

        Example:
        >>> save_network(G_rhizo, "rhizo_network.gexf", format='gexf')
        >>> save_network(G_bulk, "bulk_network.edgelist", format='edgelist')
        '''

        if format == 'gexf':
            nx.write_gexf(G, output_path)
        elif format == 'edgelist':
            nx.write_edgelist(G, output_path, data=['weight', 'sign'])
        else:
            raise ValueError(f"Unsupported format: {format}. Supported formats are 'gexf' and 'edgelist'.")
    

class NetworkVisualizer:
    '''A collection of static methods for visualizing microbial networks.
    Containing:
        - Layout generation
        - Node and edge styling
        - Legend creation
    Examples:
    >>>NetworkVisualizer.visualize_network(G)
    '''

    def _auto_spring_layout(self, G: nx.Graph, seed: int = 42):
        """
        适用于微生物组网络：自动调参的 spring layout + 去重叠（push apart）
        返回更均匀、无重叠、稳定的布局，适应不同网络规模与密度。
        """

        def push_apart_positions(pos, min_dist=0.05, iterations=50, force=0.01):
            """减少节点重叠：在原布局基础上迭代推开过近节点。
            pos: dict {node: [x,y]} from spring_layout
            min_dist: 最小允许距离（可根据节点数自动设定）
            iterations: 迭代次数
            force: 推开强度
            """
            nodes = list(pos.keys())
            for _ in range(iterations):
                moved = False
                for i, u in enumerate(nodes):
                    for v in nodes[i+1:]:
                        dx = pos[u][0] - pos[v][0]
                        dy = pos[u][1] - pos[v][1]
                        dist = np.sqrt(dx*dx + dy*dy)

                        if dist < min_dist and dist > 0:
                            # 归一化方向
                            ux, uy = dx/dist, dy/dist
                            # 推开
                            shift = force * (min_dist - dist)
                            pos[u][0] += ux * shift
                            pos[u][1] += uy * shift
                            pos[v][0] -= ux * shift
                            pos[v][1] -= uy * shift
                            moved = True
                if not moved:
                    break
            return pos

        n = G.number_of_nodes()
        m = G.number_of_edges()
        degrees = np.array([d for _, d in G.degree()])
        d_avg = degrees.mean() if n > 0 else 1

        # --------- 自动计算 spring_layout 参数 ------------
        k = 1.2 * (1/np.sqrt(n)) * (1 + d_avg/10)
        iterations = int(200 + n*2 + d_avg*10)

        pos = nx.spring_layout(G, k=k, iterations=iterations, seed=seed)

        # --------- 自动计算去重叠参数 ------------------------
        min_dist = 0.8 / np.sqrt(n)
        pos = push_apart_positions(pos, min_dist=min_dist, iterations=80, force=0.02)

        return pos

    def _get_layout_style(self, layout_type: str, G: nx.Graph = None):
        '''Generate layout for network visualization.
        Args:
            layout_type : str, type of layout ('spring', 'circular', 'shell', 'kamada_kawai', etc.)
            G : networkx.Graph, the network to layout
        '''

        if layout_type == 'spring':
            logger.info("Using spring layout for visualization.")
            # pos = nx.spring_layout(
            #     G, 
            #     k= 3/np.sqrt(G.number_of_nodes()),              # optimal distance between nodes, k_value = 4 / np.sqrt(len(G.nodes))
            #     iterations= 5000,                               # more iterations for better convergence
            #     seed= 42,         # reproducibility
            #     scale= 1,         # scale the layout for better spacing
            # )
            pos = self._auto_spring_layout(G=G)

        elif layout_type == 'force_atlas2':
            logger.info("Using force_atlas2 layout for visualization.")
            try:
                import fa2
            except ImportError:
                raise ImportError("fa2 package not found. Please install it via 'pip install fa2' to use force_atlas2 layout.")

            forceatlas2 = fa2.ForceAtlas2(
                outboundAttractionDistribution=False,
                linLogMode=False,
                adjustSizes=False,
                edgeWeightInfluence=1.0,
                jitterTolerance=1.0,
                barnesHutOptimize=True,
                barnesHutTheta=1.2,
                scalingRatio=10.0,
                strongGravityMode=False,
                gravity=1.0,
                verbose=True
            )
            pos = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=2000)

        elif layout_type == 'circular':
            logger.info("Using circular layout for visualization.")
            pos = nx.circular_layout(G)
        elif layout_type == 'shell':
            logger.info("Using shell layout for visualization.")
            pos = nx.shell_layout(G)
        elif layout_type == 'kamada_kawai':
            logger.info("Using kamada_kawai layout for visualization.")
            pos = nx.kamada_kawai_layout(G)
        else:
            raise ValueError(f"Unsupported layout type: {layout_type}")

        return pos
    
    def _scaled_node_size_by_degree(self, G: nx.Graph, 
                                 standard_type: str = 'sqrt_min_max', min_max: Tuple[int, int] = (1, 200),
        ) -> List[float]:
        '''Calculate node sizes based on their degree.
        Args:
            G : networkx.Graph, the network
        Returns:
            Tuple[List[float], List[float]]: raw node sizes and scaled node sizes
        '''

        if standard_type == 'none':
            degree_dict = dict(G.degree())
            scaled_node_sizes = [degree_dict[n] for n in G.nodes]
        
        elif standard_type == 'min_max':
            from sklearn.preprocessing import MinMaxScaler
            degree_dict = dict(G.degree())
            degrees = np.array(list(degree_dict.values())).reshape(-1, 1)
            scaler = MinMaxScaler(feature_range=min_max)    # scale node sizes between 100 and 1000
            scaled_node_sizes = scaler.fit_transform(degrees).flatten().tolist()

        elif standard_type == 'sqrt':
            degree_dict = dict(G.degree())
            scaled_node_sizes = [np.sqrt(degree_dict[n]) * 100 for n in G.nodes]  # scale by sqrt

        elif standard_type == 'sqrt_min_max':
            from sklearn.preprocessing import MinMaxScaler
            degree_dict = dict(G.degree())
            degrees = np.array([np.sqrt(degree_dict[n]) for n in G.nodes]).reshape(-1, 1)
            scaler = MinMaxScaler(feature_range=min_max)    # scale node sizes between 100 and 1000
            scaled_node_sizes = scaler.fit_transform(degrees).flatten().tolist()
        else:
            raise ValueError(f"Unsupported standard_type: {standard_type}")
        
        degree_dict = dict(G.degree())
        raw_node_sizes: List[float] = [degree_dict[n] for n in G.nodes]

        return (raw_node_sizes, scaled_node_sizes)
    
    def _draw_curved_edges_weight_color(
        self, G, pos, ax=None, curvature: float = 0.2,
        pos_color: str = "red", neg_color: str = "blue",
        width_scale: float = 1.0, alpha: float = 0.7
    ):
        
        from matplotlib.patches import FancyArrowPatch
        
        if ax is None:
            ax = plt.gca()

        # 处理多重边：对同一对节点增加不同曲率，避免重叠
        edge_counter = {}
        for i, (u, v) in enumerate(G.edges()):
            key = tuple(sorted([u, v]))
            edge_counter[key] = edge_counter.get(key, 0) + 1
            rad = curvature * edge_counter[key]

            # 获取边的 weight
            w = G[u][v].get("weight", 0)

            # 画曲线
            arrow = FancyArrowPatch(
                pos[u], pos[v],
                connectionstyle= f"arc3,rad={rad}",
                arrowstyle= '-',                                # 不要箭头
                color= pos_color if w > 0 else neg_color,       # 自动决定颜色（你也可改成 colormap）
                linewidth= abs(w) * width_scale,                # 根据 weight 调整宽度
                alpha= alpha                                    # 调整透明度
            )
            ax.add_patch(arrow)

    def _show_legends(self, ax: plt.Axes, G: nx.Graph, communities_handle: Optional = None):
        '''Generate legend for network visualization based on node and edge attributes.'''

        if ax is None:
            ax = plt.gca()

        # ------------- Node size legend (degree) --------------------
        raw_node_sizes, scaled_node_sizes = self._scaled_node_size_by_degree(G=G)

        raw_sizes = np.array(raw_node_sizes)
        sizes = np.array(scaled_node_sizes)
        if len(sizes) > 0:
            raw_size_vals = np.percentile(raw_sizes, [25, 50, 75])
            size_vals = np.percentile(sizes, [25, 50, 75])
            # vertical line markers for node sizes
            node_size_handles = [
                plt.Line2D(
                    [0], [0],
                    label=f"{int(raw_s)}",
                    marker='o',
                    color='w',
                    markerfacecolor='w',
                    markeredgecolor='black',
                    markersize=np.sqrt(s)       # scale down for legend
                )
                for raw_s, s in zip(raw_size_vals, size_vals)
            ]
        else:
            node_size_handles = []

        # Create legend 1
        # vertival line markers for node sizes
        legend1 = ax.legend(
            handles=node_size_handles,
            title="Node Degree",
            loc="upper left",
            bbox_to_anchor=(1.05, 1),
            frameon=False
        )

        ax.add_artist(legend1)   # <-- THIS IS KEY! Prevent overwrite, pretty important!

        # ------------- Edge sign legend (positive / negative) --------------------
        edge_sign_handles = [
            plt.Line2D([0], [0], color='red', lw=2, label='Positive'),
            plt.Line2D([0], [0], color='blue', lw=2, label='Negative')
        ]

        # Create legend 2
        legend2 = ax.legend(
            handles=edge_sign_handles,
            title="Edge Sign",
            loc="upper left",
            bbox_to_anchor=(1.05, 0.7),  # place below Degree legend
            frameon=False
        )

        if communities_handle is not None:
            ax.add_artist(legend2)
            # ------------- Community legend --------------------
            legend3 = ax.legend(
                handles=communities_handle, 
                title="Communities",
                loc="upper left", 
                bbox_to_anchor=(1.05, 0.4), 
                frameon=False
            )

    def show(
            self, G, 
            pos_style: str = 'spring',
            node_color: Union[str, List[str]] = "skyblue",
            node_label_fontsize: Optional[int] = None,
            edge_width_factor: float = 0.5, edge_line_style: str = 'curved', # straight / curved  
            alpha: float = 0.7,
        ):
        """visualization of microbiome network.    
        Args:
            G : NetworkX Graph
            pos_style : str, type of layout ('spring', 'circular', 'shell', etc.)
            node_color : str or List[str], color for nodes (single color or list of colors per node)
            node_label_fontsize : Optional[int], font size for node labels (None to hide labels)
            edge_width_factor : float, control the scale of edge line width
            alpha : float, transparency level for nodes and edges
        Function:
            - Node size is based on its degree.
            - Edge line width is based on the value of abs(weight).
            - Different color between positive and negative relationship.
        """

        if len(G.nodes) == 0:
            logger.error("[Warning] Empty, skill plot.")
            return

        # ---- 布局： spring等等 ----
        pos = self._get_layout_style(layout_type=pos_style, G=G)

        # ---- 节点： 大小、颜色----
        raw_node_sizes, scaled_node_sizes = self._scaled_node_size_by_degree(G=G)
        
        if isinstance(node_color, str):
            node_colors: List[str] = [node_color for _ in G.nodes]
        elif isinstance(node_color, list) and len(node_color) == len(G.nodes):
            node_colors: List[str] = node_color
        else:
            raise ValueError("ERROR: node_color must be a string or a list with length equal to number of nodes.")

        # ---- 边： 宽度（相关性强度）、颜色----
        if edge_line_style == 'straight':
            logger.info("Using straight edges for visualization.")
            edge_weights: List[float] = [abs(G[u][v]["weight"]) * edge_width_factor for u, v in G.edges]

            edge_colors: List[str] = ["red" if G[u][v]["weight"] > 0 else "blue" for u, v in G.edges]


        # [节点]：大小、颜色
        nx.draw_networkx_nodes(
            G=G, 
            pos= pos,
            node_size= scaled_node_sizes,
            node_color= node_colors,
            alpha= alpha,
            # linewidths=0.5,
            # edgecolors="skyblue" # skyblue, black
        )

        # [节点标签]：名字（微生物名）可能很多，一般不推荐显示
        if node_label_fontsize is not None:
            nx.draw_networkx_labels(
                G=G, pos=pos, 
                font_size= node_label_fontsize
            )

        # [边]：曲线、颜色、宽度
        if edge_line_style == 'curved':
            logger.info("Using curved edges for visualization.")
            self._draw_curved_edges_weight_color(
                G=G, pos=pos, curvature=0.25, 
                pos_color='red', neg_color='blue',
                width_scale= edge_width_factor, alpha= alpha,
            )
        elif edge_line_style == 'straight':
            # [边]：直线、颜色、宽度
            nx.draw_networkx_edges(
                G=G, pos=pos,
                width=edge_weights,
                edge_color=edge_colors,
                alpha= alpha
            )

        # [边标签]：颜色
        # nx.draw_networkx_edge_labels(
        #     G=G, pos=pos,
        #     edge_labels=nx.get_edge_attributes(G, 'weight'),
        #     font_color='green',
        #     font_size=8,
        # )

        self._show_legends(ax=plt.gca(), G=G)

        plt.axis("off")
        plt.title("Microbial Co-occurrence Network", fontsize='large')

    def show_community(self, G: nx.Graph, communities: Optional[List[set]] = None, **kwargs):
        '''Visualize network with nodes colored by their community membership.
        Args:
            G : NetworkX Graph
            communities : list of sets, each set contains nodes in one community
            **kwargs : additional arguments for the show() method
        '''

        # 1. calculate communities if not provided
        if communities is None:
            communities = NetworkMetrics.detect_communities(G=G, algorithm='girvan_newman')

        # 2. Generate a color map for communities
        color_map = plt.cm.get_cmap('tab20', len(communities))
        node_color: List[str] = []
        comm_dict: Dict[str, int] = {}
        for i, comm in enumerate(communities):
            for node in comm:
                comm_dict[node] = i

        for node in G.nodes:
            comm_id = comm_dict.get(node, -1)
            if comm_id >= 0:
                node_color.append(color_map(comm_id))
            else:
                node_color.append('gray')  # 未分配社区的节点为灰色

        # 3. Visualize the network
        self.show(
            G=G, 
            node_color=node_color, 
            node_label_fontsize=kwargs.get('node_label_fontsize'), 
            edge_width_factor=kwargs.get('edge_width_factor'), 
            alpha=kwargs.get('alpha', 0.7),
        )

        # 4. reruning show with legend for communities
        # legend for communities
        from matplotlib.patches import Patch
        handles = [
            Patch(color=color_map(i), label=f'Community {i+1}')
            for i in range(len(communities))
        ]
        self._show_legends(ax=plt.gca(), G=G, communities_handle=handles)


class NetworkMetrics:
    '''A collection of static methods for analyzing microbial networks.
    Containing: 
        - Basic network statistics
        - Node centrality measures
        - Community detection

    Examples:
    >>>stats = NetworkMetrics.compute_basic_metrics(G)
    >>>centrality = NetworkMetrics.compute_centrality(G)
    >>>communities = NetworkMetrics.detect_communities(G)
    '''

    @staticmethod
    def compute_basic_metrics(G: nx.Graph) -> Dict[str, float]:
        """Calculate the basic statistics of a microbial network.
        Args:
            G, NetworkX Graph
        Returns:
            stats: dict, including:
                            - num_nodes,
                            - num_edges, 
                            - avg_degree, 
                            - density, 
                            - avg_clustering, 
                            - avg_path_length

        Examples:
        >>>stats = compute_basic_stats(G)
        >>>print(stats)
        >>>plot_basic_stats(stats)
        """

        stats: Dict[str, float] = {}

        # 节点数
        stats["num_nodes"] = G.number_of_nodes()

        # 边数
        stats["num_edges"] = G.number_of_edges()

        # 平均度
        degrees = dict(G.degree())
        stats["avg_degree"] = np.mean(list(degrees.values()))

        # 网络密度
        stats["density"] = nx.density(G)

        # 平均聚类系数
        stats["avg_clustering"] = nx.average_clustering(G)

        # 连通性处理（必须对最大子图）
        if nx.is_connected(G):
            stats["avg_path_length"] = nx.average_shortest_path_length(G)
        else:
            giant = max(nx.connected_components(G), key=len)
            subgraph = G.subgraph(giant)
            stats["avg_path_length"] = nx.average_shortest_path_length(subgraph)

        return stats
    
    @staticmethod
    def compute_centrality(G: nx.Graph, level: str, metric: str = "degree"):
        '''Calculate centrality measures for nodes in the network and return as a DataFrame.
        Args:
            G : networkx.Graph, the network
            level : str, sample level (e.g., 'Bulk', 'Rhizo') in PatName
            metric : str, centrality metric to compute ('degree', 'betweenness', 'closeness', 'eigenvector')
        Returns:
            centrality_df : pd.DataFrame, with columns ['node', 'centrality', 'Level', 'Metric']
        '''

        if metric == "degree":
            c = nx.degree_centrality(G)
        elif metric == "betweenness":
            c = nx.betweenness_centrality(G)
        elif metric == "closeness":
            c = nx.closeness_centrality(G)
        elif metric == "eigenvector":
            c = nx.eigenvector_centrality(G, max_iter=1000)

        return pd.DataFrame({
            "Node": list(c.keys()),
            "Centrality": list(c.values()),
            "Level": level,
            "Metric": metric
        })
    
    @staticmethod
    def compute_centrality_by_batch(Gs: Dict[str, nx.Graph]):
        '''Calculate centrality with batch
        Args:
            Gs: Dict
            group: group level
        Returns:
            df: pd.DataFrame
        '''

        for name, G in tqdm(Gs.items(), desc="Batch Centrality Calculation", position=0, leave=True):
            taxa = name.split('_')[0]
            group = name.split('_')[1]
            level = name.split('_')[2]
            for metric in tqdm(['degree', 'betweenness', 'closeness', 'eigenvector'], desc=f"Calculating centrality for {name}", position=1, leave=False):
                df_metric = NetworkMetrics.compute_centrality(G, level=level, metric=metric)
                df_metric['Taxa'] = taxa
                df_metric['Group'] = group

                if 'df' not in locals():
                    df = df_metric
                else:
                    df = pd.concat([df, df_metric], axis=0)
        
        return df
    
    @staticmethod
    def detect_communities(G: nx.Graph, algorithm: str = 'girvan_newman') -> List[set]:
        """Detect communities in the network using the Girvan-Newman method.
        Args:
            G, NetworkX Graph
        Returns:
            communities: list of sets, each set contains nodes in one community
        
        Examples:
        >>>communities = detect_communities(G)
        >>>for i, comm in enumerate(communities):
        >>>    print(f"Community {i+1}: {comm}")
        """

        if algorithm == 'girvan_newman':
            from networkx.algorithms.community import girvan_newman

            comp = girvan_newman(G)
            first_level_communities = next(comp)
            communities = [set(c) for c in first_level_communities]

        elif algorithm == 'louvain':
            import community as community_louvain

            partition = community_louvain.best_partition(G)
            comm_dict: Dict[int, set] = {}
            for node, comm_id in partition.items():
                if comm_id not in comm_dict:
                    comm_dict[comm_id] = set()
                comm_dict[comm_id].add(node)
            communities = list(comm_dict.values())

        elif algorithm == 'greedy':
            from networkx.algorithms.community import greedy_modularity_communities
            communities = list(greedy_modularity_communities(G))
            modularity = nx.algorithms.community.modularity(G, communities)
            return (modularity, communities)

        else:
            raise ValueError(f"Unsupported algorithm: {algorithm}. Supported algorithms are 'girvan_newman' and 'louvain'.")
        
        return communities
    
    @staticmethod
    def fit(G: nx.Graph):
        stats: Dict[str, float] = NetworkMetrics.compute_basic_metrics(G)
        centrality: Dict[str, Dict[str, float]] = NetworkMetrics.compute_centrality(G)
        # communities: List[set] = NetworkMetrics.detect_communities(G)

        results = {
            **stats,
            **centrality,
            # **communities,
        }

        return results

    @staticmethod
    def fit_by_batch(Gs: Dict[str, nx.Graph]) -> pd.DataFrame:
        '''Batch compute network metrics for multiple networks.
        Args:
            Gs : dict, mapping from network names to NetworkX Graph objects
        Returns:
            metrics_df : pd.DataFrame, each row corresponds to a network and columns are metrics
        Example:
        >>> metrics_df = fit_by_batch(networks)
        >>> print(metrics_df)
        '''

        records = []
        for name, G in tqdm(Gs.items(), desc="Batch Network Metrics Calculation", position=0, leave=True):
            stats = NetworkMetrics.compute_basic_metrics(G)
            record = {
                'network': name,
                **stats,
            }
            records.append(record)

        metrics_df = pd.DataFrame(records)

        return metrics_df

    @staticmethod
    def show(metrics_df: pd.DataFrame, rank: Optional[str] = None, group: Optional[str] = None, level: Optional[str] = None,
             figsize: Tuple[int, int] = (15, 10)
        ) -> pd.DataFrame:
        '''Visualize network metrics across multiple networks.
        Args:
            metrics_df : pd.DataFrame, output from fit_by_batch()
            rank : Optional[str], filter by taxonomic rank (e.g., 'Genus', 'Species')
            group : Optional[str], filter by sample group (e.g., 'rhizo', 'bulk')
            level : Optional[str], filter by taxonomic level (e.g., 'Phylum', 'Class')
        Return:
            tem_df : pd.DataFrame, melted DataFrame used for plotting
        '''

        tem_df = metrics_df.copy()
        tem_df['rank'] = tem_df['network'].apply(lambda x: x.split('_')[0])
        tem_df['group'] = tem_df['network'].apply(lambda x: x.split('_')[1])
        tem_df['level'] = tem_df['network'].apply(lambda x: x.split('_')[2])

        if rank is not None:
            tem_df = tem_df[tem_df['rank'] == rank]
        if group is not None:
            tem_df = tem_df[tem_df['group'] == group]
        if level is not None:
            tem_df = tem_df[tem_df['level'] == level]

        tem_df = tem_df.melt(
            id_vars=['network', 'rank', 'group', 'level'],
            value_vars=['num_nodes', 'num_edges', 'avg_degree', 'density', 'avg_clustering', 'avg_path_length', ],
            var_name='metric',
            value_name='value'
        )

        # sns.barplot(
        g= sns.catplot(
            data=tem_df,
            x='level',
            y='value',
            hue='level',
            col='metric',
            col_wrap=3,
            kind='bar',
            sharex=False,
            sharey=False,
            palette='husl',
            height=3,
            aspect=1.2,
        )

        return (tem_df, g)

    @staticmethod    
    def robustness_attack(Gs: Dict[str, nx.Graph], n_rep=100, method: str = 'random', noise: float = 1e-6) -> pd.DataFrame:
        '''Calculate AUC for random attack simulations on multiple networks.
        Args:
            Gs : dict, mapping from network names to NetworkX Graph objects
            n_rep : int, number of repetitions for each network
            method : str, 'random' for random attack, 'target' for targeted attack based on degree
            noise : float, small noise to break degree ties in targeted attack
        Returns:
            aucs : List[Dict], each dict contains 'Taxa', 'Group', 'Level', and 'AUC'
        '''

        from sklearn.metrics import auc
        results: Dict[str, List] = {
            'Taxa': [],
            'Group': [],
            'Level': [],
            'Fraction of nodes removed': [],
            'Relative size of largest component': [],
            'repeat': [],
            'AUC': [],
            'Method': [],
        }

        for name, G in tqdm(Gs.items(), desc=f"Random Attack AUC Calculation", position=0, leave=True):
            for repeat_record in tqdm(range(n_rep), desc=f"  Repeat of {name}", position=1, leave=False):
                G_tmp = G.copy()

                if method == 'random':
                    # 随机打乱节点顺序
                    nodes = list(G_tmp.nodes())
                    random.shuffle(nodes) 
                elif method == 'target':
                    # 加微小噪声打破 degree ties，并按照 degree 降序排列节点
                    deg = {
                        n: G_tmp.degree(n) + random.uniform(0, noise)
                        for n in G_tmp.nodes()
                    }
                    nodes = sorted(deg, key=deg.get, reverse=True)

                x, y = [], []
                # 逐步移除节点，计算最大连通分量比例；
                # 都是同一个 G_tmp, 有累计的路径效果（即被移除节点顺序会影响最终AUC计算的结果）
                for i, node in enumerate(nodes):
                    if i > 0:
                        G_tmp.remove_node(node)

                    # 计算最大连通分量比例作为网络连通性指标
                    largest_cc = max(
                        (len(c) for c in nx.connected_components(G_tmp)),
                        default=0
                    )
                    x.append(i / len(nodes))
                    y.append(largest_cc / len(nodes))

                length = len(x)
                results['Taxa'].extend([name.split('_')[0]] * length)
                results['Group'].extend([name.split('_')[1]] * length)
                results['Level'].extend([name.split('_')[2]] * length)
                results['Fraction of nodes removed'].extend(x)
                results['Relative size of largest component'].extend(y)
                results['repeat'].extend([repeat_record + 1] * length)
                results['AUC'].extend([auc(x, y)] * length)
                results['Method'].extend([method] * length)

        results_df: pd.DataFrame = pd.DataFrame(results)

        return results_df
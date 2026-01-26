import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D


class OrdinationPlotter:
    """
    一个适用于所有 skbio 排序方法（PCA, PCoA, CCA, RDA 等）的可视化类。
    支持：
        - 样本点绘制（离散/连续变量着色）
        - 环境变量箭头
        - 物种箭头
        - Group 置信椭圆
        - 自动图例
        - 自动解释度显示
    """

    def __init__(self, ordination, metadata=None):
        """
        参数
        -------
        ordination : skbio OrdinationResults
            包括 PCA, PCoA, CCA, RDA 等的结果对象

        metadata : pandas DataFrame
            索引与 ordination.samples 一致，用于映射样本信息
        """

        self.ord = ordination
        self.metadata = metadata

        # 样本 (scores)
        self.samples = ordination.samples

        # 物种 (feature scores)，只有 CCA / RDA / CA 才有
        self.features = getattr(ordination, "features", None)

        # 环境变量箭头（biplot scores）
        self.env = getattr(ordination, "biplot_scores", None)

        # 解释度（如果存在）
        if hasattr(ordination, "proportion_explained"):
            self.explained = ordination.proportion_explained * 100
        else:
            self.explained = np.array([np.nan] * self.samples.shape[1])

    def plot(self,
             color_by=None,
             continuous_color=False,
             add_ellipse=False,
             feature_arrows=False,
             sample_labels=False,
             title=None):
        """
        参数
        -------
        color_by : str
            metadata 中用于样本着色的列名（分类或连续）

        continuous_color : bool
            若 True，则颜色以渐变方式映射连续变量

        add_ellipse : bool
            是否在分类变量组间绘制置信椭圆

        feature_arrows : bool
            是否绘制物种箭头（CCA/RDA）

        sample_labels : bool
            是否显示样本标签

        title : str
            图标题
        """

        fig, ax = plt.subplots(figsize=(5, 5))

        # ------------------------------------------------------------
        # 1. 样本点
        # ------------------------------------------------------------
        if self.metadata is not None and color_by is not None:
            values = self.metadata.loc[self.samples.index][color_by]

            if continuous_color:  # 连续变量渐变颜色
                sc = ax.scatter(
                    self.samples.iloc[:, 0],
                    self.samples.iloc[:, 1],
                    c=values,
                    s=70,
                    cmap="viridis",
                    edgecolor="black"
                )
                plt.colorbar(sc, label=color_by)

                color_map = None  # 连续变量不使用分类图例

            else:  # 分类变量
                categories = pd.Categorical(values)
                sc = ax.scatter(
                    self.samples.iloc[:, 0],
                    self.samples.iloc[:, 1],
                    c=categories.codes,
                    s=70,
                    cmap="tab20",
                    edgecolor="black"
                )

                # 分类图例
                legend1 = ax.legend(
                    handles=[
                        Line2D([0], [0],
                               marker='o', color='w',
                               markerfacecolor=plt.cm.tab20(i / len(categories.categories)),
                               label=str(cat), markersize=10)
                        for i, cat in enumerate(categories.categories)
                    ],
                    title=color_by,
                    loc="upper right",
                    frameon=False,
                    bbox_to_anchor=(1.05, 1),
                )
                ax.add_artist(legend1)

        else:
            ax.scatter(self.samples.iloc[:, 0], self.samples.iloc[:, 1],
                       color="steelblue", s=70, edgecolor="black")

        # ------------------------------------------------------------
        # 2. 置信椭圆（仅分类变量）
        # ------------------------------------------------------------
        if add_ellipse and (not continuous_color) and (color_by is not None):
            categories = self.metadata.loc[self.samples.index][color_by]
            for cat in categories.unique():
                pts = self.samples[categories == cat].iloc[:, :2].values
                if len(pts) < 3:
                    continue
                self._draw_ellipse(ax, pts)

        # ------------------------------------------------------------
        # 3. 样本标签
        # ------------------------------------------------------------
        if sample_labels:
            for name, (x, y) in self.samples.iloc[:, :2].iterrows():
                ax.text(x, y, name, fontsize=9)

        # ------------------------------------------------------------
        # 4. 环境变量箭头（CCA、RDA）
        # ------------------------------------------------------------
        if self.env is not None:
            for var, (x, y) in self.env.iloc[:, :2].iterrows():
                ax.arrow(0, 0, x, y,
                         color="red", width=0.002,
                         head_width=0.04,
                         length_includes_head=True)
                ax.text(x * 1.12, y * 1.12, var,
                        fontsize=12, color="red", fontweight="bold")

        # ------------------------------------------------------------
        # 5. 物种箭头（可选）
        # ------------------------------------------------------------
        if feature_arrows and self.features is not None:
            for var, (x, y) in self.features.iloc[:, :2].iterrows():
                ax.arrow(0, 0, x, y,
                         color="gray", alpha=0.5, width=0.001)
                ax.text(x, y, var, fontsize=6, alpha=0.6)

        # ------------------------------------------------------------
        # 6. 坐标轴 + 解释度
        # ------------------------------------------------------------
        pc1, pc2 = self.explained[:2]

        ax.axhline(0, color="gray", linewidth=0.5)
        ax.axvline(0, color="gray", linewidth=0.5)

        xlabel = f"Axis 1 ({pc1:.2f}%)" if not np.isnan(pc1) else "Axis 1"
        ylabel = f"Axis 2 ({pc2:.2f}%)" if not np.isnan(pc2) else "Axis 2"

        ax.set_xlabel(xlabel, fontsize='medium')
        ax.set_ylabel(ylabel, fontsize='medium')

        # ------------------------------------------------------------
        # 7. 标题
        # ------------------------------------------------------------
        if title is None:
            title = getattr(self.ord, "short_method_name", "Ordination")
        ax.set_title(title, fontsize='large', fontweight='bold')

        ax.set_aspect("equal", adjustable="box")

        plt.tight_layout()
        plt.show()

    # ================================================================
    # ⭕ 工具函数：绘制置信椭圆
    # ================================================================
    @staticmethod
    def _draw_ellipse(ax, points, alpha=0.3, edgecolor="black"):
        """
        使用协方差矩阵绘制 95% 置信椭圆
        """
        covariance = np.cov(points, rowvar=False)
        eigenvals, eigenvecs = np.linalg.eigh(covariance)

        order = eigenvals.argsort()[::-1]
        eigenvals, eigenvecs = eigenvals[order], eigenvecs[:, order]

        angle = np.degrees(np.arctan2(*eigenvecs[:, 0][::-1]))

        # chi2 = 5.991 对应 95% CI
        width, height = 2 * np.sqrt(eigenvals * 5.991)

        mean = points.mean(axis=0)

        ell = Ellipse(xy=mean, width=width, height=height,
                      angle=angle, edgecolor=edgecolor,
                      facecolor=edgecolor, alpha=alpha)

        ax.add_patch(ell)

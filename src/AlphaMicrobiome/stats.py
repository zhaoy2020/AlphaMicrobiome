from loguru import logger
from typing import List, Tuple, Dict, Optional, Union, Any

import pandas as pd
import numpy as np
import pingouin as pg


class AutoDiffTest:
    """
    自动差异检验类
    支持：
    - 自动正态性检验和方差齐性检验
    - 两组比较（t-test / Welch / Mann–Whitney 自动选择）
    - 多组比较（ANOVA / Welch ANOVA / Kruskal 自动选择）
    - 自动事后多重比较（Tukey / Dunn）
    """

    def __init__(self, data: pd.DataFrame, value_col: str, group_col: str):
        """
        Parameters
        ----------
        data : pd.DataFrame
            包含数值列和分组列的数据框
        value_col : str
            数值列名称（如 Alpha diversity）
        group_col : str
            分组列名称（如 Treatment）
        """
        self.data = data
        self.value_col = value_col
        self.group_col = group_col
        self.groups = data[group_col].unique()

        # 存储结果
        self.normality_results = None                   # 正态性检验结果
        self.homo_results = None                        # 方差齐性检验结果
        self.test_selected: Optional[str] = None        # 选择的统计检验方法名称
        self.main_result = None                         # 主要检验结果
        self.posthoc_result = None                      # 事后多重比较结果

    # ---------------------------------------------------------------
    # 内部函数：检查分布
    # ---------------------------------------------------------------

    def check_normality(self):
        """对每组进行正态性检验（Shapiro），返回结果表"""
        logger.info(f'Start normality check.')
        res = pg.normality(self.data, dv=self.value_col, group=self.group_col)
        self.normality_results = res
        return res

    def check_homoscedasticity(self):
        """检验各组方差齐性（Levene）"""
        logger.info(f'Start homoscedasticity check.')
        res = pg.homoscedasticity(self.data, dv=self.value_col, group=self.group_col)
        self.homo_results = res
        return res

    # ---------------------------------------------------------------
    # 核心：自动选择统计方法
    # ---------------------------------------------------------------

    def run(self):
        """
        自动执行：
        1. 正态性检验
        2. 方差齐性检验
        3. 选择正确统计方法
        4. 执行事后检验（如果多于 2 组）
        """
        self.check_normality()
        self.check_homoscedasticity()

        n_groups = len(self.groups)
        is_normal = (self.normality_results["normal"].all())
        is_homo = bool(self.homo_results["equal_var"].iloc[0])

        # ---------------------------------------------------
        # 情况 1：两组比较
        # ---------------------------------------------------
        if n_groups == 2:
            g1, g2 = self.groups

            if is_normal:
                if is_homo:
                    test = "Student t-test"
                    res = pg.ttest(
                        self.data[self.data[self.group_col] == g1][self.value_col],
                        self.data[self.data[self.group_col] == g2][self.value_col],
                        correction=False  # 不做 Welch
                    )
                else:
                    test = "Welch t-test"
                    res = pg.ttest(
                        self.data[self.data[self.group_col] == g1][self.value_col],
                        self.data[self.data[self.group_col] == g2][self.value_col],
                        correction=True  # Welch
                    )
            else:
                test = "Mann-Whitney U"
                res = pg.mwu(
                    self.data[self.data[self.group_col] == g1][self.value_col],
                    self.data[self.data[self.group_col] == g2][self.value_col]
                )

            self.test_selected = test
            self.main_result = res
            return res

        # ---------------------------------------------------
        # 情况 2：多组比较
        # ---------------------------------------------------
        else:
            if is_normal:
                if is_homo:
                    test = "ANOVA"
                    res = pg.anova(data=self.data, dv=self.value_col, between=self.group_col)
                else:
                    test = "Welch ANOVA"
                    res = pg.welch_anova(dv=self.value_col, between=self.group_col, data=self.data)
            else:
                test = "Kruskal-Wallis"
                res = pg.kruskal(data=self.data, dv=self.value_col, between=self.group_col)

            self.test_selected = test
            self.main_result = res

            # ---------------------------------------------------
            # 事后多重比较
            # ---------------------------------------------------
            if test in ["ANOVA"]:
                self.posthoc_result = pg.pairwise_tukey(
                    data=self.data, dv=self.value_col, between=self.group_col
                )
            else:
                # 参数：p-adjust = 多重比较校正方法
                self.posthoc_result = pg.pairwise_tests(
                    data=self.data,
                    dv=self.value_col,
                    between=self.group_col,
                    parametric=(test != "Kruskal-Wallis"),
                    padjust="fdr_bh"
                )

            return res

    # ---------------------------------------------------
    # 导出汇总（方便报告）
    # ---------------------------------------------------
    def summary(self):
        return {
            "normality": self.normality_results,
            "homoscedasticity": self.homo_results,
            "test_used": self.test_selected,
            "main_result": self.main_result,
            "posthoc": self.posthoc_result
        }


class DBRDA:

    def __init__(self, otu_table: pd.DataFrame, metadata_table: pd.DataFrame):
        pass


class PERMANOVA:

    def __init__(self, otu_table: pd.DataFrame, metadata_table: pd.DataFrame):
        pass


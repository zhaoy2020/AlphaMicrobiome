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


class AddLetter:
    def get_p_value(self, group1: str, group2: str, diff_matrix: pd.DataFrame):
        return diff_matrix.loc[group1, group2]

    def show_letter(self, df: pd.DataFrame, group: str, metric: str = 'shannon'):

        # =========================
        # Step 1: Kruskal-Wallis
        # =========================
        dif_df = pg.kruskal(data=df, dv=metric, between=group)

        if dif_df['p-unc'].values[0] < 0.05:

            # =========================
            # Step 2: pairwise comparisons
            # =========================
            pairwise_df = pg.pairwise_tests(
                data=df,
                dv=metric,
                between=group,
                parametric=False,
                padjust='fdr_bh'
            )

            # =========================
            # Step 3: 排序（从大到小）
            # =========================
            mean_df = df.groupby(by=group)[metric].mean().reset_index()
            mean_df = mean_df.sort_values(by=metric, ascending=False).reset_index(drop=True)
            cats = mean_df[group].tolist()

            # =========================
            # Step 4: 构建p值矩阵
            # =========================
            diff_matrix = pd.DataFrame(np.nan, index=cats, columns=cats)

            for _, row in pairwise_df.iterrows():
                g1, g2 = row['A'], row['B']
                pval = row['p-corr']  # ✅ 用校正后的p值
                diff_matrix.loc[g1, g2] = pval
                diff_matrix.loc[g2, g1] = pval

            # =========================
            # Step 5: 初始化
            # =========================
            letters = [chr(ord('a') + i) for i in range(26)]
            cat_letter_dict = {cat: "" for cat in cats}

            letter_idx = 0

            # =========================
            # Step 6: 主循环（逐组处理）
            # =========================
            for i, base_cat in enumerate(cats):

                # 如果还没标记 → 给当前字母
                if cat_letter_dict[base_cat] == "":
                    cat_letter_dict[base_cat] += letters[letter_idx]

                current_letter = letters[letter_idx]

                # =========================
                # 向下比较
                # =========================
                for j in range(i + 1, len(cats)):
                    down_cat = cats[j]

                    pval = self.get_p_value(base_cat, down_cat, diff_matrix)

                    if pd.isna(pval) or pval >= 0.05:
                        # 不显著 → 共享字母
                        if current_letter not in cat_letter_dict[down_cat]:
                            cat_letter_dict[down_cat] += current_letter
                    else:
                        # 显著 → 新字母
                        letter_idx += 1
                        new_letter = letters[letter_idx]

                        if new_letter not in cat_letter_dict[down_cat]:
                            cat_letter_dict[down_cat] += new_letter

                        # =========================
                        # 向上回溯（关键修复）
                        # =========================
                        for k in range(0, j):
                            up_cat = cats[k]
                            pval_up = self.get_p_value(down_cat, up_cat, diff_matrix)

                            # 不显著 → 应该共享新字母
                            if pd.isna(pval_up) or pval_up >= 0.05:
                                if new_letter not in cat_letter_dict[up_cat]:
                                    cat_letter_dict[up_cat] += new_letter

                        break  # 必须跳出当前向下循环

            # =========================
            # Step 7: 清理字母（排序去重）
            # =========================
            for cat in cat_letter_dict:
                cat_letter_dict[cat] = ''.join(sorted(set(cat_letter_dict[cat])))

            return diff_matrix, cat_letter_dict

        else:
            print(f"No significant differences among groups for {metric} "
                f"(p={dif_df['p-unc'].values[0]:.2e})")
            
            return None, None

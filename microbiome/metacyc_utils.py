from bs4 import BeautifulSoup
import csv


def parse_html(file_path):
    """
    解析 HTML 文件，递归提取所有层级的分类信息
    :param file_path: HTML 文件路径
    :return: 分类信息列表，每个元素包含层级的路径
    """
    with open(file_path, "r", encoding="utf-8") as file:
        soup = BeautifulSoup(file, "html.parser")
    
    # 获取根节点列表（例如第一个 <ul> 标签）
    ul_tags = soup.find_all("ul", class_="jstree-children")
    
    categories = []
    
    def extract_categories(li_tag, current_path):
        """
        递归提取每个层级的分类信息
        :param li_tag: 当前 <li> 标签
        :param current_path: 当前分类的路径
        """
        # 获取当前分类名称
        a_tag = li_tag.find("a")
        if a_tag:
            category_name = a_tag.get_text(strip=True)
            # 拼接当前路径
            new_path = current_path + [category_name]
            
            # 如果有子节点（<ul>），递归调用
            ul_tag = li_tag.find("ul")
            if ul_tag:
                for child_li in ul_tag.find_all("li", recursive=False):
                    extract_categories(child_li, new_path)
            else:
                # 没有子节点则记录当前路径
                categories.append(new_path)
    
    # 遍历所有根级别的 <ul> 标签
    for ul in ul_tags:
        for li in ul.find_all("li", recursive=False):
            extract_categories(li, [])
    
    return categories

def save_to_csv(categories, output_file):
    """
    将分类信息保存为 CSV 文件
    :param categories: 分类信息列表
    :param output_file: 输出 CSV 文件路径
    """
    # 确定 CSV 列的最大长度（即最大层级数）
    max_length = max(len(category) for category in categories)
    
    # 填充列，使所有行列数一致
    for category in categories:
        category.extend([None] * (max_length - len(category)))
    
    # 写入 CSV 文件
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        # 写入表头（列名：Level1, Level2, ...）
        writer.writerow([f"Level{i+1}" for i in range(max_length)])
        # 写入每一行分类
        writer.writerows(categories)

# 主程序
def get_metacyc_pathway_maps(html_file_path: str, output_csv_path: str):
    # html_file_path = "results/functional_analysis/deseq2_diff/annotation_files/MetaCyc Pathways.html"  # 替换为你的 HTML 文件路径
    # output_csv_path = "output_categories.csv"  # 输出的 CSV 文件路径

    # 解析 HTML 文件，获取分类层级
    categories = parse_html(html_file_path)

    # 保存为 CSV 文件
    save_to_csv(categories, output_csv_path)

    print(f"分类层级已经成功保存到 {output_csv_path}")

if __name__ == "__main__":
    get_metacyc_pathway_maps(
        html_file_path="results/functional_analysis/deseq2_diff/annotation_files/MetaCyc Pathways.html", 
        output_csv_path="output_categories.csv"
    )
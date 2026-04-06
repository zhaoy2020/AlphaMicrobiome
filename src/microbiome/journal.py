def mm_to_inch(mm):
    '''将毫米转换为英寸'''

    # 1 inch = 25.4 mm
    # width / height = 4:3
    
    return mm / 25.4


config = {
        'svg.fonttype': 'none',     # svg 矢量图字体设置,方便在AI中编辑
        'pdf.fonttype': 42,         # pdf 矢量图字体设置,方便在AI中编辑
        'ps.fonttype': 42,          # ps 矢量图字体设置,方便在AI中编辑
    }


nature_config = {
    # === Figure size ===
    'figure.figsize': (mm_to_inch(170), mm_to_inch(85 * 0.75)),  # 85mm x 63.75mm

    # === Font ===
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],  # 无衬线体
    'font.size': 7,

    # === Axes ===
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'axes.linewidth': 0.8,

    # === Ticks ===
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,

    # === Lines ===
    "lines.linewidth": 1.0,
    "lines.markersize": 3,

    # ===== Legend =====
    "legend.fontsize": 6,
    "legend.frameon": False,
    "legend.borderaxespad": 0.2,

    # ===== PDF/PS font embedding =====
    "pdf.fonttype": 42,   # TrueType
    "ps.fonttype": 42,
    'svg.fonttype': 'none',  # Keep text as text in SVG
}
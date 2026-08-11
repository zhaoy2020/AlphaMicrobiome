
from pathlib import Path
from loguru import logger

import tkinter as tk 
from tkinter import filedialog, messagebox, scrolledtext

import usearch
from functools import partial 


class UsearchPipeline:
    def __init__(self, data_dir: Path, results_dir: Path):
        self.data_dir = data_dir
        self.results_dir = results_dir 
        
        self.command_dict: dict = {
            'quality_control': partial(usearch.quality_control, data_dir, results_dir),
            'merge_paired_reads': partial(usearch.merge_paired_end_reads, data_dir, results_dir),
            'merge_all_samples': partial(usearch.merge_all_samples, results_dir),
            'filter': partial(usearch.filter_low_quality_reads, input_dir=results_dir.joinpath('003_all_samples_merged')),
            'dereplication': partial(usearch.dereplication, input_dir=results_dir.joinpath('003_all_samples_merged')),
            'remove_singletons': partial(usearch.remove_singletons, input_dir=results_dir.joinpath('003_all_samples_merged')),
            'preorder': partial(usearch.preorder, input_dir=results_dir.joinpath('003_all_samples_merged'), minisize=8),
            'cluster_to_otus': partial(usearch.cluster_to_otus, input_dir=results_dir.joinpath('003_all_samples_merged')),
            'denoising': partial(usearch.denoising, input_dir=results_dir.joinpath('003_all_samples_merged')),
            'otu_table': partial(usearch.build_otu_feature_table, input_dir=results_dir.joinpath('003_all_samples_merged')),
            'zotu_table': partial(usearch.build_zotu_feature_table, input_dir=results_dir.joinpath('003_all_samples_merged')),
        }

        return None 
    
    def run(self, command_name: str):
        self.command_dict[command_name]() # partial_function(), directly call.

        return None 
    
    def one_step_run(self):
        for command_name, command_function in self.command_dict.items():
            command_function() # partial_function(), directly call.

        return None 


class UsearchPage:
    def __init__(self, frame):
        self.frame = frame
        print('UsearchPage initialized')
        # self.pipeline = UsearchPipeline()

    def create(self):
        usearch_frame = tk.Frame(self.frame)
        usearch_frame.grid(row=0, column=0, sticky='nsew')
        # Grid configures: 6 x 3
        rows: int = 12
        columns: int = 3
        for r in range(rows):
            usearch_frame.rowconfigure(r, weight=0)
        for c in range(columns):
            usearch_frame.columnconfigure(c, weight= 3 if c == 2 else 0)

        # Widgets
        datas_button = tk.Button(usearch_frame, text='Datas dir.', command= lambda: self.load_directory(datas_label, 'data_dir')) 
        results_button = tk.Button(usearch_frame, text='Results dir.', command=lambda: self.load_directory(results_label, 'results_dir'))
        datas_label = tk.Label(usearch_frame, text='./datas')
        results_label = tk.Label(usearch_frame, text='./results')
        quality_control_button = tk.Button(usearch_frame, text='Quality control', command=lambda: self.run_command('quality_control'))
        merge_paired_reads_button = tk.Button(usearch_frame, text='Merge paired-end reads', command=lambda: self.run_command('merge_paired_end_reads'))
        merge_all_samples_button = tk.Button(usearch_frame, text='Merge all samples', command=lambda: self.run_command('merge_all_samples'))
        filter_button = tk.Button(usearch_frame, text='Filter', command=lambda: self.run_command('filter_low_quality_reads'))
        dereplication_button = tk.Button(usearch_frame, text='Dereplication', command=lambda: self.run_command('dereplication'))
        remove_singletons_button = tk.Button(usearch_frame, text='Remove singletons', command=lambda: self.run_command('remove_singletons'))
        preorder_button = tk.Button(usearch_frame, text='Preorder', command=lambda: self.run_command('preorder'))
        cluster_to_otus_button = tk.Button(usearch_frame, text='Cluster to OTUs', command=lambda: self.run_command('cluster_to_otus'))
        denoising_button = tk.Button(usearch_frame, text='Denoising', command=lambda: self.run_command('denoising'))
        build_otu_feature_table_button = tk.Button(usearch_frame, text='Build OTU feature table', command=lambda: self.run_command('build_otu_feature_table'))
        build_zotu_feature_table_button = tk.Button(usearch_frame, text='Build ZOTU feature table', command=lambda: self.run_command('build_zotu_feature_table'))
        one_step_button = tk.Button(usearch_frame, text='One step', command=lambda: self.run_command('one_step'))

        command_text = scrolledtext.ScrolledText(usearch_frame)

        # Layout with grid
        datas_button.grid(row=0, column=0, sticky='ew')
        results_button.grid(row=1, column=0, sticky='ew')
        datas_label.grid(row=0, column=1, sticky='ew')
        results_label.grid(row=1, column=1, sticky='ew')
        quality_control_button.grid(row=2, column=0, columnspan=2, sticky='ew')
        merge_paired_reads_button.grid(row=3, column=0, columnspan=2, sticky='ew')
        merge_all_samples_button.grid(row=4, column=0, columnspan=2, sticky='ew')
        filter_button.grid(row=5, column=0, columnspan=2, sticky='ew')
        dereplication_button.grid(row=6, column=0, columnspan=2, sticky='ew')
        remove_singletons_button.grid(row=7, column=0, columnspan=2, sticky='ew')
        preorder_button.grid(row=8, column=0, columnspan=2, sticky='ew')
        cluster_to_otus_button.grid(row=9, column=0, columnspan=1, sticky='ew')
        denoising_button.grid(row=9, column=1, columnspan=1, sticky='ew')
        build_otu_feature_table_button.grid(row=10, column=0, columnspan=1, sticky='ew')
        build_zotu_feature_table_button.grid(row=10, column=1, columnspan=1, sticky='ew')
        one_step_button.grid(row=11, column=0, columnspan=2, sticky='ew')

        command_text.grid(row=0, column=2, rowspan=rows, sticky='nsew')

        return usearch_frame
    
    def load_directory(self, label: tk.Label, type: str):
        if type == 'data_dir':
           data_dir = filedialog.askdirectory(title='Select data directory')
           label.config(text=data_dir)
           self.data_dir = Path(data_dir)
           logger.info(f'Selected data directory: {self.data_dir}')
        elif type == 'results_dir':
            results_dir = filedialog.askdirectory(title='Select results directory')
            label.config(text=results_dir)
            self.results_dir = Path(results_dir)
            logger.info(f'Selected results directory: {self.results_dir}')
        else:
            logger.error(f'Unknown directory type: {type}')

        return None

    def run_command(self, command_name: str):
        pipeline = UsearchPipeline(self.data_dir, self.results_dir)
        if command_name == 'one_step':
            pipeline.one_step_run()
        else:
            messagebox.showinfo(command_name, f'Running {command_name}... This feature is not implemented yet.')

        return None


class GUI:
    def __init__(self, root):
        self.root = root
        self.root.title('Microbiome Analysis Tool')
        self.root.geometry('800x600')
        self.root.minsize(800, 600)

        # Grid configuration, 2 x 2
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_rowconfigure(1, weight=0)
        self.root.grid_columnconfigure(0, weight=0)
        self.root.grid_columnconfigure(1, weight=1)

        self.create_menu()
        self.create_status_bar()
        # Create main content area
        self.create_main_content()

        return None
    
    def create_menu(self):
        # Main menu bar
        menu_bar = tk.Menu(self.root)
        # Submenu
        file_menu = tk.Menu(menu_bar, tearoff=False)
        file_menu.add_command(label='Open')
        file_menu.add_command(label='Save')
        file_menu.add_separator()
        file_menu.add_command(label='Exit', command=self.root.quit)
        menu_bar.add_cascade(label='File', menu=file_menu)
        # Submenu
        microbiome_menu = tk.Menu(menu_bar, tearoff=False)
        microbiome_menu.add_command(label='Usearch', command=lambda: self.show_page(self.usearch))
        microbiome_menu.add_command(label='Home', command=lambda: self.show_page(self.home))
        menu_bar.add_cascade(label='Amplicon', menu=microbiome_menu)
        # Submenu
        help_menu = tk.Menu(menu_bar, tearoff=False)
        help_menu.add_command(label='About', command=lambda: messagebox.showinfo('About', 'Microbiome Analysis Tool v0.1.0'))
        help_menu.add_command(label='Check for updates', command=lambda: messagebox.askyesno('Check for updates', 'This feature is not implemented yet. Do you want to check the GitHub repository?') and self.open_github())
        menu_bar.add_cascade(label='Help', menu=help_menu)

        # Bind the menu bar to the root window
        self.root.config(menu=menu_bar)

        return None
    
    def create_main_content(self):
        main_frame = tk.Frame(self.root)
        main_frame.grid(row=0, column=0, columnspan=2, sticky='nsew')

        self.usearch = self.usearch_page(main_frame)
        self.home = self.home_page(main_frame)
        for page in (self.usearch, self.home):
            page.place(relx=0, rely=0, relwidth=1, relheight=1)

        return None
    
    def home_page(self, main_frame):
        home = tk.Frame(main_frame)
        tk.Label(home, text='Welcome to the Microbiome Analysis Tool!').pack()

        return home
    
    def usearch_page(self, main_frame):
        usearch = UsearchPage(main_frame).create()

        return usearch
    
    def show_page(self, page):
        print(f'Showing page: {page}')
        page.tkraise()

        return None

    def create_status_bar(self):
        status_bar = tk.Frame(self.root, height=30, bg='red')
        status_bar.grid(row=1, column=0, columnspan=2, sticky='ew')
        tk.Label(status_bar, text='Ready').pack(side='left')

        return None

if __name__ == '__main__':
    root = tk.Tk()
    app = GUI(root)
    root.mainloop()
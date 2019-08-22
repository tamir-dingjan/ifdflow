#!/usr/bin/env python
# coding: utf-8


class Generation:
    """
    The Generation object indexes a collection of screens and provides a wrapper for functions
    that are run on all the screens.
    """
    
    def __init__(self, n_processes, members, library, n_actives, n_decoys):
        """
        Initialise a Generation object with the provided models. These can be docking output files
        (.maegz) or pregenerated grid (.zip) files.
        
        Parameters
        ----------
        n_processes : int
            Number of simultaneous processes through which to run.
        members : list
            List of screen filepaths.
        library : str
            Filepath for the screening library.
        n_actives : int
            Number of actives.
        n_decoys : int
            Number of decoys.
        """
        self.n_processes = n_processes
        
        self.members = [os.path.realpath(x) for x in members]
        self.workdir = os.path.dirname(self.members[0]) + "/"
        
        self.library = os.path.realpath(library)
        self.n_actives = n_actives
        self.n_decoys = n_decoys
        
        self.pool = pp.ProcessPool(self.n_processes)

    
    def pool_min(self, min_input):
        """
        Parallel utility to write the minimization job input.
        """
        (x, ligand) = min_input
        m = Model(x, ligand)
        # If the model can't load, there is no pose for the specified ligand.
        # If the ligand isn't found in the file, then no pose for the ligand.
        # In these cases, don't want to run the minimization job, so don't 
        # return any filename.
        m.load()
        if (m.load_status & m.ligandmatch):
            min_name = m.minimize()
            return min_name
        
    
    def write_minimize(self, ligand):
        """
        Write the minimization run files for the generation members.
        
        Parameters
        ----------
        ligand : str
            Title for the ligand to pass to Model. 
            
        """
        
        min_input = [(x, ligand) for x in self.members]
        
        if not (os.path.splitext(self.members[0])[-1] == '.zip'):
            self.pool.clear()
            self.pool.restart()
            output = self.pool.map(self.pool_min, min_input)
            self.pool.close()

            self.mins = output
            
        else:
            print("The loaded models are pregenerated grids (.zips) - No minimization files written.")
        
    
    def prepare_next_generation(self, pattern, name):
        """
        Rename the minimization run files to specify the names for the next generation of model files
        that minimization produces. If the minimization output file already exists, or if
        the required input file is not present, the minimization run file is not renamed but rather
        stored in self.min_exists or self.min_unprepared, respectively.

        Parameters
        ----------
        name : str
            Label to insert into the filename

        """
        
        self.min_exists = []
        self.min_unprepared = []
        self.next_gen = []
        
        for x in self.mins:
            if (x != None):
                ext = os.path.splitext(x)[-1]

                rename = os.path.basename(x).split(pattern)[0]+name+"-min"+ext
                output = os.path.basename(x).split(pattern)[0]+name+"-min-out.maegz"
                
                # Check if the minimization parameter file already exists
                if (len(glob.glob(self.workdir+output)) == 1):
                    self.min_exists.append(self.workdir+rename)
                # Check if the required input file is present
                elif (len(glob.glob(self.workdir+os.path.basename(x).split(pattern)[0]+"-*complex.mae")) == 0):
                    self.min_unprepared.append(x)
                else:
                    shutil.copy2(x, self.workdir+rename)
                    self.next_gen.append(self.workdir+rename)

        self.mins = self.next_gen
    
    
    def run_minimize(self):
        """
        Perform the minimization runs.
        """
        
        jobcodes = []
        
        if self.mins:
            for f in self.mins:
                submitted=False
                while not submitted:
                    if (self.get_used_licenses('Users of PSP_PLOP') < 210):
                        jobid = subprocess.run(["/home/data/Schrodinger2019_1/prime",
                                                "-HOST",
                                                "dragon_serial_mashaq",
                                                f], cwd = self.workdir, text=True, capture_output=True)
                        jobcodes.append(jobid.stdout.split()[-1])
                        submitted = True
                        time.sleep(5)
                    else:
                        time.sleep(10)
        
        self.mins_jobcodes = jobcodes
        

    def pool_split_by_mol(self, x):
        """
        Parallel worker function to split the provided minimization output file
        by molecule.
        Splits the filename by 'min', such that X-min-out.maegz is written 
        as a multi-entry file to X-min.maegz.
        
        NOTE: This process RENAMES the ligand entry in the minimization output
        file, to the same name as the receptor. So, the ligand in these split
        files CANNOT be directly selected for grid generation.
        
        Parameters
        ----------
        x : str
            Filepath to output of minimization run. Must contain the characters
            'min' for filename manipulations.
        
        """
        subprocess.run(["/home/data/Schrodinger2019_1/run",
                        "/home/data/mashas/tamir/T2R14_gmeiner/automation/split_by_mol.py",
                        self.workdir+os.path.basename(x).split('min')[0]+'min-out.maegz',
                        self.workdir+os.path.basename(x).split('min')[0]+'min.maegz'],
                       cwd=self.workdir)
        
        
    def mins_to_members(self):
        """
        Split the minimization job output files into protein & ligand entires.
        Then transfer the list of split files to self.members and tests whether
        each output file can be globbed.
        """
        
        self.pool.clear()
        self.pool.restart()
        output = self.pool.map(self.pool_split_by_mol, self.mins)
        self.pool.close()
        
        self.members = [self.workdir+os.path.basename(x).split('min')[0]+"min.maegz" for x in self.mins]
        
        for x in self.members:
            if (glob.glob(x) == ''):
                self.members.remove(x)
        
        print("Loaded "+str(len(self.members))+" screens from minimization.")
    
    
    def pool_sc(self, func_input):
        """
        Parallel worker function to run ScreenConstructor methods.
        
        Parameters
        ----------
        func_input : list
            Input parameters for the ScreenConstructor method:
            Element 0: File path to a .maegz for which the screen will be constructed.
            Element 1: String to identify the ligand for model loading. Can also be "noligand".
            Element 2: File path to screening library file.
        
        """
        
        sc = ScreenConstructor(func_input[0], func_input[1], func_input[2])
        grid_command = sc.gridgen()
        dock_command = sc.dock()
        
        return [grid_command, dock_command]
        

    def dock_from_grids(self, func_input):
        """
        A parallel utility to write dock.in files using pregenerated grid .zip files instead of 
        .maegz docking output files.
        
        Parameters
        ----------
        func_input : list
            Input parameters for the utility:
            Element 0: File path to a .maegz for which the screen will be constructed.
            Element 1: String to identify the ligand for model loading. Can also be "noligand".
            Element 2: File path to screening library file.
        
        """
        
        options = ["GRIDFILE\t", "./" + os.path.basename(func_input[0]),
                   "\nLIGANDFILE\t", "./" + os.path.basename(func_input[2]),
                   "\nPRECISION\tSP",
                   "\n"]
        
        dock_name = self.workdir + os.path.basename(func_input[0]).split('grid')[0] + "dock.in"
        
        with open(dock_name, 'w') as f:
            f.write(''.join(options))
    
        return [func_input[0], dock_name]
    
        
    def construct_screen(self, ligand):
        """
        Runs the ScreenConstructor.gridgen() and .dock() methods on each member of the generation.
        Alternatively, if the first member of self.members is a .zip file, this method will
        write dock.in files using the members as pre-generated grids.
        
        NOTE: If using split files from self.mins_to_members, the ligand titles will be renamed.
        Must use 'noligand' to allow automatic selection of the ligand.
        
        Parameters
        ----------
        ligand : str
            Title of the ligand to use for grid centering.
        
        """
        func_input = [[x, ligand, self.library] for x in self.members]
        
        if (os.path.splitext(self.members[0])[-1] == '.zip'):
            self.pool.clear()
            self.pool.restart()
            output = self.pool.map(self.dock_from_grids, func_input)
            self.pool.close()
        else:    
            self.pool.clear()
            self.pool.restart()
            output = self.pool.map(self.pool_sc, func_input)
            self.pool.close()
        
        self.grids = [x[0] for x in output]
        self.docks = [x[1] for x in output]
        

    def get_used_licenses(self, label):
        """
        Query the Schrodinger license server to get the number of currently occupied
        Glide licenses.
        
        """
        checklic = subprocess.run(["/home/data/Schrodinger2019_1/licadmin", "STAT"], text=True, capture_output=True)
        for l in checklic.stdout.splitlines():
            if label in l:
                print(l)
                used_licenses = l.split('licenses')[1].split()[-1]
        return int(used_licenses)
    
    
    def run_glide(self, input_files):
        """
        Submit the provided grid.in or dock.in files to the Schrodinger job control
        queueing system.

        Configured to use the "dragon_serial_mashaq" queue. Will submit jobs if fewer
        than 100 of the "GLIDE_SUITE_22JUN2017" licenses are in use.

        Parameters
        ----------
        dockuns : list
            File paths for grid.in or dock.in files
            
        """     
        
        jobcodes = []
        
        # Note that the label must include the 'Users of' prefix
        # to avoid matching subsequent entries for the same category
        
        glide_label = 'Users of GLIDE_SUITE_17JUN2019'
        #glide_label = 'Users of GLIDE_SUITE_22JUN2017'
        #glide_label = 'Users of GLIDE_SUITE_05JUN2019'

        for f in input_files:
            submitted = False
            while not submitted:
                if (self.get_used_licenses(glide_label) < 350):
                    jobid = subprocess.run(["/home/data/Schrodinger2019_1/glide",
                                            "-HOST",
                                            "dragon_serial_mashaq",
                                            f], cwd = self.workdir, text=True, capture_output=True)
                    jobcodes.append(jobid.stdout.split()[-1])
                    submitted = True
                    time.sleep(5)
                else:
                    #print("Number of submitted jobcodes: " + str(len(jobcodes)))
                    time.sleep(1*60) # 1min sleep

        return jobcodes
    
    
    def run_grids(self):
        """
        Perform the grid generation runs.
        
        """
        if self.grids:
            self.grids_jobcodes = self.run_glide(self.grids)
        else:
            print("No known grid.in files - no grids generated.")
        
    def run_docking(self):
        """
        Perform the docking runs.
        
        """
        if self.docks:
            self.docks_jobcodes = self.run_glide(self.docks)
        else:
            print("No known dock.in files - no docking runs started.")
    
    
    def load_docking_results(self):
        """
        Generate the filepaths for the docking results, and test if each result file actually exists.
        
        """
        self.screens = [x.split('dock')[0]+"dock_pv.maegz" for x in self.docks]
        
        for x in self.screens:
            if (glob.glob(x) == ''):
                self.screens.remove(x)
        
        print("Loaded "+str(len(self.screens))+" screens from docking.")
        
    def pool_hitlist_writer(self, x):
        """
        A parallel utility to write the hitlist file for the provided screen.
        
        Parameters
        ----------
        screen : str
            Filepath for a docking result file.
            
        """
        se = ScreenEvaluator(x, 'noligand')
        se.write_hitlist()
    
    def analyse_enrichment(self, name):
        """
        Analyse the enrichment metrics for the provided screens, and save the collected
        results to the provided filename.
        
        Parameters
        ----------
        name : str
            Filename for the enrichment metrics CSV file. 
        """
        
        data = []
        
        # First, write hitlists in parallel
        self.pool.clear()
        self.pool.restart()
        output = self.pool.map(self.pool_hitlist_writer, self.screens)
        self.pool.close()
        
        for s in self.screens:
            se = ScreenEvaluator(s, 'nomodel')
            se.calc_enrichment(self.n_actives, self.n_decoys)
            se.precision_recall(se.roclist)
            df = pd.DataFrame(se.metrics, index=[0])
            df['screen'] = s
            data.append(df)
            
        self.metrics = pd.concat(data)
        self.metrics.to_csv(self.workdir+name)
        
    
    def load_reference_screens(self, refs):
        """
        Load the enrichment metric data for the given screens.
        
        Parameters
        ----------
        refs : list
            Filepaths for reference screens to load.
        """
        self.refs = refs
        data = []
        
#         # First, write hitlists in parallel
#         self.pool.clear()
#         self.pool.restart()
#         output = self.pool.map(self.pool_hitlist_writer, self.screens)
#         self.pool.close()
        
        for s in self.refs:
            se = ScreenEvaluator(s, 'nomodel')
            se.write_hitlist()
            se.calc_enrichment(self.n_actives, self.n_decoys)
            df = pd.DataFrame(se.metrics, index=[0])
            df['screen'] = s
            data.append(df)
            
        self.refs = pd.concat(data)
        
    
    def metric_distance(self, row):
        """
        Calculate the distance metric from maximum enrichment.
        
        """
        auc_dist = (self.metrics['AUC'].max() - row['AUC'])/self.metrics['AUC'].max()
        ef10_dist = (self.metrics['EF10'].max() - row['EF10'])/self.metrics['EF10'].max()
        return auc_dist + ef10_dist
    
    
    def fetch_top_screens(self, n):
        """
        Return the names of the n-top enriching screens, as measured
        by the distance metric calculated by "self.metric_distance".
        
        Parameters
        ----------
        n : int
            The number of top screens to fetch.
            
        """
        self.metrics['metric_distance'] = self.metrics.apply(lambda row: self.metric_distance(row), axis=1)
        self.topmodels = self.metrics.nsmallest(n, 'metric_distance')
        return self.topmodels['screen'].tolist()
    
    
    def save_convergence_poses(self, name, screens, ligands):
        """
        Saves a file containing the receptors and specified ligands from the given screens.
        The output file is titled 'name_convergence.maegz'.

        Parameters
        ----------
        name : str
            A label to be prefixed to the saved convergence pose file.
        screens : list
            List of screens to merge and filter.
        ligands : list
            List of title for ligands selected to save out poses for.
        
        """
        
        with open(self.workdir+"screens.merge", 'w') as f:
            for s in screens:
                f.write(s+"\n")
        f.close()    

        subprocess.run(["/home/data/Schrodinger2019_1/utilities/glide_merge",
                        "-o",
                        "merge.maegz",
                        "-epv",
                        "-f",
                        "screens.merge"],
                       cwd = self.workdir)

        subprocess.run(["/home/data/Schrodinger2019_1/utilities/proplister",
                        "-c",
                        "-p",
                        "s_m_title",
                        "merge.maegz",
                        "-o",
                        "filter_titles.csv"],
                       cwd = self.workdir)

        df = pd.read_csv(self.workdir+"filter_titles.csv")
        
        # Detect receptor entries - once minimized, these are set to '0'
        #receptors = [x for x in df['s_m_title'].tolist() if ':' in x]
        receptors = ['0']
        retrieve = receptors+ligands

        with open(self.workdir+"filter.titles", 'w') as f:
            for t in retrieve:
                f.write(t+"\n")
        f.close()

        subprocess.run(["/home/data/Schrodinger2019_1/utilities/maesubset",
                        "-title_list_file",
                        "filter.titles",
                        "merge.maegz",
                        "-o",
                        name+"convergence.maegz"],
                       cwd = self.workdir)



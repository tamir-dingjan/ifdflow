#!/usr/bin/env python
# coding: utf-8

class Model:
    """The protein-ligand complex coordinates, with a minimize method.
        
    Attributes
    ----------
    pose : mdtraj.Trajectory object 
    
    """
    def __init__(self, filename, ligand):
        """
        Initialise the Model with an MDTraj Trajectory object created
        from the given Maestro file.
        
        The input file must contain multiple entries, with the first entry
        containing just the protein and subsequent entries containing ligands.
        This is the default structure for docking output .maegz files.
        
        The protein-only entry is extracted using 'maesubset', and is used in 
        grid generation.
        
        The MDTraj Trajectory object, however, must contain a ligand to enable
        centroid position calculations. So, the complex must be assembled before
        loading into MDTraj - this is accomplished by 'structcat'. To do this,
        the desired ligand file has to be extracted from the docking output .maegz.
        This file is identified by the Maestro Title string, extracted using 
        'maesubset' and merged. Currently, the ligand chosen is just the second 
        entry.
        
        Parameters
        ----------
        filename : str
            Filename for the ligand-protein coordinates. The file must
            be a Maestro format .maegz file containing a single protein molecule
            as the first entry, and any number of ligand molecules in 
            subsequent entries.
        ligand : str
            Title of the ligand to be used for assembling the complex.
            This ligand will feature in the model minimization, and so
            is an important choice as it will guide the model architecture
            during refinement. The string must correspond to the Maestro 
            title of the ligand entry in the input pose-viewer file. This field
            is termed 's_m_title' in Maestro's internal referencing available
            with the proplister tool. If the string 'noligand' is provided, 
            the model will be setup with the first ligand listed in the docking 
            input file. 

            
        """
        
        # Cast the provided filename as an absolute path to enable processing in 
        # the same directory
        self.input_file = os.path.realpath(filename)
        self.name = os.path.splitext(os.path.basename(self.input_file))[0]
        self.output_dir = os.path.dirname(self.input_file) + "/"
        self.load_status = False
        
        # The Model is useless without a ligand. So, check if the specified ligand
        # is available within the input.
        # However, when assembling a new model (i.e., not one from docking output)
        # a complex with the desired ligand may not be available. In this case,
        # use the 'noligand' title as a manual override of sorts to force
        # model loading.
        
        self.ligand = ligand
        
        if (self.ligand == 'noligand'):
            
            # If the self.ligand attribute is set to 'noligand', then this model
            # does not come from a docking run containing the desired convergence
            # ligand. Therefore, all we really need from the ligand is the know
            # the binding site location. The first ligand in the file will do.
            # Use the proplister function to get this information, and update
            # the self.ligand attribute to have the new ligand name.


            entry_titles = subprocess.run(["/home/data/Schrodinger2018_3/utilities/proplister",
                                           "-p",
                                           "s_m_title",
                                           self.input_file], text=True, capture_output=True)
            self.ligand = entry_titles.stdout.rstrip("\n").split('\n')[3].strip()
            
            

        checkligand = subprocess.check_output(["/home/data/Schrodinger2018_3/utilities/maesubset",
                                                 "-title",
                                                 self.ligand,               
                                                 self.input_file],
                                                cwd = os.path.dirname(self.input_file))

#         checkligand_p = subprocess.Popen(["/home/data/Schrodinger2018_3/utilities/maesubset",
#                                                  "-title",
#                                                  self.ligand,               
#                                                  self.input_file],
#                                          stdout = checkligand,
#                                          cwd = os.path.dirname(self.input_file))
#         checkligand_p.wait()
    
        if (checkligand == ''):
            self.ligandmatch = False
        else:
            self.ligandmatch = True

        
    def load(self):
        """
        Run the Model's file-conversion routines and load the MDTraj object.
        
        This will perform multiple operations:
            1. The protein structure is extracted to <self.name>-protein.mae
            2. All entries matching the provided ligand name are saved to 
               <self.name>-matchedligands.mae.
            3. The ligand structure is extracted to <self.name>-ligand.mae
            4. The protein-ligand complex is assembled to <self.name>-merge.mae.
               This file contains two entries in the Maestro file.
            5. The protein-ligand complex file is merged to a single Maestro entry,
               saved to <self.name>-complex.mae.
            6. The merged complex is saved to PDB to allow reading by MDTraj, saved
               to <self.name>-complex.pdb.
        
        The resulting MDTraj Trajectory object is stored in self.pose.
        
        """
        
        # Extract the protein from the pose-viewer format multi-entry .maegz
        # ASSUMPTION: This step assumes that the first entry in the pose-viewer file
        #             is the receptor. Usually, this will be the case.
        
        if (self.ligandmatch == False):
            print("No ligand match for this model - cancel loading.")
            return
        
        
        with open(self.output_dir + self.name + "-protein.mae", 'w') as f:
            extractprotein = subprocess.Popen(["/home/data/Schrodinger2018_3/utilities/maesubset",
                                                 "-n",
                                                 "1",               
                                                 self.input_file],
                                                stdout = f,
                                                cwd = os.path.dirname(self.input_file))
        

                
        
        with open(self.output_dir + self.name + "-matchedligands.mae", 'w') as f:
            fetchligands = subprocess.Popen(["/home/data/Schrodinger2018_3/utilities/maesubset",
                                               "-title",
                                               self.ligand,               
                                               self.input_file],
                                              stdout = f,
                                              cwd = os.path.dirname(self.input_file))
        
        # With no delay, the script doesn't seem able to read the -matchedligands.mae file.
        # Subsequent execution with the above block commented out works as normal.
        # Perhaps it just takes a little bit of time to write the file?
        # Try using a 5 sec delay. Works. 1 second also works, but sometimes fails.
        #time.sleep(3)
        
        # Alternatively, try waiting for the process to terminate
        fetchligands.wait()
    
        with open(self.output_dir + self.name + "-ligand.mae", 'w') as f:
            extractligand = subprocess.Popen(["/home/data/Schrodinger2018_3/utilities/maesubset",
                                                "-n",
                                                "1",
                                                self.output_dir + self.name + "-matchedligands.mae"],
                                               stdout = f,
                                               cwd = os.path.dirname(self.input_file))

        merge = subprocess.Popen(["/home/data/Schrodinger2018_3/utilities/structcat",
                                    "-imae",
                                    self.output_dir + self.name + "-protein.mae",
                                    "-imae",
                                    self.output_dir + self.name + "-ligand.mae",
                                    "-omae",
                                    self.output_dir + self.name + "-complex_pv.mae"])
        
        complexation = subprocess.call(["/home/data/Schrodinger2018_3/run",
                                          "pv_convert.py",
                                          "-mode",
                                          "merge",
                                          self.output_dir + self.name + "-complex_pv.mae"])
    
        writepdb = subprocess.call(["/usr/people/mashagr/tamird/bin/silico1.14/bin/mol_combine",
                                      self.output_dir + self.name + "-protein.mae",
                                      self.output_dir + self.name + "-ligand.mae",
                                      "-O",
                                      self.output_dir + self.name + "-complex.pdb",
                                      "-o",
                                      "pdb",
                                      "-dr"])
        
        
        try:
            self.pose = md.load(self.output_dir + self.name + "-complex.pdb")
        except IOError:
            print("Could not load file: " + self.output_dir + self.name + "-complex.pdb")
            
        self.complexfile = self.output_dir + self.name + "-complex-out_complex.mae"
        self.protfile = self.output_dir + self.name + "-protein.mae"
        self.load_status = True
        
    
    def get_ligand_resid(self):
        """Use the MDTraj Trajectory object to automatically identify the ligand
        residue name. This process relies on the ligand consisting of a single
        molecule possessing a single residue name. Single-atom molecules are 
        explicitely ignored for this step, as some protein hydrogen atoms are
        not always correctly identified as protein. The ligand residue is then 
        selected as the residue possessing the smallest number of atoms.
        
        Returns
        -------
        str
            A string containing the residue name of the ligand
            
        """
        molecules = list(nx.connected_components(self.pose.topology.to_bondgraph()))
        molecules.sort(key=len)
        largemols = [m for m in molecules if (len(m)>1)]
        ligand_resid = str(list(set([atom.residue.resSeq for atom in largemols[0]]))[0])
        return ligand_resid
    
    def get_molecules(self):
        """Return the molecules detected by MDTraj.
        
        Returns
        -------
        list
            The length of the list corresponds to the number of molecules 
            detected. Each element of the list is the number of atoms in each 
            molecule
        """
        molecules = list(nx.connected_components(self.pose.topology.to_bondgraph()))
        molecules.sort(key=len)
        return molecules
    
    def get_ligand_centroid(self):
        """Use the MDTraj Trajectory object to find the centroid of the ligand.
        
        Returns
        -------
        list
            A three-element list containing the X, Y, and Z coordinates of the 
            position corresponding to the centroid of the ligand atoms.
        
        """
        ligand_atoms = self.pose.topology.select("resSeq "+self.get_ligand_resid())
        centroid = np.mean(self.pose.xyz[0][ligand_atoms], axis=0)*10
        return list(centroid)
    
    def minimize(self):
        """Write the Schrodinger command to minimize the ligand-protein 
        complex contained in the model object. 
        
        This method tries to automatically identify the ligand. For this, 
        the PDB-export version of the original Maestro file is used,
        as MDTraj cannot read Maestro formats.
        
        
        Parameters
        ----------
        name : str
            The name of the Maestro cordinate file to be used, excluding
            the file path and extension. The minimization command will use
            this file, and the options file is also named using this string.
            
        Returns
        -------
        string
            A shell command to use the Schrodinger suit to minimize the 
            model. This command can be executed any time, using the coordinate
            file written. The command references an options file, which is
            written to the disk using the name parameter
            
        """
        
        # Write the minimization job input file
        # This can be executed using the following in the terminal:
        # >>$SCHRODINGER/prime jobname
           
        # Identify the ligand
        ligand_resid = self.get_ligand_resid()
        
        options = ["STRUCT_FILE\t", os.path.basename(self.complexfile),
                   "\nJOB_TYPE\tREFINE",
                   "\nPRIME_TYPE\tSITE_OPT",
                   "\nSELECT\tasl = fillres within 5.000000 ( res.num ", ligand_resid, ")",
                   "\nLIGAND\tasl = res.num ", ligand_resid,
                   "\nNPASSES\t1",
                   "\nINCRE_PACKING\tno",
                   "\nUSE_CRYSTAL_SYMMETRY\tno",
                   "\nUSE_RANDOM_SEED\tno",
                   "\nSEED\t0",
                   "\nOPLS_VERSION\tOPLS3e",
                   "\nEXT_DIEL\t80.00",
                   "\nUSE_MEMBRANE\tno"]
        
        jobname = self.name + "-min"
        
        with open(self.output_dir + jobname + ".inp", 'w') as f:
            f.write(''.join(options))

        #command = "/home/data/Schrodinger2018_3/prime " + jobname
        
        return self.output_dir + jobname + ".inp"
    



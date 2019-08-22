#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import glob
import subprocess
import time
import math
import shutil
import pathos.pools as pp
import numpy as np
import pandas as pd
import mdtraj as md
import networkx as nx
from sklearn import metrics 
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns


# # Notes

# This notebook will contain a class and methods that enable the automatic refinement of homology models for virtual screening. The process will require the following:
# 1. A retrospective screening library containing known actives and decoys.
# 2. Selection criteria for the pose to converge upon, consisting of a ligand identity and the number of top-ranked models required to converge. When this pose is formed by the requested number of the top-ranked models, the refinement ends. 
# 3. Binding site location.
# 
# The following methods will be defined:
# 1. Screening - to conduct the screen of the specified library at a model in the selected binding site location.
# 2. Minimization - to minimize the binding site residues around the specified convergence pose.
# 3. Enrichment - to calculate the enrichment metrics for a completed screen.
# 4. Convergence - to detect whether the top-ranked models from a generation have converged.
# 5. Refine - a wrapper method to perform the minimization and screening methods at the generation level, since the entire generation is refined simultaneously.
# 
# The classes used to structure this process:
# 
# "Generation" - a collection of "Screens". The Convergence method operates at the "Generation" level.
# 
# "Screen" - a single screening result, consisting of a protein with a list of hits retrieved by the protein.  Each screening result will have its enrichment metrics associated with it, so that the top-ranked models can be easily selected. The "Screening" and "Enrichment" methods operate at the "Screen" level. In particular, the "Screening" method populates a "Screen" by running a virtual screening calculation on the "Screen" object's protein model.
# 
# "Model" - contains a single protein conformation. This object will have the "Minimization" method as well as methods for exporting the protein coordinates.
# 

# # Thoughts

# Do I need to build a job monitor to make this work? Some kind of process that can schedule jobs and report the results back when the jobs finish? Otherwise, the whole script could be stuck waiting for a single docking job to finish.
# 
# Definitely need a way to schedule multiple jobs simultaneously. Could do this with a mix of manual and automated - for example, could write out a big file of Schrodinger terminal commands and then just manually give that to the Schrodinger server in Bash. Then, could re-run this notebook on the output of those jobs when they finish - that way, the notebook does the heavy lifting of checking for convergence and writing job submission files, but doesn't have to be aware of when the jobs finish. An added bonus, the notebook won't be able to automatically submit jobs so things won't get out of hand by accident.
# 
# This means that a Screen object should be able to be created from a .maegz file, for easy input of the completed jobs when they finish. Simple way to do this is to have two kinds of Screen object - a ScreenConstructor, and a ScreenEvaluator

# I'd like to be able to do the convergence check within MDTraj, without needing to run the RMSD calculations through Maestro. But, this might cause problems with the binding site location when starting a new generation, because MDTraj and Maestro might shift around coordinate spaces / residue numbering. Watch out for this. So far the coordinate space seems to remain consistent.

# File Nomenclature:
# 
# Models are named by Generation, Screen, and Model: "0-1-2.pdb"
# 
# Screens are named by Generation and Screen: "0-1_{grid/dock}.maegz"

# In[ ]:





# In[2]:


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
    


# In[35]:


testfile = "/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/icb_1-out_pv.mae"
x = Model(testfile, "noligand")


# In[36]:


x.ligand


# In[27]:


entry_titles = subprocess.run(["/home/data/Schrodinger2018_3/utilities/proplister",
                               "-p",
                               "s_m_title",
                               testfile], text=True, capture_output=True)


# In[28]:


type(entry_titles.stdout)


# In[29]:


entry_titles.stdout


# In[31]:


entry_titles.stdout.rstrip("\n").split('\n')[3].strip()


# In[ ]:


self.ligand = entry_titles.rstrip("\n").split('\n')[3].strip()
            


# In[ ]:





# In[3]:


class Grid:
    def __init__(self, filename):
        """
        In some cases, protein structure files may be generated without
        any bound or docked ligands. Generating a grid from these structures
        requires identifying a binding site independantly of any bound 
        ligand. 
        
        This can be done by specifying a set of residues defining the borders
        of the site and defining an equidistant center.
        
        Parameters
        ----------
        filename : str
            PDB file containing a single protein without any ligands. 
        """
        
        self.input_file = os.path.realpath(filename)
        self.name = os.path.splitext(self.input_file)[0].split('/')[-1]
        self.output_dir = os.path.dirname(self.input_file) + "/"
        self.maename = self.name.split('Q9NYV8.B99990')[-1]
        
        
        try:
            self.t = md.load(self.input_file)
        except IOError:
            print("Could not load file: " + self.input_file)
        
        
    def find_site_center(self, selected_residues):
        """
        Identify the center of the provided residues.
        
        Parameters
        ----------
        selected_residues : lst
            Collection of residue index numbers drawn from 
            the topology of the file.
        """
        
        try:
            self.residue_centers = [np.mean(self.t.xyz[0][self.t.topology.select("resSeq "+str(i))], axis=0) for i in selected_residues]
            self.site_center = np.mean(np.asarray(self.residue_centers), axis=0) * 10
        except IndexError:
            print("Unable to select residue "+str(i)+". Site not found.")
        
    
    def convert_to_mae(self):
        """
        Convert the input file to a MAE format file.
        """
        
        convert = subprocess.call(["/home/data/Schrodinger2019_1/utilities/structconvert",
                                     "-ipdb",
                                     self.input_file,
                                     "-omae",
                                     self.output_dir + self.name + ".mae"])
    
    
    def generate_grid(self):
        """
        Write the grid.in file for running the grid generation job.
        """
        
        options = ["GRID_CENTER\t", str(self.site_center[0]), ",", str(self.site_center[1]), ",", str(self.site_center[2]), ",",
                   "\nGRIDFILE\t", self.name + "-grid.zip",
                   "\nINNERBOX\t10, 10, 10",
                   "\nOUTERBOX\t25, 25, 25",
                   "\nRECEP_FILE\t", "out_" + str(self.maename) + ".mae", 
                   "\n"]

        grid_file = self.output_dir + self.name + "-grid.in"
        
        with open(grid_file, 'w') as f:
            f.write(''.join(options))
        
        return grid_file
        


# In[177]:


g = Grid("/home/data/mashas/tamir/T2R14_gmeiner/antagonists/inactive_state_homology/models/multi_align_2/Q9NYV8.B99990001.pdb")


# In[178]:


g.find_site_center([69, 89, 175])


# In[179]:


g.residue_centers


# In[180]:


g.site_center


# In[87]:


g.convert_to_mae()


# In[88]:


g.generate_grid()


# In[181]:


models = glob.glob("/home/data/mashas/tamir/T2R14_gmeiner/antagonists/inactive_state_homology/models/multi_align_2/*pdb")


# In[182]:


for m in models:
    g = Grid(m)
    g.find_site_center([89, 69, 175])
    g.generate_grid()


# In[92]:


gridjobs = glob.glob("/home/data/mashas/tamir/T2R14_gmeiner/antagonists/inactive_state_homology/screening/multi_align_2/*grid*in")


# In[171]:


t = md.load("/home/data/mashas/tamir/T2R14_gmeiner/antagonists/inactive_state_homology/models/multi_align_2/Q9NYV8.B99990001.pdb")


# In[175]:


for i in t.topology.select("resSeq 89 69 175"):
    print t.topology.atom(i).residue


# In[ ]:





# In[3]:


class IFDProcessor:
    def __init__(self, filename):
        """Given a Maestro file containing the results of an Induced Fit Docking calculation
        (i.e., multiple entries each containing a protein and ligand), produce a "docking output"
        style of Maestro file (i.e., containing a single protein molecule as the first entry,
        followed by another entry for the ligand).
        
        Parameters
        ----------
        name : str
            Filename for the IFD results file. The file must
            be a Maestro format .maegz file containing
            a single entry, consisting of a ligand-protein complex.
            
        """
        
        self.input_file = os.path.realpath(filename)
        self.output_file = "No output yet"
        
        
    def run(self):
        """Execute the conversion of the IFD docking results file into a Maestro Pose Viewer
        file.
        
        """
        
        # Split the IFD file into protein and ligand
        # Use subprocess.Popen to allow specifying of the directory
        split = subprocess.Popen(["/home/data/Schrodinger2019_1/run",
                                    "pv_convert.py",
                                    "-mode",
                                    "split_pv",
                                    self.input_file,
                                    "-lig_last_mol"],
                                   cwd=os.path.dirname(self.input_file))
        
        output_name = os.path.splitext(os.path.basename(self.input_file))[0] + "-out_pv" + os.path.splitext(os.path.basename(self.input_file))[1]
        self.output_file = os.path.dirname(self.input_file) + "/" + output_name
        
        
        


# In[12]:


ifds = glob.glob("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/*mae")


# In[13]:


len(ifds)


# In[18]:


x = IFDProcessor(ifds[0])


# In[19]:


x.run()


# In[20]:


x.output_file


# In[21]:


for i in ifds:
    x = IFDProcessor(i)
    x.run()


# In[ ]:





# In[4]:


class ScreenConstructor:
    def __init__(self, model_path, ligand, library_path):
        """Initialise the ScreenConstructor with the paths to the ligand-protein
        model and retrospective screening library files. The model path is 
        used to add a Model object.
        
        Parameters
        ----------
        model_name : str
            Name of the docking output .maegz file, beginning with an entry 
            consisting of a single molecule protein and followed by at least one
            entry consisting of a single molecule ligand
        ligand : str
            Label for the ligand used to load the model.
        library_path : str
            Path to the retrospective screening library, containing the actives
            and decoys to be screened.
        
        """
        self.model = Model(model_path, ligand)
        self.library_path = os.path.realpath(library_path)
        self.output_dir = os.path.dirname(os.path.realpath(model_path)) + "/"        
    
    def gridgen(self):
        """
        Write the Schrodinger input files and return the name of the grid.in file.
        
        """
        # Load the Model's coordinates - includes a check for ligand presence
        if (self.model.load_status == False):
            self.model.load()
        
        # Identify the ligand centroid
        ligand_centroid = self.model.get_ligand_centroid()
        
        # Grid generation options      
        options = ["GRID_CENTER\t", str(ligand_centroid[0]), ",", str(ligand_centroid[1]), ",", str(ligand_centroid[2]), ",",
                   "\nGRIDFILE\t", self.model.name + "-grid.zip",
                   "\nINNERBOX\t10, 10, 10",
                   "\nOUTERBOX\t25, 25, 25",
                   "\nRECEP_FILE\t", os.path.basename(self.model.protfile), 
                   "\n"]

        grid_file = self.output_dir + self.model.name + "-grid.in"
        
        with open(grid_file, 'w') as f:
            f.write(''.join(options))
        
        return grid_file
    
    def dock(self):
        """Write the Schrodinger input files and return the execution command to perform
        screening, using the created grid and the specified library file.
        
        """
        options = ["GRIDFILE\t", self.model.name + "-grid.zip",
                   "\nLIGANDFILE\t", os.path.basename(self.library_path),
                   "\nPRECISION\tSP",
                   "\n"]
        
        dock_file = self.output_dir + self.model.name + "-dock.in"

        with open(dock_file, 'w') as f:
            f.write(''.join(options))
        
        return dock_file


# In[5]:


class ScreenEvaluator:
    def __init__(self, model_path, ligand):
        """Initialise the ScreenEvaluator with the screening results file. This populates 
        the self.results object. The Model object is not populated upon initialization, to
        allow rapid inspection of ScreenEvaluator objects.
        
        The self.results object contains the hitlist from the screening experiment as a 
        Pandas dataframe, specifying the Maestro Title, Docking Score, and Decoy status
        of each hit.
        
        Parameters
        ----------
        model_path : str
            Name of the docking output .maegz file, beginning with an entry 
            consisting of a single molecule protein and followed by at least one
            entry consisting of a single molecule ligand
        ligand : str
            Name of the ligand used to load the model. This is passed directly to
            the Model() class method. Options include any string corresponding
            to a ligand title, and 'noligand', which triggers automatic identification
            of the ligand. Alternatively, to bypass model loading, can pass 'nomodel'
        
        """
        # Extract the project table results from the docking output file
        self.input_file = os.path.realpath(model_path)
        self.name = os.path.splitext(os.path.basename(self.input_file))[0]
        self.output_dir = os.path.dirname(self.input_file) + "/"
        
        # Model loading is a slow step, and not required if a hitlist
        # has already been written (usually)
        if not (ligand == 'nomodel'):
            self.model = Model(model_path, ligand)
        
        self.results = ''
        self.metrics = {}
    
    
    def write_hitlist(self):
        """Write the docking hitlist to disk. The hitlist is saved as a .csv
        with the following fields:
        s_m_title : The Maestro title field for each entry.
        r_i_docking_score : Docking score for each hit.
        s_user_RetroScreenRole : type of compound (active/decoy)
        
        This data is loaded into self.results as a Pandas dataframe.
        
        """
        
        with open(self.output_dir + self.name + "-results.csv", 'w') as f:
            extractprotein = subprocess.call(["/home/data/Schrodinger2019_1/utilities/proplister",
                                                "-c",
                                                "-p",
                                                "s_m_title",               
                                                "-p",
                                                "r_i_docking_score",
                                                "-p",
                                                "s_user_RetroScreenRole",
                                                "-p",
                                                "s_user_CompoundSet",
                                                self.input_file],
                                               stdout = f)
                
            
    def compare_score(self, active_score, decoy_score):
        """Function used in AUC calculation to build a score tally based on
        the relationship between an active score and a decoy score. 
        
        Currently, this function is set to assume that lower scores are better.
        Accordingly, the AUC is increased when an active has a lower score 
        than a decoy.
        
        Parameters
        ----------
        active_score : float
            A docking score for an active compound.
        decoy_score : float
            A docking score for a decoy compound.
            
        """
        if active_score < decoy_score:
            return float(1)
        elif active_score > decoy_score:
            return float(0)
        elif active_score == decoy_score:
            return float(0.5)
        
        
    def calc_auc(self, actives_scores, decoys_scores, n_actives, n_decoys):
        """Calculate the ROC AUC for a set of actives and decoys.
        
        This routine scales the AUC number depending on how many of the actives
        and decoys present in the input are present in the score lists.
        The user provides the number of actives and decoys originally used
        as input to the screen, and the number of items in the actives_scores
        and decoys_scores lists are used to determine the scaling of the AUC.
        
        Parameters
        ----------
        actives_scores : list
            List of docking scores for actives
        decoys_scores : list
            List of docking scores for decoys
        n_actives : int
            Number of active compounds used as input to this screen
        n_decoys : int
            Number of decoy compounds used as input to this screen
            
        Returns
        -------
        auc : float
            A scaled AUC that reflects the probability that an active will
            score better than a decoy.
            
        """
        
        n_a_found = float(len(actives_scores))
        n_d_found = float(len(decoys_scores))
        n_a_input = float(n_actives)
        n_d_input = float(n_decoys)

        score_compare_tally = float(0)

        for a in actives_scores:
            for d in decoys_scores:
                score_compare_tally += self.compare_score(a, d)

        auc_raw = float(score_compare_tally / (n_a_found * n_d_found))

        auc_scaled = float((auc_raw * (n_a_found / n_a_input) * (n_d_found / n_d_input)) + ((1 - (n_d_found/n_d_input)) * (n_a_found/n_a_input)))
        
        return auc_scaled

    
    def calc_ef_and_sd(self, threshold, n_actives, n_decoys):
        """Calculate the enrichment factor at the selected threshold.
        The EF is calculated as the ratio between an arbitrarily-selected
        fraction of inactives and the fraction of actives retrieved at that
        fraction of inactives. This definition also corresponds to the 
        gradient of a line constructured from the origin to the ROC-curve value
        coordinate at the selected inactive fraction.
        
        Parameters
        ----------
        threshold : float
            The fraction of inactives at which to determine the EF.
            I.e., EF1% corresponds to a threshold of 0.01.
        n_actives : int
            Number of actives used as input to the screen
        n_decoys : int
            Number of decoys used as input to the screen.
        
        Returns
        -------
        ef : float
            The enrichment factor.
        sd : float
            The analytical error of the enrichment factor
            
        """
        
        for (i, a) in zip(self.roclist_x, self.roclist_y):
            if i > threshold:
                ef = float(a/i)
                break
        
        s = (a*np.log(a)) / (i*np.log(i))
        sd = (((a*(1-a))/float(n_actives)) + ((s**2)*i*(1-i))/n_decoys) / (i**2)
        
        return (ef, sd)
    
    
    def calc_ci95(self, variance):
        """Calculates a 95% confidence interval for any metric, with the 
        assumption that the distribution is normal. 
        
        Caveat from the original paper:
        "This will break down for systems with very large errors (generally
        due to having very few actives or decoys) ... When large errors do
        arise the exact magnitude of the error is generally less relevant than
        the fact that the error is large"
        
        Parameters
        ----------
        variance : float
            The variance of any metric
        
        Returns
        -------
        ci95 : float
            The 95% confidence interval. Adding and subtracting this value 
            from the mean value of the metric will give the values that define
            the borders of 95% confidence for that value.
            
        """
        ci95 = 1.96 * variance**0.5
        return ci95
    
    
    def precision_recall(self, hitlist):
        """Calculate the precision-recall profile and associated summary statistics
        (F1 score, Matthews correlation quality, Youden's J, PR-AUC).
        
        The summary statistics are defined for every point along the precision-recall
        curve, so will report the maximal value.
        
        Parameters
        ----------
        hitlist : list
            A list containing 'active' and 'decoy' entries. Can be generated from the self.results.role
            Pandas sheet column.
        
        """
        
        self.precision = []
        self.recall = []
        
        self.f1 = []
        self.mcc = []
        self.youdenj = []
        
        
        # Move the threshold through the hitlist
        for i in range(1, len(hitlist)):
            toplist = hitlist[0:i]
            bottomlist = hitlist[i:]
            
            tp = 0
            fp = 0
            tn = 0
            fn = 0
        
            # At each step, calculate the precision and recall
            for x in toplist:
                if x == 'active':
                    tp += 1
                elif x == 'decoy':
                    fp += 1
            for x in bottomlist:
                if x == 'decoy':
                    tn += 1
                elif x == 'active':
                    fn += 1
            
            #print len(toplist), len(bottomlist), "TP: ", tp, "FP: ", fp, "TN: ", tn, "FN: ", fn
            
            precision = float(tp) / float(tp + fp)
            recall = float(tp) / float(tp + fn)
            
            self.precision.append(precision)
            self.recall.append(recall)
            
            # Calculate the summary statistics
            f1 = np.reciprocal(((np.reciprocal(recall)) + (np.reciprocal(precision)))/2)
            
            self.f1.append(f1)
            
            mcc = float((tp * tn) - (fp * fn)) / float(math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
            self.mcc.append(mcc)
            
            youdenj = (float(tp)/float(tp+fn)) + (float(tn)/float(tn+fp)) - 1
            self.youdenj.append(youdenj)
            
        
        # Calculate PRAUC
        self.metrics['prauc'] = metrics.auc(self.recall, self.precision)
        
        # Store the summary statistics
        self.metrics['f1'] = max(self.f1)
        self.metrics['mcc'] = max(self.mcc)
        self.metrics['youden_j'] = max(self.youdenj)
        
        # Store the score cutoff associated with the maximum Youden's J index
        self.metrics['youden_j_score'] = self.results['score'][self.youdenj.index(max(self.youdenj))]
        
    
    def calc_enrichment(self, n_actives, n_decoys):
        """Calculate enrichment metrics for this screen. The metrics are stored in the 
        self.metrics dictionary object.
        
        Parameters
        ----------
        n_actives : int
            The number of actives in the input library used to perform the screen.
        n_decoys : int
            The number of decoys in the input library used to perform the screen.
            
        """
        # Load hitlist
        try:
            self.results = pd.read_csv(self.output_dir + self.name + "-results.csv", names = ['title', 'score', 'role', 'compoundset'], skiprows = [0,1])
        except IOError:
            print("File not found: " + self.output_dir + self.name + "-results.csv")
            return
        
        # Remove duplicate items from the hitlist.
        self.results.drop_duplicates(subset='title', keep="first", inplace=True)
        
        # Extract score lists for actives and decoys
        self.actives_scores = self.results[self.results["role"] == 'active'].score.tolist()
        self.decoys_scores = self.results[self.results["role"] == 'decoy'].score.tolist()
        self.roclist = self.results.role.tolist()
        
        # Prepare ROC lists
        x = [0]
        y = [0]
        d_count = 0
        a_count = 0
        for i in self.roclist:
            if i == 'decoy':
                d_count += 1
            elif i == 'active':
                a_count += 1
            x.append(d_count)
            y.append(a_count)

        #Scale the coordinates by the total numbers
        self.roclist_x = np.asarray(x, dtype=float)/float(d_count)
        self.roclist_y = np.asarray(y, dtype=float)/float(a_count)
        
        # Populate the self.enrichment object with the metrics
        self.metrics['AUC'] = self.calc_auc(self.actives_scores, self.decoys_scores, n_actives, n_decoys)
        
        q1 = self.metrics['AUC'] / float(2-self.metrics['AUC'])
        q2 = (2*(self.metrics['AUC'])**2) / (1 + self.metrics['AUC'])
        self.metrics['AUC_SD'] = (self.metrics['AUC']*(1-self.metrics['AUC']) + ((n_actives-1)*(q1-(self.metrics['AUC'])**2)) + ((n_decoys-1)*(q2-(self.metrics['AUC'])**2))) / float(n_actives * n_decoys)
        
        self.metrics['EF1'], self.metrics['EF1_SD'] = self.calc_ef_and_sd(0.01, n_actives, n_decoys)
        self.metrics['EF2'], self.metrics['EF2_SD'] = self.calc_ef_and_sd(0.02, n_actives, n_decoys)
        self.metrics['EF5'], self.metrics['EF5_SD'] = self.calc_ef_and_sd(0.05, n_actives, n_decoys)
        self.metrics['EF10'], self.metrics['EF10_SD'] = self.calc_ef_and_sd(0.10, n_actives, n_decoys)
        
        self.metrics['AUC_CI95'] = self.calc_ci95(self.metrics['AUC_SD'])
        self.metrics['EF1_CI95'] = self.calc_ci95(self.metrics['EF1_SD'])
        self.metrics['EF2_CI95'] = self.calc_ci95(self.metrics['EF2_SD'])
        self.metrics['EF5_CI95'] = self.calc_ci95(self.metrics['EF5_SD'])
        self.metrics['EF10_CI95'] = self.calc_ci95(self.metrics['EF10_SD'])
    
    
    def plot_roc(self, color):
        """Plot the ROC curve and save to a file named <self.name>_roc.png.
        
        """

        x = [0]
        y = [0]
        d_count = 0
        a_count = 0
        for i in self.roclist:
            if i == 'decoy':
                d_count += 1
            elif i == 'active':
                a_count += 1
            x.append(d_count)
            y.append(a_count)

        #Scale the coordinates by the total numbers
        x_scale = np.asarray(x, dtype=float)/float(d_count)
        y_scale = np.asarray(y, dtype=float)/float(a_count)

        plt.figure(figsize=(3,3))

        sns.set(context="notebook", style="ticks")

        plt.plot([0,1], [0,1], color=sns.xkcd_rgb["grey"], linestyle='--')
        plt.plot(x_scale, y_scale, sns.xkcd_rgb[color])

        plt.ylim(0,1)
        plt.ylabel("Actives")

        plt.xlim(0,1)
        plt.xlabel("Inactives")

        plt.savefig(self.output_dir + self.name + "_roc.png", dpi=300, format="png", bbox_inches="tight", transparent=True)

    def plot_roc_marked(self, color, watch_compound, watch_role):
        """Plot the ROC curve and add text labels to compounds
        specified.
        
        Parameters
        ----------
        color : str
            The colour to plot the ROC curve in.
        watch_compound : str
            Text label to match against the "compoundset" property for
            selecting the marked compounds.
        watch_role : str
            Text label to match against the "role" property for selecting
            the marked compounds.
        """

        x = [0]
        y = [0]
        marked_y = []
        marked_x = []
        marked_labels = []

        d_count = 0
        a_count = 0

        for i in range(0, len(self.results)):
            if (self.results.iloc[i]['role'] == 'decoy'):
                d_count += 1
            elif (self.results.iloc[i]['role'] == 'active'):
                a_count += 1
            if (self.results.iloc[i]['compoundset'] == watch_compound) and (self.results.iloc[i]['role'] == watch_role):
                marked_y.append(a_count)
                marked_x.append(d_count)
                marked_labels.append(self.results.iloc[i]['title'])

            x.append(d_count)
            y.append(a_count)


        #Scale the coordinates by the total numbers
        x_scale = np.asarray(x, dtype=float)/float(d_count)
        y_scale = np.asarray(y, dtype=float)/float(a_count)
        marked_x_scale = np.asarray(marked_x, dtype=float)/float(d_count)
        marked_y_scale = np.asarray(marked_y, dtype=float)/float(a_count)

        plt.figure(figsize=(6,6))

        sns.set(context="notebook", style="ticks")

        plt.plot([0,1], [0,1], color=sns.xkcd_rgb["grey"], linestyle='--')
        plt.plot(x_scale, y_scale, sns.xkcd_rgb[color])

        plt.scatter(marked_x_scale, marked_y_scale, color='k', marker='x', s=50)

        texts = []
        align = 'right'
        for index,(i,j) in enumerate(zip(marked_x_scale, marked_y_scale)):
            if (align == 'right'):
                texts.append(plt.text(i+.02, j-0.01, marked_labels[index]))
                align = 'left' 
            elif (align == 'left'):
                texts.append(plt.text(i-0.02*len(marked_labels[index]), j-0.01, marked_labels[index]))
                align = 'right'



        plt.ylim(0,1)
        plt.ylabel("Actives")

        plt.xlim(0,1)
        plt.xlabel("Inactives")

        adjust_text(texts)

        plt.savefig(self.output_dir + self.name + "_roc_" + watch_role + "_" + watch_compound + ".png", dpi=300, format="png", bbox_inches="tight", transparent=True)
        plt.close()
        
        
        
    def plot_pr(self, color):
        """Plot the Precision-Recall curve.
        
        """
        
        plt.figure(figsize=(7,7))

        sns.set(context="notebook", style="ticks")
        
        f_scores = np.linspace(0.2, 0.8, num=4)
        for f_score in f_scores:
            x = np.linspace(0.01, 1)
            y = f_score * x / (2 * x - f_score)
            l, = plt.plot(x[y >= 0], y[y>=0], color='gray', alpha=0.2)
            plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

        plt.plot(self.recall, self.precision, sns.xkcd_rgb[color])

        plt.ylim(0,1)
        plt.ylabel("Precision")

        plt.xlim(0,1)
        plt.xlabel("Recall")

        plt.savefig(self.output_dir + self.name + "_pr.png", dpi=300, format="png", bbox_inches="tight", transparent=True)
        plt.close()
        


# In[41]:


x = "/home/data/mashas/tamir/T2R14_gmeiner/library3/ifd_1743-out_pv-dock_pv.maegz"


# In[42]:


s = ScreenEvaluator(x, "nomodel")


# In[43]:


n_actives = 193
n_decoys = 731


# In[44]:


s.calc_enrichment(n_actives, n_decoys)
s.precision_recall(s.roclist)
df = pd.DataFrame(s.metrics, index=[0])
df['screen'] = x


# In[45]:


s.results['score'][s.youdenj.index(max(s.youdenj))]


# In[46]:


max(s.youdenj)


# In[48]:


df.to_csv("/home/data/mashas/tamir/T2R14_gmeiner/fastrocs/top1-1743-data.csv")


# # Putting it all together

# A generation object can be a wrapper for all of the functions that need to be run on all the screens in a generation. These are:
# 1. Make receptor grids
# 2. Dock the library
# 3. Save enrichment metrics
# 4. Minimize receptors in complex with convergence ligands
# 
# One important thing to include is a way to check which steps finish successfully for each member. When there's an error, it would be convenient to have a simple way to re-run a step only on the ones that failed.
# 
# Also, a better way to schedule jobs with a "check-in, check-out" system. At the moment, there's an arbitrarily coded delay time to make sure the license server has time to update the number of jobs licenses currently using a license to ensure that we don't submit a job when there are no licenses available. However, this doesn't allow for fast and accurate submission, so the load-balancing is inefficient (especially for minimization jobs). A system which keeps track of which jobs have started, are running, and have finished would be able to more finely manage this process.

# In[59]:


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


# In[60]:


initial_members = glob.glob("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/*out_pv.mae*")


# In[61]:


len(initial_members)


# In[62]:


g1 = Generation(20, initial_members, "/home/data/mashas/tamir/T2R10_cucurbitacins/ligands/bitterdb_tp_tn_dude_decoys.mae", 45, 408)


# In[63]:


#g1.construct_screen('noligand')


# In[64]:


grids = glob.glob("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/*grid.in")


# In[65]:


g1.grids = grids


# In[56]:


#g1.run_grids()


# In[66]:


grids = glob.glob("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/*grid.zip")


# In[67]:


len(grids)


# In[68]:


# Assemble self.docks
docks = glob.glob("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/*dock.in")


# In[69]:


docks[0]


# In[70]:


len(docks)


# In[71]:


'139' in docks[0]


# In[72]:


# Remove completed runs
completed = glob.glob("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/*dock*maegz")
print(len(completed))

for x in docks:
    name = x.split('icb')[-1].split('out_pv')[0]
    for i in completed:
        if name in i:
            docks.remove(x)


# In[73]:


len(docks)


# In[74]:


name


# In[75]:


g1.docks = docks


# In[76]:


len(g1.docks)


# In[77]:


g1.run_docking()


# In[78]:


docks = glob.glob("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/*dock_pv.maegz")


# In[79]:


len(docks)


# In[80]:


# Remove from 


# In[81]:


g1.docks = docks


# In[82]:


g1.load_docking_results()


# In[83]:


g1.analyse_enrichment('g1_allmetrics.csv')


# In[ ]:


g1.metrics = pd.read_csv(g1.workdir+'g1_allmetrics.csv')


# In[ ]:


g1.workdir


# In[ ]:


g1.fetch_top_screens(10)


# In[ ]:


g1.metrics.to_csv(g1.workdir+'g1_allmetrics_with_rank.csv')


# In[ ]:


sns.set(context="notebook", style="ticks")

plt.figure(figsize=(9,6))

sns.scatterplot(x='AUC', 
                y='prauc',
                marker='.',
                s=60,
                alpha=1,
                #hue="dataset",
                #style="dataset",
                edgecolor='',
                #palette=cmap,
                color=sns.xkcd_rgb['navy'],
                data=g1.metrics)


#plt.xlim(0.4, 0.8)
plt.xlabel("ROC AUC")

#plt.ylim(0, 6.5)
plt.ylabel("Precision-Recall AUC")


plt.savefig("/home/data/mashas/tamir/T2R10_cucurbitacins/receptor_models/ifd_cluster_bitterdb/individuals/g1_enrichment.png", format='png', dpi=300, bbox_inches='tight')

plt.show()
plt.close()


# In[ ]:





# In[74]:


def plot_pr_multi(screens, color, name):
    """Plot the Precision-Recall curves for multiple provided datasets.

    """
    
    precision = []
    recall = []
    
    for s in screens:
        se = ScreenEvaluator(s, "nomodel")

        #se.write_hitlist()

        n_actives = 218
        n_decoys = 1121

        se.calc_enrichment(n_actives, n_decoys)
        se.precision_recall(se.roclist)
        
        precision.append(se.precision)
        recall.append(se.recall)

        
        
    # Plot
    plt.figure(figsize=(8,8))

    sns.set(context="poster", style="ticks")

    f_scores = np.linspace(0.2, 0.8, num=4)
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y>=0], color='gray', alpha=0.2)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.8, y[45] + 0.02))

    for (p, r, c) in zip(precision, recall, color):
        plt.plot(r, p, sns.xkcd_rgb[c])

    plt.ylim(0,1)
    plt.ylabel("Precision")

    plt.xlim(0,1)
    plt.xlabel("Recall")
    
    plt.hlines([float(10)/float(84)], 0, 1, linestyles='dashed')

    plt.savefig(name, dpi=300, format="png", bbox_inches="tight", transparent=True)
    plt.show()
    plt.close()


# In[75]:


import random


# In[76]:


n = 10

plot_pr_multi(g1.fetch_top_screens(n), random.sample(sns.xkcd_rgb.keys(), n), "/home/data/mashas/tamir/T2R14_gmeiner/antagonists/inactive_state_homology/screening/multi_align_2_ifd/library4_g1_top10_pr-multi-curve.png")


# ### Report PRAUC and MCC

# In[24]:


n = 10

for s in g1.fetch_top_screens(n):
    se = ScreenEvaluator(s, "nomodel")

    #se.write_hitlist()

    n_actives = 10
    n_decoys = 74

    se.calc_enrichment(n_actives, n_decoys)
    se.precision_recall(se.roclist)
    
    print se.metrics['AUC'], se.metrics['AUC_SD'], se.metrics['prauc'], se.metrics['mcc'], se.metrics['f1']


# In[167]:


se.metrics


# In[ ]:





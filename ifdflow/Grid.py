#!/usr/bin/env python
# coding: utf-8

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
        



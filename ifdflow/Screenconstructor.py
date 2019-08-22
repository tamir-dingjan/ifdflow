#!/usr/bin/env python
# coding: utf-8


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



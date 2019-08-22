#!/usr/bin/env python
# coding: utf-8


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
        
        

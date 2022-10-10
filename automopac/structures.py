from abc import ABC, abstractmethod
from collections.abc import Iterable

import mpmath as mpm

from symmetry import LineGroup, PointGroup, SrewAxis


class Abstract_Structure(ABC):
    '''Prototype for all structure types (Abstract factory)

    TODO:
        * Prototype for all types of structures (3D, 2D, 1D, 0D)
        * Contains information about structure
        * Interface to manipulate with structure (Supercell, monomer scrolling, N_merization)
    '''
    def __init__(self):
        pass

    @abstractmethod
    def monomer():
        pass

    @abstractmethod
    def symcell(self):
        pass

    @abstractmethod
    def xyz(self):
        pass

    @abstractmethod
    def cell(self):
        pass

    @abstractmethod
    def group(self):
        pass

class Atom:
    '''Interface for atoms for easy manipulating with structures
    '''
    def __init__(self, atom:str, coords:Iterable):
        self.atom = atom
        self.coordinates = tuple(mpm.mpf(i) for i in coords)

    @classmethod
    def from_string(line):
        """Turn XYZ-type line into atoms

        Args:
            line (srt): _description_
        """
        atom, coords = line.split()[0], line.split()[0:]
        coords = tuple(mpm.mpf(i) for i in coords)
        return Atom(atom=atom, coords=coords)

class Structure1D(Abstract_Structure):
    def __init__(self, symmetry: LineGroup, symcell: tuple) -> None:
        pass

    @property
    def xyz(self):
        pass

    def get_monomer(self):
        pass

    def get_symcell(self):
        pass

    def get_group(self):
        pass

    def get_cell():
        pass

class Molecule(Abstract_Structure):
    def __init__(self, symmetry: PointGroup, symcell: tuple) -> None:
        self.group = symmetry.group
        self.symcell = symcell 

    @property
    def xyz(self):
        pass

    def get_monomer(self):
        pass

    def get_symcell(self):
        pass

    def get_group(self):
        pass

    def get_cell():
        pass

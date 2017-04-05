from HexagonalCPSurface import HexagonalCPSurface
from HexagonalSurface import HexagonalSurface
from SquareSurface import SquareSurface
from ChainSegment import ChainSegment
from Bead import Bead
from DummyMolecule import DummyMolecule


molecule_name_map = {}
molecule_name_map['HexagonalSurface'] = HexagonalSurface
molecule_name_map['HexagonalCPSurface'] = HexagonalCPSurface
molecule_name_map['SquareSurface']    = SquareSurface
molecule_name_map['ChainSegment']     = ChainSegment

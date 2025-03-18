import json
from turtle import update

from attrs import define, field

from src.Basic.Atom import Atom
from src.Basic.utilites import make_gui_f1, ATOMS
from src.Symmetry.LineGroups import LineGroup
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.ScrewAxis import ScrewAxis


@define(kw_only=True, slots=True)
class Structure1DBase:
    # TODO: Docstring
    # TODO: Test it

    LG: LineGroup
    symcell: tuple[Atom] = field(converter=tuple)

    status: str = field(default=None)
    monomer: tuple[Atom] = field(repr=False)

    properties: dict[str, float] = field(repr=False, factory=dict)
    container: dict[str, tuple[float]] = field(repr=False, factory=dict)
    atomic_prop_container: dict[str, tuple[float]] = field(repr=False, factory=dict)


    monomer_xyz: tuple[Atom] = field(repr=False, default=None)
    structure_xyz: tuple[Atom] = field(repr=False, default=None)

    @property
    def q_p_string(self):
        return f"{self.LG.SA.q}_{self.LG.SA.p}"

    @property
    def A(self):
        return self.LG.A

    @property
    def Q(self):
        return self.LG.Q

    @property
    def f(self):
        return self.LG.f

    @monomer.default
    def _make_monomer(self):
        return self.LG.PG.apply(self.symcell)

    @property
    def structure(self):
        return self.LG.apply(self.symcell)

    @property
    def monomer_atoms(self):
        return [self.structure.index(atom) for atom in self.monomer]

    @property
    def symcell_atoms(self):
        return [self.structure.index(atom) for atom in self.symcell]


class Structure1D(Structure1DBase):
    def get_orbits(self) -> dict[int, Atom]:
        # TODO: Docstring
        # TODO: Test it

        orbits = {}
        for atom in self.structure:
            if orbits.get(atom._orbit):
                orbits[atom._orbit].append(atom)
            else:
                orbits[atom._orbit] = []
                orbits[atom._orbit].append(atom)

        return orbits

    def to_xyz(self, path):
        # TODO:
        structure = self.structure
        with open(path, "w") as wf:
            wf.write(f"{len(structure)}\n\n")
            for atom in structure:
                wf.write(str(atom) + "\n")

            x = self.A if self.LG.SA.axis == "x" else 0
            y = self.A if self.LG.SA.axis == "y" else 0
            z = self.A if self.LG.SA.axis == "z" else 0

            wf.write(f"tv\t{x:^20.12E} {y:^20.12E} {z:^20.12E}")

    def to_dict(self):
        # TODO:
        d = dict(
            LG=self.LG.to_dict(),
            symcell=[str(atom) for atom in self.symcell],
            status=self.status,
            properties=self.properties,
            container=self.container,
            monomer=[str(atom) for atom in self.monomer],
            monomer_xyz=[str(atom) for atom in self.monomer_xyz] if self.monomer_xyz else [],
            structure_xyz=[str(atom) for atom in self.structure_xyz] if self.structure_xyz else [],
            atomic_prop_container=self.atomic_prop_container,
        )
        return d

    @classmethod
    def from_dict(cls, parameters: dict):
        # TODO:
        _ = cls(
            LG=LineGroup.from_dict(parameters.get("LG")),
            symcell=[
                Atom.from_string(atom)
                for atom in parameters.get("symcell", [None])
                if atom
            ],
            status=parameters.get("status"),
            properties=parameters.get("properties"),
            container=parameters.get("container"),
            monomer_xyz=[
                Atom.from_string(atom)
                for atom in parameters.get("monomer_xyz", [None])
                if atom
            ],
            structure_xyz=parameters.get("structure_xyz"),
            monomer=[
                Atom.from_string(atom)
                for atom in parameters.get("monomer", [None])
                if atom
            ],
            atomic_prop_container=parameters.get("atomic_prop_container"),
        )
        return _

    def update(self, properties:dict = None, atomic_prop_container:dict=None, status:str=None, structure_data:dict=None):
        # TODO: Переписать этот метод, чтобы он просто обновлял данные вещества
        # Нужно задать новую лайн группу и обновить свойства
        # symcell_type = C1|PG|No
        # parameters, container, atomic_prop_container
        self.properties = properties or {}
        self.atomic_prop_container = atomic_prop_container or {}
        self.status = status or None
        structure_data = structure_data or None
        self.structure_xyz = structure_data.get('structure_xyz', [])
        A = structure_data.get('A', self.LG.A)
        new_LG = LineGroup(PG=self.LG.PG, SA=ScrewAxis(q=self.LG.SA.q, p=self.LG.SA.p, A=A))
        self.LG = new_LG

    def to_json(self, path):
        # TODO:
        with open(path, "w") as wf:
            json.dump(fp=wf, obj=self.to_dict())

    @classmethod
    def from_json(cls, path):
        # TODO:
        with open(path, "r") as fr:
            d = json.load(fp=fr)
        return cls.from_dict(d)

    def to_f34(self, path):
        import re
        srch = re.compile(r'([A-Za-z]+)')

        q, p, A = self.LG.SA.q, self.LG.SA.p, self.LG.SA.A

        coords = []

        for i in self.symcell:
            result = re.search(srch, str(i))
            atom = result.group(1)
            coords.append(str(i).replace(atom, str(ATOMS[atom])).split())

        make_gui_f1(int(q), int(p), A, coords, path)

    def symcell_from_xyz(self):
        # Нахожу длинну каждой орбиты
        # В XYZ нахожу атомы из орбит
        # Беру первый представитель каждой орбиты из XYZ
        return [self.structure_xyz[i] for i in self.symcell_atoms]

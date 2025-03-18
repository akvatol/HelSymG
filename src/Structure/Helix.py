import json
from pathlib import Path

import numpy as np
from attrs import define, field

from src.Basic.utilites import all_LG1F_generator, choose_groups
from src.QChem.Calculators.Template import CalculatorTemplate
from src.QChem.Parsers.Template import ParserTemplate
from src.Symmetry.LineGroups import LineGroup
from src.Symmetry.ScrewAxis import ScrewAxis

from .Structure1D import Structure1D


@define(slots=True)
class HelixBase:
    compounds: dict = field(repr=True) # Key Here is f'{Q}'
    init_compound: Structure1D = field(repr=True, default=None)

    @compounds.validator
    def _check_compunds(instance, attributes, value):
        if len(value) == 0:
            raise Exception('Helix must include at least one compound')

    @classmethod
    def from_compound(cls, compound:Structure1D, Q_interval:tuple[float], q_max:int, torsion_interval:tuple[float]):
        # Генерируем все LG в заданом диапазоне
        # По ним создаём объекты в compounds
        line_groups_ = all_LG1F_generator(q_max=q_max, Q_interval=Q_interval)
        line_groups = choose_groups(groups=line_groups_, interval_torsion=torsion_interval)
        # compounds
        pg = compound.LG.PG
        f = compound.f
        symcell = compound.symcell
        compounds = {i.Q: Structure1D(symcell=symcell, LG=LineGroup(PG=pg, SA=ScrewAxis(q=i.q, p=i.p, A=f*i.q))) for i in line_groups}
        compounds[compound.Q].status = 'INIT'
        return cls(compounds=compounds, init_compound=compound)

class Helix(HelixBase):

    def run(self, calculator:CalculatorTemplate, parser:ParserTemplate, root_dir:Path, template:str=None, exec:Path = None, status:str = 'Done'):
        sides = HelixIteratorInit(self.compounds, self.init_compound.Q, direction='right'), HelixIteratorInit(self.compounds, self.init_compound.Q, direction='left')

        for side in sides:
            self._run_side(calculator=calculator(), parser=parser(), side=side, exec=exec, root_dir=root_dir, template=template)

    def _run_side(self, calculator, side, parser, root_dir:Path, exec:Path=None, status:str = 'Done', template=None):
        symcell = None
        f = None
        for compound in side:
            if compound.status == status:
                continue
            else:
                if symcell:
                    # Обновляю данные чтобы расчет не заваливался из-за неправильной структуры
                    compound.symcell = symcell
                    lg_dict = compound.LG.to_dict()
                    lg_dict['SA']['A'] = compound.LG.SA.q*f
                    compound.LG = LineGroup.from_dict(lg_dict)
                # Калькулятор выдал путь к аут файлу
                path_to_out = calculator.run(root_dir=root_dir, compound=compound, exec=exec, template=template)
                # Распарсил данные
                # XXX Переделать try except
                # Это сделано на тот случай, когда мономер искозился, и больше оно работать не будет
                # ! TODO: Хранить последний удачный мономер, чтобы перезапускать расчет с него
                try:
                    data = parser.run(path_to_out)
                    data.update({'status': status})
                    compound.update(**data)
                    symcell = compound.symcell_from_xyz()
                    f = compound.LG.SA.f
                    compound.to_f34(path_to_out.parent)
                except:
                    break

    def to_dict(self):
        return {i: self.compounds[i].to_dict() for i in self.compounds}

    @classmethod
    def from_dict(cls, data:dict):
        compounds = {float(Q):Structure1D.from_dict(data[Q]) for Q in data}
        return cls(compounds)

    def to_json(self, path:Path):
        # TODO:
        with open(path, "w") as wf:
            json.dump(fp=wf, obj=self.to_dict(), cls=NpEncoder)

    @classmethod
    def from_json(cls, path:Path):
        # TODO:
        with open(path, "r") as fr:
            d = json.load(fp=fr)
        return cls.from_dict(d)

    def update(self, Q_interval:tuple[float], q_max:int):
        pass


class HelixIteratorInit:
    # Итератор для первоначального обхода по словарю. Сперва берет начальную точку. Потом ближаейшие к ней и т.д.

    def __init__(self, helix:dict, start:float, direction:str='right'):
        # TODO: Переименовать helix
        self.helix = helix
        self._position = 0
        self.helix_Q = [(key, helix[key]) for key in helix.keys()]
        # Сортирует соединения по расстоянию от начального Q

        self.helix_Q = list(filter(lambda x: start <= x[0] if direction == 'right' else start >= x[0], sorted(self.helix_Q, key= lambda xs: abs(xs[0] - start))))

    def __next__(self):
        try:
            value = self.helix_Q[self._position]
            self._position += 1
        except IndexError:
            raise StopIteration()
        return value[1]

    def __iter__(self):
        self._position = 0
        return self
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

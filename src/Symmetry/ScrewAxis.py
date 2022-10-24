from attrs import field, frozen
from src.Basic.symmetry_element import symmetry_element

from src.Basic.Templates import GroupTemplate
from src.Basic.Atom import Atom
from src.Basic.utilites import _make_srew_axis, _positive_validator, find_rp
import mpmath as mpm

mpm.mp.dps = 100


@frozen(slots=True)
class ScrewAxisBase(GroupTemplate):
    q: int = field(kw_only=True, default=1, validator=_positive_validator)
    p: int = field(kw_only=True, default=1, validator=_positive_validator)
    A: int = field(kw_only=True, validator=_positive_validator)
    axis: str = field(kw_only=True, default='x')
    r: int = field(init=False)
    Q: int = field(init=False)
    
    __generators: dict = field(init=False, repr=False)
    __group: frozenset = field(init=False, repr=False)

    @p.validator
    def _p_validation(instance, attribute, value):
        if instance.p > instance.q:
            raise ValueError('p value cannot be greater than q')

    @r.default
    def _find_r(self):
        r = find_rp(self.q, self.p)
        if r:
            return r
        else:
            raise ValueError('q and p must be coprime numbers')

    @Q.default
    def _set_Q(self) -> float:
        return self.q/self.r

    @__generators.default
    def _srew_axis_generator(self) -> symmetry_element:
        return _make_srew_axis(q=self.q, p=self.p, A=self.A, axis=self.axis)

    @__generators.validator
    def _screw_axis_validation(instance, attribute, value):
        if value.get('q').order != instance.q:
            raise ValueError(f'Generated screw axis does not constent\n q = {instance.q}, screw axis order = {value.order}')
    
    @__group.default
    def _make_screw_axis_group(self) -> frozenset:
        return self.__generators['q'].get_all_powers()

    @property
    def group(self):
        return self.__group

    @property
    def generators(self):
        return self.__generators


class ScrewAxis(ScrewAxisBase):

    def reduce(self, index: int):
        if index < self.Q:
            # TODO: Перевести
            raise ValueError("Группа не может быть меньше чем порядок оси Q")
        elif index == self.Q:
            SA = ScrewAxis(q=1, p=1, A=self.A)
        else:
            SA = self._reduce_screw_axis(index)
            SA = ScrewAxis(q=SA[0], p=SA[1], A=self.A, )
        return SA

    def _reduce_screw_axis(self, index: int ) -> list[int]:
        if self.q % index == 0:
            new_q = self.q / index
        else:
            new_q = self.q * index
            new_p = find_rp(self.q, self.r)
        return new_q, new_p

    def apply(self, atoms:list[Atom]) -> frozenset:
        strucure = frozenset([SE.apply(atom) for SE in self.group for atom in atoms])
        return strucure

    def get_stabilizer(self, atom:Atom):
        # ! Stabilizer of screw axis is L1 Group
        # ? Do i need that?
        return ScrewAxis(q=1, p=1, A=self.A)

    def get_orbit(self, atom: Atom):
        """Возвращает словарь вида {Atom:[Atom, Atom1, Atom2 ..., AtomN]}, где Atom1, Atom2 ..., AtomN получены из Atom преобразованияями симметрии

        Returns:
            dict[Atom:list[Atom]]: Словарь содержащий все орбиты 
        """
        orbit = {atom:[]}
        for elements in self.group:
            new_atom = elements.apply(atom)
            if new_atom in orbit[atom]:
                pass
            else:
                orbit[atom].append(new_atom)

        return orbit

    def to_dict(self):
        return {'q': self.q, 'p':self.p, 'A':self.A}

    @classmethod
    def from_dict(cls, parameters:dict[str, int]):
        return cls(q=parameters.get('q', 1), p=parameters.get('p', 1), A=parameters.get('A'), axis=parameters.get('axis', 'x'))

    def copy(self):
        return self.from_dict(self.to_dict())
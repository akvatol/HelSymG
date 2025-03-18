import math
import re
from collections import namedtuple
from dataclasses import dataclass
from pathlib import Path

import mpmath as mpm
import numpy as np
from mpmath import mp

from src.Basic.SymmetryElement import SymmetryElement

mpm.mp.mpds = 100

ATOMS = {'H': 1,
 'He': 2,
 'Li': 3,
 'Be': 4,
 'B': 5,
 'C': 6,
 'N': 7,
 'O': 8,
 'F': 9,
 'Ne': 10,
 'Na': 11,
 'Mg': 12,
 'Al': 13,
 'Si': 14,
 'P': 15,
 'S': 16,
 'Cl': 17,
 'Ar': 18,
 'K': 19,
 'Ca': 20,
 'Sn': 21,
 'Ti': 22,
 'V': 23,
 'Cr': 24,
 'Mn': 25,
 'Fe': 26,
 'Co': 27,
 'Ni': 28,
 'Cu': 29,
 'Zn': 30,
 'Ga': 31,
 'Ge': 32,
 'As': 33,
 'Se': 34,
 'Br': 35,
 'Kr': 36,
 'Rb': 37,
 'Sr': 38,
 'Y': 39,
 'Zr': 40}

class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args, **kwargs)
        return cls._instances[cls]

#TODO: переписать под Attrs
@dataclass(slots=True, frozen=True)
class SymmetryElements(metaclass=Singleton):
    # ? Возможно стоит перенести этот класс в отдельный файл
    I = SymmetryElement(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]))
    sigma_xz = SymmetryElement(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]]))
    sigma_yz = SymmetryElement(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    sigma_yx = SymmetryElement(rotation=mpm.matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]))
    C2x = SymmetryElement(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, -1]]))
    C2y = SymmetryElement(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]))
    C2z = SymmetryElement(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]))

    @staticmethod
    def make_cn_z(n: int) -> SymmetryElement:
        """Support function for make Cn symmetry element

        Args:
            n (int): axis order

        Returns:
            SymmetryElement: Cn group generator
        """
        angle = mpm.radians(360 / n)
        rot = mpm.matrix(
            [[mpm.cos(angle), -mpm.sin(angle), 0.0], 
            [mpm.sin(angle), mpm.cos(angle), 0.0],
            [0, 0, 1]],
        )
        trans = mpm.matrix(3, 1)
        return SymmetryElement(rotation=rot, translation=trans)

    @staticmethod
    def make_cn_x(n: int) -> SymmetryElement:
        """Support function for make Cn symmetry element

        Args:
            n (int): axis order

        Returns:
            SymmetryElement: Cn group generator
        """
        angle = mpm.radians(360 / n)
        rot = mpm.matrix(
            [
                [1, 0, 0],
                [0.0, mpm.cos(angle), -mpm.sin(angle)],
                [0.0, mpm.sin(angle), mpm.cos(angle)],
            ]
        )
        trans = mpm.matrix(3, 1)
        return SymmetryElement(rotation=rot, translation=trans)

    @staticmethod
    def make_cn_y(n: int) -> SymmetryElement:
        """Support function for make Cn symmetry element

        Args:
            n (int): axis order

        Returns:
            SymmetryElement: Cn group generator
        """
        angle = mpm.radians(360 / n)
        rot = mpm.matrix(
            [
                [mpm.cos(angle), 0, -mpm.sin(angle)],
                [0, 1, 0],
                [mpm.sin(angle), 0.0, mpm.cos(angle)],
            ]
        )
        trans = mpm.matrix(3, 1)
        return SymmetryElement(rotation=rot, translation=trans)

# TODO: refactor this
# http://www.pci.tu-bs.de/aggericke/PC4e/Kap_IV/Matrix_Symm_Op.htm
PRESETS = {
    "x": {
        "v": SymmetryElements.sigma_xz,
        "I": SymmetryElements.I,
        "h": SymmetryElements.sigma_yz,
        "U": SymmetryElements.C2y,
    },
    "y": {
        "v": SymmetryElements.sigma_yz,
        "I": SymmetryElements.I,
        "h": SymmetryElements.sigma_xz,
        "U": SymmetryElements.C2z,
    },
    "z": {
        "v": SymmetryElements.sigma_xz,
        "I": SymmetryElements.I,
        "h": SymmetryElements.sigma_yx,
        "U": SymmetryElements.C2x
    },
}

def all_LG1F_generator(q_max: int, Q_interval: tuple, direction='left', n=1) -> tuple:
    """Generate all LG with q <= q_max, n = n, and Q_interval[0] <= Q <= Q_interval[1].

    Args:
        q_max (int): Maximal values of q_tilde (q_tilde = q/n)
        n (int): Index n of Cn group of monomer
        Q_interval (tuple): Q_interval - tuple with minimal and maximal Q

    Returns:
        tuple: list of sorted line groups (srew-axis if n = 1)
    """

    # (r*p)%q = +-1 original formula from Damnjanovich
    q_max = q_max + 1

    LineGroup1 = namedtuple("LineGroup1", "Q q p r torsion")

    q = np.arange(0, q_max)
    r = np.arange(0, q_max - 1)
    data = []

    rp_mod_q = np.mod.outer(np.multiply.outer(r, r), q)
    indexes = np.nonzero(np.where(rp_mod_q == 1, rp_mod_q, 0))


    # q here - q_tilde, p - p_tilde
    for r, p, q in zip(indexes[0], indexes[1], indexes[2]):
        if q < r or q < p:
            continue
        q = q*n/math.gcd(q, n)
        Q = q / r
        if Q >= Q_interval[0] and Q <= Q_interval[1]:
            A = 360/Q
            data.append(LineGroup1(Q, q, p*n % q, r, round(A, 5)))

    if direction == 'left':
      pass
    else:
        new_data = []
        for i in data:
            new_data.append(LineGroup1(i.Q, i.q, p=i.q-i.p, r=i.r, torsion=round(360/i.Q, 5)))
        data = new_data

    ans = tuple(sorted(data, key=lambda x: x.Q, reverse=True))

    return ans

def choose_next_group(group, groups:tuple, interval_torsion:tuple):
    t = group.torsion

    possible_groups = []

    for i, j in enumerate(groups):
        if interval_torsion[0] < abs(j.torsion - t) < interval_torsion[1]:
            possible_groups.append((i, j))

    if possible_groups:
        return min(possible_groups, key=lambda x: x[1].q)
    else:
        return (None, None) 

def choose_groups(groups:tuple, interval_torsion:tuple):
    data = []
    i_max = len(groups)
    data.append(groups[0])

    i_current = 1
    while i_current:
        current_groups = data[-1]
        i_current_, target_group = choose_next_group(current_groups, groups[i_current:], interval_torsion)
        if i_current_:
            i_current += i_current_
            data.append(target_group)
        else:
            break

    if groups[-1] in data:
        pass
    else:
        data.append(groups[-1])

    return data


def find_rp(q: int, p: int) -> int:
    """Return r or p depends on given arguments: if f(q, p) --> r else f(q, r) --> p;

    Args:
        q (int): q - is a q from line group theory (order of a screw-axis)
        p (int): in this case it may be r or p from LGT

    Returns:
        int: r or p depends on input
    """

    if (q == 1) and (p == 1):
        return 1

    for r in range(1, 100):
        mod_q = r * p % q
        if mod_q == 1:
            return int(r)


def point_group_symbol_parser(symbol: str) -> tuple:
    """Takes Schönflies symbol and return it`s readable representation

    Args:
        symbol (str): Schönflies symbol of point group

    Returns:
        tuple: Schönflies symbol representation
    """

    regex = r"([CSD]?)(\d+)([VDH]?)|(C[IS])"
    symbol = symbol.upper()
    search = re.search(regex, symbol)
    if search:
        if search.group(4):
            return search.group(4), "", ""
        else:
            s1 = search.group(1)
            s2 = search.group(2)
            s3 = search.group(3)
            return s1, s2, s3

def make_generators(parameters: dict) -> dict[str, SymmetryElement]:
    # TODO: DOCSTRING
    """_summary_

    Returns:
        _type_: _description_
    """

    axis = parameters.get('axis', 'x')
    preset = PRESETS.get(axis)

    if axis == 'x':
        make_cn = SymmetryElements.make_cn_x(n=parameters.get('n', 1))
    elif axis == 'z':
        make_cn = SymmetryElements.make_cn_z(n=parameters.get('n', 1))
    else:
        make_cn = SymmetryElements.make_cn_y(n=parameters.get('n', 1))

    preset.update(dict(n=make_cn))

    generators_ = {name: preset.get(name) for name in parameters if parameters.get(name) and name != 'axis'}

    return generators_

def make_group(generators: dict) -> tuple:
    # TODO Docstring
    # TODO: regularize orders of group elements
    """_summary_

    Args:
        generators (frozenset): _description_

    Returns:
        frozenset: _description_
    """
    generators_elements = [generators[element].get_all_powers() for element in generators]
    generators_elements.sort(key=len, reverse=True)
    group = list(generators_elements[0])

    # G = x1*H + x2*H + ... + xn*H

    # optimize it
    for subgroup in generators_elements[1:]:
        for s_element in subgroup:
            if s_element in group:
                continue
            else:
                for element in group:
                    new_element = element*s_element
                    if new_element in group:
                        continue
                    else:
                        group.append(new_element)

    return tuple(group)

def detect_group(group: frozenset, axis: str) -> dict:
    # TODO: make it
    preset = PRESETS.get(axis)
    group_parameters = {}
    for element in preset:
        if element == 'n':
            pass

def _positive_validator(instance, attribute, value):
    # TODO; Docstring
    # Validator for point and line groups parameters
    if value <= 0:
        raise ValueError('Value vector must be positive')

def _make_srew_axis(q:int, p:int, A:float, axis:str = 'x') -> dict[str, SymmetryElement]:
    # TODO: Docstring
    rotation_ = make_generators({'n':q, 'axis':axis})
    rotation = rotation_['n'].rotation
    translation_part = mpm.fmul(A, (mpm.fdiv(p, q)))
    screw_generators = {
        'x':SymmetryElement(rotation=rotation, translation=mpm.matrix([translation_part, 0, 0]), translation_vector=mpm.matrix([A, 0, 0]),),
        'y':SymmetryElement(rotation=rotation, translation=mpm.matrix([0, translation_part, 0]), translation_vector=mpm.matrix([0, A, 0]),),
        'z':SymmetryElement(rotation=rotation, translation=mpm.matrix([0, 0, translation_part]), translation_vector=mpm.matrix([0, 0, A]),),
        }

    return {'q':screw_generators[axis]}

def list_to_str(list_):
    """Переводит mpf в строку"""
    # a = map(list_, lambda x: str(x))
    return list(map(lambda x: f'{float(x):16.12f}', list_))

def make_gui_f1(q: int, p: int, A: float, coords: list, folder_path: str, name: str = None):
    """
    Создает файл .gui в указанной папке.

    :param q: Параметр q.
    :param p: Параметр p.
    :param A: Параметр A.
    :param coords: Список координат.
    :param folder_path: Путь к папке, где будет создан файл.
    :param name: Имя файла (опционально). Если не указано, используется стандартное имя.
    """
    # Преобразуем строку пути в объект Path
    folder = Path(folder_path)

    # Убедимся, что папка существует
    folder.mkdir(parents=True, exist_ok=True)

    # Если имя файла не указано, используем стандартное имя
    if name is None:
        name = f'{q}_{p}'
    else:
        name = f'{q}_{p}_{name}'

    # Полный путь к файлу
    file_path = folder / f'{name}.gui'  # Используем оператор / для объединения путей

    # Вычисления
    angle = mp.mpf(2 * mp.pi / q)
    A = mp.mpf(A)

    template = (
        f" 1  1  1   E=0.000\n"
        f"{float(A):16.12f} {0.0:16.12f} {0.0:16.12f}\n"
        f"{0.0:16.12f} {500.0:16.12f} {0.0:16.12f}\n"
        f"{0.0:16.12f} {0.0:16.12f} {500.0:16.12f}\n"
        f" {int(q)}\n"
    )

    with open(file_path, 'w') as fw:
        fw.write(template)
        for i in range(q):
            fw.write(' '.join([' '] + list_to_str([1, 0, 0])) + '\n')
            fw.write(' '.join([' '] + list_to_str([
                0.,
                mp.cos(angle * i),
                -mp.sin(angle * i)
            ])) + '\n')
            fw.write(' '.join([' '] + list_to_str([
                0.,
                mp.sin(angle * i),
                mp.cos(angle * i)
            ])) + '\n')
            fw.write(' '.join([' '] + list_to_str([
                A * mp.mpf(i * p / q % 1),
                0.,
                0.
            ])) + '\n')
        fw.write(f' {len(coords)}\n')
        for i in coords:
            fw.write(f" {float(i[0]):.0f} {float(i[1]):14.8f} {float(i[2]):14.8f} {float(i[3]):14.8f}\n")
        fw.write('   0    0\n')  # для правильного формата файла
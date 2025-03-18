"""Module contain Calculator for MOPAC programm.

See:
* https://github.com/openmopac/mopac
* http://openmopac.net/manual/index.html
"""

import os
import subprocess
from pathlib import Path

import mpmath as mpm

from src.Structure.Structure1D import Structure1D

from .Template import CalculatorTemplate

BASIC_TEMPLATE_MOPAC = 'AUX LARGE CHARGE=0 SINGLET GEO-OK NOREOR PDBOUT PRTXYZ ITRY=500 PM6-D3H4X NOSYM RECALC=5 MOZYME PRECISE=500'
class MopacCalculator():
    """Простой калькулятор для запуска расчетов в прогамме MOPAC. Запуск расчетов осуществляется методом run.
    """

    def _make_template(self, compound:Structure1D, template: str = BASIC_TEMPLATE_MOPAC) -> str:
        """Создаёт текст для входного файла MOPAC

        Args:
            compound (Structure1D): Вещество.
            template (str, optional): Первая строка входного файла (Параметры расчета). Defaults to BASIC_TEMPLATE_MOPAC.

        Returns:
            str: _description_
        """
        
        if template:
            pass
        else:
            template = BASIC_TEMPLATE_MOPAC

        template =[template, '\n', '\n\n']
        for atom in compound.structure:
            atom_ = str(atom).split()
            _ = f'{atom_[0]} {atom_[1]} 1 {atom_[2]} 1 {atom_[3]} 1\n'
            template.append(_)
        A = compound.A
        tv = f'tv {str(A) + " 1" if compound.LG.SA.axis == "x" else "0.0 0"} {str(A) + " 1" if compound.LG.SA.axis == "y" else "0.0 0"} {str(A) + " 1" if compound.LG.SA.axis == "y" else "0.0 0"}\n'
        template.append(tv)
        return ''.join(template)

    def _save_input(self, template:str, compound:Structure1D, root_dir:Path) -> Path:
        """Создаёт входной файл для Mopac.

        Args:
            template (str): Содержание входного файла.
            compound (Structure1D): Вещество.
            root_dir (Path): Путь к папке где запускатся расчеты.

        Returns:
            new_path(Path): путь выходному файлу.
        """        
        angle = mpm.nstr(mpm.fdiv(360, compound.LG.SA.Q), n=8)
        angle_q_p = f'{angle}_{compound.LG.SA.q}_{compound.LG.SA.p}'
        # This calc dir
        new_path = root_dir / angle_q_p

        # Создали папку
        if new_path.is_dir():
            pass
        else:
            os.mkdir(new_path)

        # Cохранили папку
        with open(new_path / 'input.mop', 'w') as fw:
            fw.write(template)

        return new_path # calc_dir

    def _run_task(self, calc_dir:Path, exec:str):
        """Запускает расчет.

        Args:
            calc_dir (Path): Путь к папке со входным файлом. (Имя входного файла предпологается как input.mop)
            exec (str): Путь к выполняемому файлу Mopac.

        Returns:
            (str): путь к выходному файлу.
        """        
        subprocess.call([exec, calc_dir/'input.mop'])
        return calc_dir/'input.out'

    def run(self, root_dir:Path, compound:Structure1D, exec:str = None, template:str = BASIC_TEMPLATE_MOPAC) -> str:
        """Пайплайн, который запускает квантовохимический расчет.

        Args:
            root_dir (Path): _description_
            compound (Structure1D): _description_
            exec (str, optional): _description_. Defaults to None.

        Returns:
            (str): путь к выходному файлу.
        """
        # Make intput
        calc_dir = self._save_input(self._make_template(compound=compound, template=template), compound=compound, root_dir=root_dir)
        # Run it
        exec = exec or 'mopac'
        path_to_out = self._run_task(calc_dir=calc_dir, exec=exec)
        return path_to_out

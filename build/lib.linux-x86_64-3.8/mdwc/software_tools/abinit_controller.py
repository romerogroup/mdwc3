#!/usr/bin/env python

import subprocess
import time
import re
import os.path
import numpy as np

masses = [
    0.0,
    1.00794,
    4.002602,
    6.941,
    9.012182,
    10.811,
    12.011,
    14.00674,
    15.9994,
    18.9984032,
    20.1797,
    22.989768,
    24.3050,
    26.981539,
    28.0855,
    30.973762,
    32.066,
    35.4527,
    39.948,
    39.0983,
    40.078,
    44.955910,
    47.88,
    50.9415,
    51.9961,
    54.93805,
    55.847,
    58.93320,
    58.69,
    63.546,
    65.39,
    69.723,
    72.61,
    74.92159,
    78.96,
    79.904,
    83.80,
    85.4678,
    87.62,
    88.90585,
    91.224,
    92.90638,
    95.94,
    98.9062,
    101.07,
    102.9055,
    106.42,
    107.8682,
    112.411,
    114.82,
    118.710,
    121.753,
    127.60,
    126.90447,
    131.29,
    132.90543,
    137.327,
    138.9055,
    140.115,
    140.90765,
    144.24,
    147.91,
    150.36,
    151.965,
    157.25,
    158.92534,
    162.50,
    164.93032,
    167.26,
    168.93421,
    173.04,
    174.967,
    178.49,
    180.9479,
    183.85,
    186.207,
    190.2,
    192.22,
    195.08,
    196.96654,
    200.59,
    204.3833,
    207.2,
    208.98037,
    209.0,
    210.0,
    222.0,
    223.0,
    226.0254,
    230.0,
    232.0381,
    231.0359,
    238.0289,
    237.0482,
    242.0,
    243.0,
    247.0,
    247.0,
    249.0,
    254.0,
    253.0,
    256.0,
    254.0,
    257.0,
    260.0,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
]


def get_nat_mass_latvec_in_strten_in(path_to_file):
    data = open(path_to_file).read()
    nat = int(re.findall('natom\s+([0-9]+)', data)[0])
    typat_str_0 = '\s+typat\s+'
    typat_srt_1 = '.+'
    typat_str = typat_str_0 + '(' + typat_srt_1 + ')'
    typat = list(map(int, re.findall(typat_str, data)[0].split()))
    # znucl= map(float, re.findall('\s+znucl((?:\s+\d+.\d+\s+)+)',data))
    znucl = list(
        map(int, list(map(float, re.findall('\s+znucl\s+(.+)', data)[0].split())))
    )
    while len(typat) < nat:
        typat_srt_1 = typat_srt_1 + '\n.+'
        typat_str = typat_str_0 + '(' + typat_srt_1 + ')'
        # x= re.findall('\s+typat\s+(.+)',data)
        # print x
        typat = list(map(int, re.findall(typat_str, data)[0].split()))
    mass = []
    for i in typat:
        mass.append(masses[znucl[i - 1]])
    mass = np.array(mass)
    a1 = list(
        map(
            float,
            re.findall('R.1.=\s*(.\d+...\d+\s+.\d+...\d+\s+.\d+...\d+)', data)[
                0
            ].split(),
        )
    )
    a2 = list(
        map(
            float,
            re.findall('R.2.=\s*(.\d+...\d+\s+.\d+...\d+\s+.\d+...\d+)', data)[
                0
            ].split(),
        )
    )
    a3 = list(
        map(
            float,
            re.findall('R.3.=\s*(.\d+...\d+\s+.\d+...\d+\s+.\d+...\d+)', data)[
                0
            ].split(),
        )
    )
    latvec_in = np.array([a1, a2, a3]).T
    latvec_in.astype('float64')
    strten_in = []
    strten_in.append(np.float64(re.findall('sigma.1\s+1.=(\s+.\d+.\d+..\d+)', data)[0]))
    strten_in.append(np.float64(re.findall('sigma.2\s+2.=(\s+.\d+.\d+..\d+)', data)[0]))
    strten_in.append(np.float64(re.findall('sigma.3\s+3.=(\s+.\d+.\d+..\d+)', data)[0]))
    strten_in.append(np.float64(re.findall('sigma.3\s+2.=(\s+.\d+.\d+..\d+)', data)[0]))
    strten_in.append(np.float64(re.findall('sigma.3\s+1.=(\s+.\d+.\d+..\d+)', data)[0]))
    strten_in.append(np.float64(re.findall('sigma.2\s+1.=(\s+.\d+.\d+..\d+)', data)[0]))
    strten_in = np.array(strten_in)
    return nat, mass, latvec_in, strten_in


def get_xred_fcart(path_to_file, nat):
    # 1 Ha/Bohr3 = 29421.02648438959 GPa
    data = open(path_to_file).readlines()
    for n, line in enumerate(data):
        if re.findall('reduced\s+coordinates\s+.array\s+xred', str(line)):
            xred_temp = data[n + 1 : n + 1 + nat]
            xred = np.array(
                [list(map(float, i.split('\n')[0].split())) for i in xred_temp]
            ).T
            xred.astype('float64')
        elif re.findall('cartesian\s+forces\s+.hartree.bohr', str(line)):
            fcart_temp = data[n + 1 : n + 1 + nat]
            fcart = np.array(
                [list(map(float, i.split('\n')[0].split())) for i in fcart_temp]
            )[:, 1:]
            fcart = fcart.T
            fcart.astype('float64')
        elif re.findall('>>>>>>>>>\s+Etotal=\s+.\d+', str(line)):  # hartree
            ener = re.findall('>>>>>>>>>\s+Etotal=(\s+.\d+.\d+..\d+)', str(line))
            ener = np.float64(ener[0])
        elif re.findall('Pressure=\s+\d+.\d+..\d+', str(line)):  # this preassure in GPa
            pressure = re.findall('Pressure=(\s+\d+.\d+..\d+)', str(line))
            pressure = np.float64(pressure[0])
    return xred, fcart


def get_md_parameters(path_to_file):
    data = open(path_to_file).read()
    Qmass = float(re.findall('\s*Qmass\s+(\d+.\d*)', data)[0])
    # temp= float(re.findall('\s*temp\s+(\d+.\d*)', data)[0])
    bmass = float(re.findall('\s*bmass\s+(\d+.\d*)', data)[0])
    Pressure = float(re.findall('\s*Pressure\s+(\d+.\d*)', data)[0])
    dt = float(re.findall('\s*dt\s+(\d+.\d*)', data)[0])
    correct_spteps = int(re.findall('\s*correct_spteps\s+(\d+)', data)[0])
    md_steps = int(re.findall('\s*md_steps\s+(\d+)', data)[0])
    abinit_steps = int(re.findall('\s*abinit_steps\s+(\d+)', data)[0])
    # return Qmass, temp, bmass, Pressure, dt, correct_spteps, md_steps, abinit_steps
    return Qmass, bmass, Pressure, dt, correct_spteps, md_steps, abinit_steps


def temp_data_reader(path_to_file, md_total):
    data_file = open(path_to_file, 'r')
    data_file = data_file.readlines()

    patt_cons = re.compile('^temp_cons')
    n_cons = 0
    patt_line = re.compile('^temp_line')
    n_line = 0
    patt_plat = re.compile('^temp_plat')
    n_plat = 0
    for string in data_file:
        if patt_cons.search(string):
            n_cons = 1
        if patt_line.search(string):
            n_line = 1
        if patt_plat.search(string):
            n_plat = 1

    data_file = str(data_file)
    indicator = n_cons + n_line + n_plat
    if indicator == 0:
        print('select a type of temperature control')
        temp_arra = []
        return temp_arra
    elif indicator == 1:
        if n_cons == 1:
            temp = float(re.findall('temp_cons([ 0-9 \. \s ,]*)', data_file)[0])
            arra = np.ones(md_total, dtype=float)
            temp_arra = temp * arra
            return temp_arra
        elif n_line == 1:
            temp = re.findall('temp_line([ 0-9 \. \s ,]*)', data_file)[0]
            temp = temp.split(',')
            temp = list(map(float, temp[: len(temp) - 1]))
            delta = (temp[1] - temp[0]) / md_total
            temp_arra = np.arange(temp[0], temp[1], delta)
            return temp_arra
        elif n_plat == 1:
            temp = re.findall('temp_plat([ 0-9 \. \s ,]*)', data_file)[0]
            temp = temp.split(',')
            temp = list(map(float, temp[: len(temp) - 1]))

            step = re.findall('temp_step([ 0-9 \. \s ,]*)', data_file)[0]
            step = step.split(',')
            step = list(map(int, step[: len(step) - 1]))

            for i, t in enumerate(temp):
                if i == 0:
                    temp_arra = t * np.ones(step[i], dtype=float)
                else:
                    temp_arra = np.concatenate(
                        (temp_arra, t * np.ones(step[i], dtype=float)), axis=0
                    )
            return temp_arra
    elif indicator == 2:
        print('only one type of temperature control')
        temp_arra = []
        return temp_arra
    elif indicator == 3:
        print('only one type of temperature control')
        temp_arra = []
        return temp_arra


def get_md_constrains(path_to_file):
    data = open(path_to_file).readlines()
    # number_constrains= []
    bond_const = []
    angl_const = []
    cell_para_const = []
    cell_angl_const = []
    atom_fix_const = []
    atom_fix_valu = []
    atom_fix_cord = []
    for line in data:
        if re.match('number_bond_constrains', line):
            numb_bond_cons = int(
                re.findall('\s*number_bond_constrains\s+(\d+)', line)[0]
            )
            # number_constrains.append(numb_bond_cons)
        elif re.match('bond_constrains', line):
            if numb_bond_cons != 0:
                data = line.split('#')
                data = re.findall('\s*bond_constrains\s+(.+)', data[0])[0].split(',')
                for const in data[:numb_bond_cons]:
                    bond_const.append(list(map(int, const.split())))
                bond_const = np.array(bond_const)
            else:
                bond_const = np.array([[1, 2]])
        elif re.match('bond_distance', line):
            if numb_bond_cons != 0:
                bond_valu = re.findall('\s*bond_distance\s+(.+)', line)[0].split(',')
                bond_valu = list(map(float, bond_valu))
                bond_valu = np.array(bond_valu)
                bool_bond_cons = 1
            else:
                bond_valu = np.array([1.00])
                numb_bond_cons = 1
                bool_bond_cons = 0
        # ***************************************
        if re.match('number_atom_fix_constrains', line):
            numb_atom_fix_cons = int(
                re.findall('\s*number_atom_fix_constrains\s+(\d+)', line)[0]
            )
            # number_constrains.append(numb_bond_cons)
        elif re.match('atom_fixed_constrains', line):
            if numb_atom_fix_cons != 0:
                data = line.split('#')
                data = re.findall('\s*atom_fixed_constrains\s+(.+)', data[0])[0].split(
                    ','
                )
                for const in data[:numb_atom_fix_cons]:
                    atom_fix_const.append(list(map(int, const.split())))
                atom_fix_const = np.array(atom_fix_const)
            else:
                atom_fix_const = np.array([[1]])

        elif re.match('atom_fixed_position', line):
            if numb_atom_fix_cons != 0:
                data = str(re.findall('\s*atom_fixed_position\s+(.+)', line)[0]).split(
                    ','
                )
                for const in data[:numb_atom_fix_cons]:
                    atom_fix_valu.append(list(map(float, const.split())))
                atom_fix_valu = np.array(atom_fix_valu)
                bool_atom_fix_cons = 1
            else:
                atom_fix_valu = np.array([[1.00, 1.00, 1.00]])
                numb_atom_fix_cons = 1
                bool_atom_fix_cons = 0

        elif re.match('atom_fix_coordinate', line):
            if numb_atom_fix_cons != 0:
                data = str(re.findall('\s*atom_fix_coordinate\s+(.+)', line)[0]).split(
                    ','
                )
                for const in data[:numb_atom_fix_cons]:
                    atom_fix_cord.append(
                        list(map(float, list(map(int, const.split()))))
                    )
                atom_fix_cord = np.array(atom_fix_cord)
            else:
                atom_fix_cord = np.array([[1.00, 1.00, 1.00]])
        # ***************************************
        elif re.match('number_angle_constrains', line):
            numb_angl_cons = int(
                re.findall('\s*number_angle_constrains\s+(\d+)', line)[0]
            )
            # number_constrains.append(numb_angl_cons)
        elif re.match('angle_constrains', line):
            if numb_angl_cons != 0:
                data = line.split('#')
                data = re.findall('\s*angle_constrains\s+(.+)', data[0])[0].split(',')
                for const in data[:numb_angl_cons]:
                    angl_const.append(list(map(int, const.split())))
                angl_const = np.array(angl_const)
            else:
                angl_const = np.array([[1, 2, 3]])
        elif re.match('value_cosine_angle', line):
            if numb_angl_cons != 0:
                angl_valu = re.findall('\s*value_cosine_angle\s+(.+)', line)[0].split(
                    ','
                )
                angl_valu = list(map(float, angl_valu))
                angl_valu = np.array(angl_valu)
                bool_angl_cons = 1
            else:
                angl_valu = np.array([1.00])
                numb_angl_cons = 1
                bool_angl_cons = 0

        elif re.match('number_cell_parameter_constrain', line):
            numb_cell_para_cons = int(
                re.findall('\s*number_cell_parameter_constrain\s+(\d+)', line)[0]
            )
            # number_constrains.append(numb_cell_para_cons)
        elif re.match('cell_parameter_constrain', line):
            if numb_cell_para_cons != 0:
                data = line.split('#')
                data = re.findall('\s*cell_parameter_constrain\s+(.+)', data[0])[
                    0
                ].split(',')
                for const in data[:numb_cell_para_cons]:
                    cell_para_const.append(list(map(int, const.split())))
                cell_para_const = np.array(cell_para_const)
            else:
                cell_para_const = np.array([[1]])
        elif re.match('cell_parameter_value', line):
            if numb_cell_para_cons != 0:
                cell_para_valu = re.findall('\s*cell_parameter_value\s+(.+)', line)[
                    0
                ].split(',')
                cell_para_valu = list(map(float, cell_para_valu))
                cell_para_valu = np.array(cell_para_valu)
                bool_cell_para_cons = 1
            else:
                cell_para_valu = np.array([1.00])
                numb_cell_para_cons = 1
                bool_cell_para_cons = 0

        elif re.match('number_cell_angle_constrain', line):
            numb_cell_angl_cons = int(
                re.findall('\s*number_cell_angle_constrain\s+(\d+)', line)[0]
            )
            # number_constrains.append(numb_cell_angl_cons)
        elif re.match('cell_angle_constrain', line):
            if numb_cell_angl_cons != 0:
                data = line.split('#')
                data = re.findall('\s*cell_angle_constrain\s+(.+)', data[0])[0].split(
                    ','
                )
                for const in data[:numb_cell_angl_cons]:
                    cell_angl_const.append(list(map(int, const.split())))
                cell_angl_const = np.array(cell_angl_const)
            else:
                cell_angl_const = np.array([[1, 2]])
        elif re.match('value_cosine_cell_angle', line):
            if numb_cell_angl_cons != 0:
                cell_angl_valu = re.findall('\s*value_cosine_cell_angle\s+(.+)', line)[
                    0
                ].split(',')
                cell_angl_valu = list(map(float, cell_angl_valu))
                cell_angl_valu = np.array(cell_angl_valu)
                bool_cell_angl_cons = 1
            else:
                cell_angl_valu = np.array([1.00])
                numb_cell_angl_cons = 1
                bool_cell_angl_cons = 0
        elif re.match('volume_constrain', line):
            volu_cons = int(re.findall('\s*volume_constrain\s+(\d+)', line)[0])
            # number_constrains.append(volu_cons)
        elif re.match('volume_value', line):
            if volu_cons != 0:
                volu_valu = re.findall('\s*volume_value\s+(.+)', line)[0].split(',')
                volu_valu = list(map(float, volu_valu))
                volu_valu = np.array(volu_valu)
                bool_volu = 1
            else:
                volu_valu = np.array([1.00])
                volu_cons = 1
                bool_volu = 0
    return (
        numb_bond_cons,
        bool_bond_cons,
        numb_angl_cons,
        bool_angl_cons,
        numb_cell_para_cons,
        bool_cell_para_cons,
        numb_cell_angl_cons,
        bool_cell_angl_cons,
        volu_cons,
        bool_volu,
        numb_atom_fix_cons,
        bool_atom_fix_cons,
        bond_const,
        angl_const,
        cell_para_const,
        cell_angl_const,
        atom_fix_const,
        bond_valu,
        angl_valu,
        cell_para_valu,
        cell_angl_valu,
        volu_valu,
        atom_fix_valu,
        atom_fix_cord,
    )


def check_out_exit_and_complet(path_to_output):
    if os.path.exists(path_to_output):
        # time.sleep(30)
        output = open(path_to_output, 'r')
        data = output.read()
        match = re.search('Calculation\s+completed', data)
        if match:
            output.close()
            n = 1
            return n
        output.close()
        time.sleep(30)


def create_directories(name, xred, h, numb, pwd=None):
    if pwd == None:
        pwd = subprocess.check_output(['pwd'])
    new_dir = pwd + '/' + name + str(numb)
    subprocess.check_call(['mkdir', new_dir])
    path_to_prot_in = pwd + '/' + name + '.in'
    path_to_prot_files = pwd + '/' + name + '.files'
    path_to_new_in = new_dir + '/' + name + str(numb) + '.in'
    path_to_new_files = new_dir + '/' + name + str(numb) + '.files'
    from_prototype_in_to_in(path_to_prot_in, path_to_new_in, xred, h)
    from_prototype_file_to_file(path_to_prot_files, path_to_new_files, numb)
    return new_dir


def from_prototype_in_to_in_step0(path_to_prot_in, path_to_new_in):
    # get in from prototipe in
    data_prot = open(path_to_prot_in).readlines()
    data_new = open(path_to_new_in, 'w')

    for n, line in enumerate(data_prot):
        data_new.write(line)
    data_new.close()
    return None


def from_prototype_in_to_in(path_to_prot_in, path_to_new_in, xred, h):
    # get in from prototipe in
    pattern1 = re.compile("^[a-z]")
    pattern2 = re.compile("^#")
    data_prot = open(path_to_prot_in).readlines()
    data_new = open(path_to_new_in, 'w')
    a1 = np.linalg.norm(h[:, 0])
    a2 = np.linalg.norm(h[:, 1])
    a3 = np.linalg.norm(h[:, 2])
    a1v = np.divide(h[:, 0], a1)
    a2v = np.divide(h[:, 1], a2)
    a3v = np.divide(h[:, 2], a3)

    for n, line in enumerate(data_prot):
        if pattern1.match(line) or pattern2.match(line):
            if re.findall('\s*acell', str(line)):
                string = 'acell %0.5f %0.5f %0.5f\n' % (a1, a2, a3)
                data_new.write(string)

            elif re.findall('\s*rprim', str(line)):
                string = 'rprim     %0.5f %0.5f %0.5f\n' % (a1v[0], a2v[0], a3v[0])
                for i in range(2):
                    string = string + '        %0.5f %0.5f %0.5f\n' % (
                        a1v[i + 1],
                        a2v[i + 1],
                        a3v[i + 1],
                    )
                data_new.write(string)
            elif re.findall('\s*xred', str(line)):
                string = 'xred   \n'
                for red_coor in xred.T:
                    string = string + '%0.5f %0.5f %0.5f\n' % (
                        red_coor[0],
                        red_coor[1],
                        red_coor[2],
                    )
                data_new.write(string)
            else:
                data_new.write(line)
    data_new.close()
    return None


def from_prototype_file_to_file(path_to_prot_files, path_to_new_files, numb):
    data_prot = open(path_to_prot_files).readlines()
    data_new = open(path_to_new_files, 'w')
    for n, line in enumerate(data_prot):
        if re.findall('/', line):
            data_new.write(line)
        else:
            string = line.split('*')
            if len(string) == 2:
                data_new.write(string[0] + str(numb) + string[1])
    data_new.close()
    return None

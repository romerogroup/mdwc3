Qmass          0.5
temp_cons      250.89
#temp_line      250.6, 275.6,
#temp_plat      250.0, 280.0, 290.5,
#temp_step      5,   2,   3,
bmass          10.0
Pressure       0.00025
dt             0.01
correct_spteps 8
md_steps       2
abinit_steps   5

number_atom_fix_constrains 0
atom_fixed_constrains 6, 7, #costrains over atoms 1 and 2
atom_fixed_position 14.96 4.57 0.00, 14.96 0.00 4.77,
atom_fix_coordinate 1 1 1, 1 1 1,

number_bond_constrains 0
bond_constrains 1 2, 1 3, #costrains between atoms 1 and 2, 7 and 8
bond_distance 6.76275, 6.89948

number_angle_constrains 0
angle_constrains 1 2 3, 4 5 6,#constrain of angle formed 3 4 5 with vertex at 3
value_cosine_angle 0.532, 0.915

number_cell_parameter_constrain 1
cell_parameter_constrain 1, 2, #constrain over the length of the first cell vector
cell_parameter_value 5.4447483764417779, 5.4607816818371040

number_cell_angle_constrain 0
cell_angle_constrain 1 2, 2 3, #constrain over the angle between the 2 and 3 cell vectors
value_cosine_cell_angle 0.0000, 0.0000

volume_constrain 0
volume_value 1738.634

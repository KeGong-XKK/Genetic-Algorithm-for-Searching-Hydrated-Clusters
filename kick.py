import numpy as np
from MOPAC_opt import write_mopac_input,run_mopac
from extract_energy_coord_MOPAC import extract_data_from_mopac_output,extract_data_from_mopac_output_after_mutation
import os
def create_rotation_matrix():
    """生成一个随机旋转矩阵"""
    theta = np.random.uniform(0, 2*np.pi)
    phi = np.random.uniform(0, np.pi)
    z = np.random.uniform(-1, 1)

    r = np.sqrt(1 - z**2)
    x, y, z = r * np.cos(phi), r * np.sin(phi), z

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    one_minus_cos_theta = 1 - cos_theta

    rotation_matrix = np.array([
        [cos_theta + x*x*one_minus_cos_theta, x*y*one_minus_cos_theta - z*sin_theta, x*z*one_minus_cos_theta + y*sin_theta],
        [y*x*one_minus_cos_theta + z*sin_theta, cos_theta + y*y*one_minus_cos_theta, y*z*one_minus_cos_theta - x*sin_theta],
        [z*x*one_minus_cos_theta - y*sin_theta, z*y*one_minus_cos_theta + x*sin_theta, cos_theta + z*z*one_minus_cos_theta]
    ])
    return rotation_matrix

def apply_transformation(coordinates, rotation_matrix, displacement):
    """应用旋转和位移变换"""
    transformed = np.dot(coordinates, rotation_matrix) + displacement
    return transformed

def write_xyz_file(filename, atoms):
    """将原子坐标写入XYZ文件"""
    with open(filename, 'a') as file:
        for atom in atoms:
            file.write(f"{atom[0]} {atom[1][0]:.3f} {atom[1][1]:.3f} {atom[1][2]:.3f}\n")

def mutation(file_path,index_kick,number_atom):
    file_MOP = f"cluster_{index_kick}.mop"
    file_xyz = f"cluster_{index_kick}.xyz"
    with open(file_MOP,'a') as file:
        file.write("PM6-D3H4 precise CHARGE=-1 SINGLET")
        file.write('\n')
        file.write('\n')
        file.write('\n')
    file.close()

    with open(file_xyz,'a') as file:
        file.write(str(number_atom) + '\n')
        file.write(str(index_kick) + '\n')
    file.close()

    with open(file_path, 'r') as file:
        lines = file.readlines()
    file.close()

    water_coord   = []
    OH_coord      = []
    atom_types    = ['O', 'H', 'H']
    atom_types_OH = ['O', 'H']
    for i in range(2,len(lines)):
        line = lines[i].split()
        if line[0] == 'O':
            line1 = lines[i+1].split()
            line2 = lines[i+2].split()
            if line1[0] == 'H' and line2[0] == 'H':
                coordinates = np.array([
                    [float(line[1]),  float(line[2]),  float(line[3])],
                    [float(line1[1]), float(line1[2]), float(line1[3])],
                    [float(line2[1]), float(line2[2]), float(line2[3])]])
                rotation_matrix = create_rotation_matrix() # 生成旋转矩阵和位移向量
                displacement = np.random.normal(0, 2, 3)   # 变换水分子
                transformed_coordinates = apply_transformation(coordinates, rotation_matrix, displacement)
                water_coord.append(i)
                transformed_atoms = [(atom_types[i], transformed_coordinates[i]) for i in range(3)]
                write_xyz_file(file_MOP, transformed_atoms)
                write_xyz_file(file_xyz, transformed_atoms)
            elif line1[0] == 'H' and line2[0] == 'O':
                OH_coord.append(i)
                transformed_coordinates = np.array([
                    [float(line[1]), float(line[2]), float(line[3])],
                    [float(line1[1]), float(line1[2]), float(line1[3])]])
                transformed_atoms = [(atom_types_OH[i], transformed_coordinates[i]) for i in range(2)]
                write_xyz_file(file_MOP, transformed_atoms)
                write_xyz_file(file_xyz, transformed_atoms)
        else:
            pass

    run_mopac(file_MOP)
    file_name = f"cluster_{index_kick}.out"
    file_name_all_opt = 'Mutation_opt_MOPAC_structure.xyz'
    energy = extract_data_from_mopac_output_after_mutation(file_name, index_kick, number_atom, file_name_all_opt)
    if energy == 0:
        os.remove(file_MOP)
        os.remove(file_xyz)
        return energy
    else:
        return energy

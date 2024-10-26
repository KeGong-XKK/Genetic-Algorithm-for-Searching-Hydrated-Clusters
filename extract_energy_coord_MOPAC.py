import re
import numpy as np

def check_atomic_overlaps(coordinates, min_distance=0.8):
    """检查原子间的重叠。

    参数:
    coordinates -- Nx3 的 NumPy 数组，包含 N 个原子的坐标。
    min_distance -- 两原子间的最小允许距离。

    返回:
    重叠的原子对列表。
    """
    num_atoms = coordinates.shape[0]
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = np.linalg.norm(coordinates[i] - coordinates[j])
            if distance < min_distance:
                return 0
            else:
                pass
    return 1


def extract_data_from_mopac_output(file_path,file_name_used_for_xyz,number_atoms,file_name_all_opt):
    """从MOPAC输出文件中提取优化后的坐标和能量
    file_path 为MOPAC优化后得到的out文件
    file_name_used_for_xyz 为用于xyz文件中的cluster编号，1-200
    number_atoms 原子数
    file_name_all_opt 为用于区分100个parents优化后得到的总结构的xyz，还是Child优化后得到的总结构xyz
    """
    energy_pattern = re.compile(r'FINAL HEAT OF FORMATION\s+=\s+([\-\d\.]+)')
    energy = None

    coord_opt = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        energy_match = energy_pattern.search(' '.join(lines))
        if type(energy_match) == str or energy_match is None:
            energy = 0
            return energy
        else:
            line_numer = []
            for i in range(len(lines)):
                if 'CARTESIAN COORDINATES' in lines[i]:
                    line_numer.append(i)
            start_coord_number = int(line_numer[-1]) + 2
            for i in range(start_coord_number, start_coord_number + number_atoms):
                line_coord = lines[i].split()
                coord_opt.append([float(line_coord[2]),float(line_coord[3]),float(line_coord[4])])
            coord_opt = np.array(coord_opt)
            energy = check_atomic_overlaps(coord_opt)
            if energy == 0:
                return energy
            else:
                energy = float(energy_match.group(1))
                with open(file_name_all_opt, 'a') as f:
                        f.write(str(number_atoms) + '\n')
                        f.write(str(file_name_used_for_xyz) + '\n')
                        for i in range(start_coord_number, start_coord_number + number_atoms):
                            line_coord = lines[i].split()
                            f.write(str(line_coord[1]) + ' ' + str(float(line_coord[2])) + ' ' + str(float(line_coord[3])) + ' ' + str(float(line_coord[4])) + '\n')

                filename = f"cluster_{file_name_used_for_xyz}_opt_MOPAC.xyz"
                with open(filename,'w') as f:
                    f.write(str(number_atoms) + '\n')
                    f.write(str(file_name_used_for_xyz) + '\n')
                    for i in range(start_coord_number, start_coord_number + number_atoms):
                        line_coord = lines[i].split()
                        f.write(str(line_coord[1]) + ' ' + str(float(line_coord[2])) + ' ' + str(float(line_coord[3])) + ' ' + str(float(line_coord[4])) + '\n')
                return  energy

def extract_data_from_mopac_output_after_mutation(file_path,file_name_used_for_xyz,number_atoms,file_name_all_opt):
    """从MOPAC输出文件中提取优化后的坐标和能量
    file_path 为MOPAC优化后得到的out文件
    file_name_used_for_xyz 为用于xyz文件中的cluster编号，1-200
    number_atoms 原子数
    file_name_all_opt 为用于区分100个parents优化后得到的总结构的xyz，还是Child优化后得到的总结构xyz
    """
    energy_pattern = re.compile(r'FINAL HEAT OF FORMATION\s+=\s+([\-\d\.]+)')
    energy = None

    coord_opt = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        energy_match = energy_pattern.search(' '.join(lines))
        if type(energy_match) == str or energy_match is None:
            energy = 0
            return energy
        else:


            line_numer = []
            for i in range(len(lines)):
                if 'CARTESIAN COORDINATES' in lines[i]:
                    line_numer.append(i)
            start_coord_number = int(line_numer[-1]) + 2
            for i in range(start_coord_number, start_coord_number + number_atoms):
                line_coord = lines[i].split()
                coord_opt.append([float(line_coord[2]), float(line_coord[3]), float(line_coord[4])])
            coord_opt = np.array(coord_opt)
            energy = check_atomic_overlaps(coord_opt)
            if energy == 0:
                return energy
            else:
                energy = float(energy_match.group(1))
                with open(file_name_all_opt, 'a') as f:
                        f.write(str(number_atoms) + '\n')
                        f.write(str(file_name_used_for_xyz) + '\n')
                        for i in range(start_coord_number, start_coord_number + number_atoms):
                            line_coord = lines[i].split()
                            f.write(str(line_coord[1]) + ' ' + str(float(line_coord[2])) + ' ' + str(float(line_coord[3])) + ' ' + str(float(line_coord[4])) + '\n')

                filename = f"cluster_{file_name_used_for_xyz}_opt_MOPAC.xyz"
                with open(filename,'w') as f:
                    f.write(str(number_atoms) + '\n')
                    f.write(str(file_name_used_for_xyz) + '\n')
                    for i in range(start_coord_number, start_coord_number + number_atoms):
                        line_coord = lines[i].split()
                        f.write(str(line_coord[1]) + ' ' + str(float(line_coord[2])) + ' ' + str(float(line_coord[3])) + ' ' + str(float(line_coord[4])) + '\n')
                return  energy


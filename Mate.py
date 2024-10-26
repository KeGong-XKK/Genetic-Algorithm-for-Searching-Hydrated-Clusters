import random
from MOPAC_opt import write_mopac_input,run_mopac
from extract_energy_coord_MOPAC import extract_data_from_mopac_output

def crossover(parent1, parent2,number_atom,cluster_id):
    # 确定拼接点，这里假设使用团簇的一半
    ''' cluster_id 为形成的child的cluster的编号，从101-200'''

    parent1 = f"cluster_{parent1}_opt_MOPAC.xyz"
    parent2 = f"cluster_{parent2}_opt_MOPAC.xyz"
    with open(parent1,'r') as f:
        lines_parent1 = f.readlines()
    with open(parent2,'r') as f:
        lines_parent2 = f.readlines()

    O_index = []
    for i in range(len(lines_parent1)):
        line = lines_parent1[i].split()
        if line[0] == 'O':
            O_index.append(i)
    O_index = O_index[1:]

    cross_site = random.choice(O_index)
    # 从parent1中获取左半部分
    left_part = lines_parent1[2:cross_site]
    # 从parent2中获取右半部分
    right_part = lines_parent2[cross_site:]

    filename = f"cluster_{cluster_id}.mop"
    with open(filename, 'w') as file:
        file.write("PM6-D3H4 precise CHARGE=-1 SINGLET")
        file.write('\n')
        file.write('\n')
        file.write('\n')
        for line in left_part:
            file.write(line)
        for line in right_part:
            file.write(line)
    run_mopac(filename)
    file_name = f"cluster_{cluster_id}.out"
    file_name_all_opt = 'Child_opt_MOPAC_structure.xyz'
    energy1 = extract_data_from_mopac_output(file_name, cluster_id, number_atom,file_name_all_opt)

    # 从parent1中获取左半部分
    right_part = lines_parent1[cross_site:]
    # 从parent2中获取右半部分
    left_part = lines_parent2[2:cross_site]

    filename = f"cluster_{cluster_id + 1}.mop"
    with open(filename, 'w') as file:
        file.write("PM6-D3H4 precise CHARGE=-1 SINGLET")
        file.write('\n')
        file.write('\n')
        file.write('\n')
        for line in left_part:
            file.write(line)
        for line in right_part:
            file.write(line)
    run_mopac(filename)
    file_name = f"cluster_{cluster_id + 1}.out"
    file_name_all_opt = 'Child_opt_MOPAC_structure.xyz'
    energy2 = extract_data_from_mopac_output(file_name, cluster_id + 1, number_atom, file_name_all_opt)
    if energy1 == 0 or energy2 == 0:
        return energy1,energy2
    else:
        filename = f"cluster_{cluster_id}.xyz"
        left_part = lines_parent1[2:cross_site]
        right_part = lines_parent2[cross_site:]
        with open(filename, 'w') as file:
            file.write(str(number_atom) + '\n')
            file.write(str(cluster_id)  + '\n')
            for line in left_part:
                file.write(line)
            for line in right_part:
                file.write(line)

        filename = f"cluster_{cluster_id + 1}.xyz"
        right_part = lines_parent1[cross_site:]
        left_part = lines_parent2[2:cross_site]
        with open(filename, 'w') as file:
            file.write(str(number_atom) + '\n')
            file.write(str(cluster_id + 1) + '\n')
            for line in left_part:
                file.write(line)
            for line in right_part:
                file.write(line)
        return energy1,energy2





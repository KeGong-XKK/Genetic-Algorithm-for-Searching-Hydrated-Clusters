import numpy as np
import os

# 路径前缀，假设脚本与文件夹在同一目录下
base_path = './'

# 存储所有能量值的列表
all_energies    = []
unsorted_energy = []
# 遍历文件夹从1到30
for i in range(1, 31):
    folder_name = f"{base_path}{i}_generation"
    file_path = os.path.join(folder_name, 'all_energies.data')

    # 检查文件是否存在
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for id in range(100):
            energy = float(lines[id].split()[0])
            formatted_energy = round(energy, 1)

            all_energies.append(formatted_energy)
            unsorted_energy.append(energy)
            with open('all_energies_GA.data','a') as f:
                f.write(str(formatted_energy) + ' ' + str(i) + ' ' + str(id+1) + '\n')
            f.close()
            with open('all_energies_GA_un_formatted.data','a') as f:
                f.write(str(energy) + ' ' + str(i) + ' ' + str(id+1) + '\n')
            f.close()
# 使用numpy去重并排序

unique_sorted_energies, indices = np.unique(all_energies, return_index=True)
print(unique_sorted_energies)
with open('all_energies_GA.data', 'r') as f:
    lines = f.readlines()
f.close()

pwd = os.getcwd()
os.system('mkdir GA_results')
new_path = pwd + '/GA_results'
os.chdir(new_path)
for i in range(len(indices)):
    generation = lines[indices[i]].split()[1]
    individual = lines[indices[i]].split()[2]
    source_file = 'cp ../' + f'{generation}_generation' + '/' + f'cluster_{individual}_opt_MOPAC.xyz' + ' ' + './' + f'{i+1}.xyz'
    os.system(source_file)
    with open(f'{i+1}.xyz','r') as f:
        lines_num = f.readlines()
        number_atom = int(lines_num[0].split()[0])
    with open('all_sorted.xyz','a') as f:
        f.write(str(number_atom) + '\n')
        f.write(str(i+1) + '\n')
        for q in range(2,len(lines_num)):
            f.write(lines_num[q])

    with open('all_energies_GA_sorted.data', 'a') as f:
        f.write(str(unsorted_energy[indices[i]]) + ' ' + str(i+1)  + '\n')
    f.close()


import os

import numpy as np
import random
import subprocess
from read_clusters import read_clusters_from_xyz
from MOPAC_opt import write_mopac_input,run_mopac
from extract_energy_coord_MOPAC import extract_data_from_mopac_output
from Mate import crossover
import kick

def run_genmer():
    process = subprocess.Popen(['genmer'],
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True)
    output, errors = process.communicate(input='\n')  # 发送回车键

    # 输出处理结果
    print("Output:", output)
    if errors:
        print("Errors:", errors)

def calculate_fitness(energies):
    """计算归一化倒置适应度"""
    max_energy = np.max(energies)
    fitnesses = np.exp(max_energy - energies)
    return fitnesses

def select_parents(fitnesses):
    """根据适应度加权随机选择两个父代"""
    return np.random.choice(len(fitnesses), 2, p=fitnesses / np.sum(fitnesses), replace=False)

'''               运行前先修改genmer.ini文件中的分子数目                   '''
number_H2O      = 5
number_OH_lizi  = 1
number_atom     = number_H2O*3 + number_OH_lizi*2
number_clusters = 100
n_GA_iter       = 30
number_molcule = number_H2O + number_OH_lizi

for gi in range(n_GA_iter):
    path = os.getcwd()  # 当前在最外层
    new_path = path + '/' + str(gi+1) + '_generation'
    os.mkdir(new_path)
    while True:
        '''                       生成100个团簇结构                              '''
        if gi == 0:
            run_genmer()
        '''                      读取traj.xyz中的团簇                            '''
        clusters = read_clusters_from_xyz("traj.xyz",number_atom)

        '''               为每个团簇生成MOPAC输入文件并运行优化                      '''
        for index, cluster in enumerate(clusters):
            write_mopac_input(cluster, index,number_atom)
            input_file = f"cluster_{index + 1}.mop"
            run_mopac(input_file)
            print(f"Optimization for Cluster {index + 1} completed.")

        '''                 统计100个团簇的MOPAC输出文件                           '''
        all_energies = []
        for i in range(1, number_clusters+1):
            file_name = f"cluster_{i}.out"
            file_name_all_opt = 'Parents_opt_MOPAC_structure.xyz'
            energy = extract_data_from_mopac_output(file_name,i,number_atom,file_name_all_opt)
            all_energies.append(energy)
            with open('Parents_opt_MOPAC_energy.data','a') as f:
                f.write(str(i) + ' ' + str(energy) + '\n')
            f.close()
        if 0 in all_energies:
            os.system('rm -rf *.data *.xyz *.mop *.out *.arc')
        else:
            break
    # 计算100个团簇的适应度及筛选合适的parent进行交叉##########################################
    fitness = calculate_fitness(all_energies)
    # 每次选择两个高适应度结构进行交叉，随后使用MOPAC优化，如果报错则重新选择组合优化#################
    for i in range(101,200,2):
        while True:
            parents = select_parents(fitness)
            energy1, energy2 = crossover(int(parents[0])+1,int(parents[1])+1,number_atom,i)
            if energy1 == 0 or energy2 == 0:
                print('Choose Again')
            else:
                all_energies.append(energy1)
                all_energies.append(energy2)
                with open('Child_opt_MOPAC_energy.data', 'a') as f:
                    f.write(str(i) + ' ' + str(energy1) + '\n')
                    f.write(str(i+1) + ' ' + str(energy2) + '\n')
                f.close()
                break
    ##############################################################################
    #对parents及child共同组成的energy进行排序，选择能量较高前30%作为后续变异的选择
    # 获取排序后的索引（np.argsort 返回的是从小到大的索引）
    sorted_indices = np.argsort(all_energies)
    start_index = int(len(all_energies) * 0.7)
    # 从排序后的列表中选取后30%的索引（这里取最低的30%）
    selected_indices = sorted_indices[start_index:]

    #接下来对选择的百分之三十的结构进行变异处理
    mutation_index_start = 201
    mutation_indivi      = []
    for index_kick in selected_indices:
        mutation_indivi.append(int(index_kick)+1)
        while True:
            file_name_kick = f"cluster_{int(index_kick) + 1}_opt_MOPAC.xyz"
            energy = kick.mutation(file_name_kick,mutation_index_start,number_atom)
            if energy == 0:
                print('Mutation Again')

            else:
                all_energies.append(energy)
                with open('Mutation_opt_MOPAC_energy.data', 'a') as f:
                    f.write(str(mutation_index_start) + ' ' + str(energy) + '\n')
                f.close()
                break
        mutation_index_start += 1

    '''对现在所有的结构能量进行排序，找到能量最稳定的100个构型合并为traj.xyz并放到2 generation文件夹中'''
    # 获取排序后的索引（np.argsort 返回的是从小到大的索引）
    sorted_indices = np.argsort(all_energies)
    with open('all_energies.data','w') as f:
        for i in range(len(all_energies)):
            f.write(str(all_energies[i]) + '\n')
    f.close()

    selected_indices = sorted_indices[:100]
    next_generation_index_start = 1
    for index_kick in selected_indices:
        file_name_next_indivi = f"cluster_{int(index_kick) + 1}_opt_MOPAC.xyz"
        with open('next_indivi.data','a') as f:
            f.write(str(next_generation_index_start) + ' ' + str(int(index_kick)+1) + '\n')
        f.close()

        with open(file_name_next_indivi,'r') as file:
            lines = file.readlines()
        file.close()

        with open('traj_next_generation.txt', 'a') as f:
            f.write(str(number_atom) + '\n')
            f.write(str(next_generation_index_start) + '\n')
            for i in range(2,len(lines)):
                f.write(lines[i])
        f.close()
        next_generation_index_start += 1
    command = 'mv *.data *.xyz *.mop *.out *.arc ' + str(gi+1) + '_generation/'
    os.system(command)
    os.system('mv traj_next_generation.txt traj.xyz')


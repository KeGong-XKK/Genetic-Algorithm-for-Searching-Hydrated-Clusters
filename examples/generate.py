import numpy as np
import random

def read_xyz(filename):
    """读取XYZ文件并返回原子类型和坐标"""
    with open(filename, 'r') as file:
        lines = file.readlines()[2:]  # 跳过前两行
        atoms = []
        for line in lines:
            parts = line.split()
            atom_type = parts[0]
            coordinates = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            atoms.append((atom_type, coordinates))
    return atoms

def generate_cluster(h2o_atoms, oh_atoms, num_h2o, R=10.0):
    """生成一个包含多个水分子和一个OH离子的团簇"""
    cluster = []
    # 添加OH离子
    cluster.append(oh_atoms[0])  # 假设OH离子文件中只有一个OH离子
    # 随机添加水分子
    for _ in range(num_h2o):
        # 随机选择一个水分子
        selected_h2o = random.choice(h2o_atoms)
        # 为水分子生成随机位移
        displacement = np.random.uniform(-R, R, 3)
        # 添加位移后的水分子到团簇中
        new_position = selected_h2o[1] + displacement
        cluster.append((selected_h2o[0], new_position))
    return cluster

def write_clusters_to_xyz(clusters, base_filename):
    """将所有团簇写入XYZ文件"""
    for i, cluster in enumerate(clusters):
        filename = f"{base_filename}_{i+1}.xyz"
        with open(filename, 'w') as file:
            file.write(f"{len(cluster)}\n")
            file.write("Generated cluster\n")
            for atom in cluster:
                file.write(f"{atom[0]} {atom[1][0]} {atom[1][1]} {atom[1][2]}\n")

# 读取水分子和OH离子的坐标
h2o_atoms = read_xyz("H2O.xyz")
oh_atoms = read_xyz("OHlizi.xyz")

# 生成100个团簇模型
clusters = [generate_cluster(h2o_atoms, oh_atoms, 4) for _ in range(100)]

# 保存生成的团簇到XYZ文件
write_clusters_to_xyz(clusters, "cluster")
def read_clusters_from_xyz(filename,num_atoms):
    """从XYZ文件中读取所有团簇"""
    with open(filename, 'r') as file:
        lines = file.readlines()
    clusters = []
    i = 0
    while i < len(lines):
        if lines[i].strip().isdigit() and int(lines[i].strip()) == num_atoms:
            cluster = lines[i + 2:i + 2 + num_atoms]  # 跳过原子数行和注释行
            clusters.append(cluster)
            i += num_atoms + 2
        else:
            i += 1
    return clusters
def write_mopac_input(cluster, index,num_atom):
    """为每个团簇写入一个MOPAC输入文件"""
    filename = f"cluster_{index + 1}.mop"
    with open(filename, 'w') as file:
        file.write("PM6-D3H4 precise CHARGE=-1 SINGLET")
        file.write('\n')
        file.write('\n')
        file.write('\n')
        for line in cluster:
            file.write(line)

    filename = f"cluster_{index + 1}.xyz"
    with open(filename,'w') as file:
        file.write(str(num_atom) + '\n')
        file.write(str(index+1)  + '\n')
        for line in cluster:
            file.write(line)


def run_mopac(input_file):
    """调用MOPAC进行优化"""
    import subprocess
    command = f"MOPAC2016.exe {input_file}"
    subprocess.run(command, shell=True)
import subprocess
import os
import numpy as np

def runc():
    
    subprocess.run(["g++", "rapid.cpp", "-o", "rapid"])
    subprocess.run(["./rapid"])

def read_only_lines(filepath, start, finish):
    f = open(filepath)

    for ii,line in enumerate(f):
        
        if ii>=start and ii<finish:
            yield line
        
        elif ii>=finish:
            return
        

def convert_str_to_num(str):
    num = ''
    for i in str:

        if i=='A':
            num+='0'

        elif i=='C':
            num+='1'

        elif i=='G':
            num+='2'

        elif i=='T':
            num+='3'

    return num


def check_memory(file_path):
    max_size = 1500000
    file_info = os.stat(file_path)
    file_size = file_info.st_size
    
    if file_size<max_size:
        return 0
    return 1


def run():
    file_path1 = 'C:\Users\yashn\Desktop\PIM\Fin\output\seq1file1.txt'
    file_path2 = 'C:\Users\yashn\Desktop\PIM\Fin\output\seq1file2.txt'
    k = len(enumerate(open(file_path1)))
    
    for i in range(0,k,8):
        str1 = read_only_lines(file_path1,i,i+10)
        str2 = read_only_lines(file_path2,i,i+10)
        
        if check_memory(file_path1) or check_memory(file_path2):
            print("File Size exceed Iteration", i)
        
        file_path1 = "seq1file2_small"
        file_path2 = "seq1file2_small"
        
        with open(file_path1, "w") as file:
            file.write(str1)
        with open(file_path2, "w") as file:
            file.write(str2)
        
        
        runc()


if __name__ == "__main__":
    run()

import subprocess
import sys

assert len(sys.argv) > 1

for cpp_file in sys.argv[1:]:
    cpp_short = cpp_file[:-4].split('-')[0]
    cmd_compile = ['g++', '-O3', '-std=c++2a', cpp_file, '-o', cpp_short]
    subprocess.run(cmd_compile)

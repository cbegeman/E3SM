import shutil, time, sys, os, subprocess, fileinput

with open('db_debug.log') as f:
    lines = f.readlines()
    for line in lines:
       if 'SSDB' in line:
           SSDB = line.split('=')[1]
os.environ['SSDB'] = SSDB
print(f'SSDB={os.environ.get("SSDB")}')

with fileinput.FileInput(f'environ.sh', inplace=True,
               backup='.bak') as f:
    for line in f:
        print(line.replace('export SSDB=None', f'export SSDB={SSDB}'), end='')

#!/usr/bin/python
import os, subprocess, sys
subprocess.call(['python', 'virtualenv.py', 'hase'])

f=open('requirements.txt')
p=f.readlines()
bin = 'bin'
for i in [j.split('\n')[0] for j in p ]:
	subprocess.call([os.path.join('hase', bin, 'pip'), 'install', i])
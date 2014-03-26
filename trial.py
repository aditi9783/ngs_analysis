#!/usr/bin/python

import subprocess

subprocess.call('ls');
subprocess.call(['bash', '-c', 'module load GNU/4.4.5']);


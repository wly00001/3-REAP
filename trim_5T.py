

############ add by Luyang Wang, the code is originally from Dinghai Zheng, https://github.com/DinghaiZ. I made a few changes for the purpose of 3'REAP.


import argparse
##arg setting parts###
parser = argparse.ArgumentParser(description="3'READS Python")
parser.add_argument("--rawfastq_dir", action="store", dest='rawfastq_dir',default='', metavar='<rawfastq_dir>', help="define your rawfastq_dir")
parser.add_argument("--project_dir", action="store", dest='project_dir',default='', metavar='<project_dir>', help="define your project_dir")


args=parser.parse_args()

rawfastq_dir = args.rawfastq_dir
project_dir = args.project_dir

############

import os
import glob
import re
import csv
import sys
import math
import subprocess
from subprocess import check_output
import argparse
import unittest
from Bio.Seq import Seq

# set up folders, parameters etc.
os.chdir(project_dir) # This is the working directory

# directory for results including read number stats.
result_dir = project_dir + '/result'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
   
################################################################################

class Fastq(object):
    """
    Fastq record
    """
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
    
    def __str__(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, 
            '+%s' % self.name, self.qual])
    
    def as_simple(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, '+', self.qual])
    
    def get_name(self):
        return self.name
    
    def clean_name(self):
        self.name = self.name.split(" ")[0]
    
    def get_seq(self):
        return self.seq
    
    def get_qual(self):
        return self.qual
    
    def get_ascii(self):
        return map(ord, self.qual)
    
    def get_length(self):
        """Return length"""
        return len(self.seq)

    def get_num_N(self, first=28):
        """Return number of Ns"""
        return self.seq.count('N', 0, first)

    def set_name(self, name):
        """Reset name"""
        self.name = name
    
    def remove_tail_N(self):
        """Remove N at the end"""
        i = len(self.seq) - 1
        while i >= 0:
            if self.seq[i] == 'N':
                i -= 1
            else:
                break
        self.seq = self.seq[:i+1]
        self.qual = self.qual[:i+1]

################################################################################
def count_fastq_in_folder(infolder, outfolder,  
                          pattern = '', outfile = 'read_num.csv'):
    """
    Counts the number of fastq records in fastq files
    """
    from subprocess import check_output
    import glob
    
    pattern1 = '*' + pattern + '.fastq'  
    fls = glob.glob(os.path.join(infolder, pattern1))
    if fls: # not empty
        cmd = 'wc -l ' + os.path.join(infolder, pattern1)
        cmd_out = check_output(cmd, shell=True)
        cmd_out = cmd_out.decode() ### modified by lw
        cmd_out = cmd_out.split('\n')[:-1] # remove  ''
        cmd_out = [x for x in cmd_out if not x.endswith(' total')] # remove 'total'
        rec_nums = [str(int(int(x.split(' /')[0])/4)) for x in cmd_out]
        #file_names = [x.split('../')[1].split('/')[-1].split('.')[0] for x in cmd_out]
        file_names = [x.split('/')[-1].split('.')[0] for x in cmd_out]
        records = zip(file_names, rec_nums)
        
        with open(os.path.join(outfolder, outfile), 'w') as f:
            f.write('Sample,Count\n')
            for rec in records:
                f.write('%s,%s\n' % (rec[0],rec[1])) ### modified by lw

################################################################################
def fastq_trim_Ts(infile, outfile, random_NT_len=6):

    import re
    pattern = re.compile('([ATCGN]{%d})(T*)'% random_NT_len)#precompile
    with open(outfile, 'w') as outhandle:
        c=0 # counting
        for fq in reader_fastq(infile):
            seq = fq.get_seq()
            qual = fq.get_qual()
            
            match1 = re.match(pattern,seq)
            seq = seq[match1.end():]
            qual = qual[match1.end():]
            T_length1 = len(match1.groups()[1])
            read_name = fq.get_name()
            if len(seq) >= 18:
                c+=1 # counting
                outhandle.write('\n'.join([read_name, seq, '+', qual]) + '\n')






### count number of reads after fastx_clipper
count_fastq_in_folder(infolder = rawfastq_dir, pattern = 'clipped', 
                      outfolder = result_dir, outfile = 'clipped_read_num.csv')


### trim 5'Ts 
#fastq_files = !ls {rawfastq_dir + '/' + '*clipped.fastq'}
for fastq_file in sorted(glob.glob(os.path.join(rawfastq_dir, '*clipped.fastq'))):
    print ("Trimming 5'Ts in ", fastq_file)
    fastq_trim_Ts(infile = fastq_file, outfile = fastq_file.replace('.fastq', 
    '.trimmed.fastq'), random_NT_len=0)

import shlex
import subprocess
import os
import signal

from distutils.spawn import find_executable


class bam_mpileup:
    '''
    Use samtools via subprocess and return an iterable
    object.
    '''

    def __init__(self, bam, fasta, q=20, Q=20,
                 samtools_bin='samtools', regions=[]):
        samtools_bin = find_executable(samtools_bin)
        if samtools_bin is None:
            raise Exception((
                'Please install samtools or provide samtools path in the '
                'appropriate argument'))
        samtools_view = "%s view -T %s -u %s" % (samtools_bin, fasta, bam)
        samtools_view = shlex.split(samtools_view) + regions
        samtools_pile = "%s mpileup -f %s -q %i -Q %i -" % (
            samtools_bin, fasta, q, Q)
        samtools_pile = shlex.split(samtools_pile)
        self.proc1 = subprocess.Popen(
            samtools_view, stdout=subprocess.PIPE, bufsize=4096)
        self.proc2 = subprocess.Popen(samtools_pile, stdin=self.proc1.stdout,
                                      stdout=subprocess.PIPE, bufsize=4096)

    def __iter__(self):
        while True:
            try:
                yield next(self.proc2.stdout).decode('utf-8')
            except StopIteration:
                break

    def close(self):
        self.proc1.stdout.close()
        # self.proc2.communicate()
        kill_subproc(self.proc1)
        kill_subproc(self.proc2)
        self.proc2.stdout.close()


class indexed_pileup:
    '''
    Use tabix via subprocess to slice the pileup data and return
    an iterable object
    '''

    def __init__(self, pileup, tabix_bin='tabix', regions=[]):
        tabix_bin = find_executable(tabix_bin)
        if tabix_bin is None:
            raise Exception((
                'Please install tabix or provide tabix path in the '
                'appropriate argument'))
        tabix = "%s %s" % (tabix_bin, pileup)
        tabix = shlex.split(tabix) + regions
        self.proc = subprocess.Popen(tabix, stdout=subprocess.PIPE,
                                     bufsize=4096)

    def __iter__(self):
        for line in self.proc.stdout:
            yield line.decode('utf-8')

    def close(self):
        kill_subproc(self.proc)
        self.proc.stdout.close()


def tabix_seqz(file_name, tabix_bin='tabix', seq=1, begin=2, end=2, skip=1):
    '''
    Index a seqz file with tabix
    '''
    tabix_bin = find_executable(tabix_bin)
    if tabix_bin is None:
        raise Exception((
            'Please install tabix or provide tabix path in the '
            'appropriate argument'))
    tabix = "%s -f -s %s -b %s -e %s -S %s %s" % (
        tabix_bin, seq, begin, end, skip, file_name)
    tabix = shlex.split(tabix)
    proc1 = subprocess.Popen(tabix)
    proc1.communicate()


def program_version(program):
    '''
    Parse tabix or samtools help message in attempt to
    retrieve the software version: return format: [major, minor, *]
    '''
    program_bin = find_executable(program)
    if program_bin is None:
        raise Exception((
            'Please install %s or provide %s path in the '
            'appropriate argument') % (program, program))
    else:
        proc1 = subprocess.Popen([program_bin], stderr=subprocess.PIPE)
        for line in proc1.stderr:
            if line.startswith(b'Version:'):
                return map(int, line.rstrip().split(b' ')[1].split(b'.'))


def kill_subproc(proc):
    try:
        proc.kill()
        proc.wait()
    except AttributeError:
        os.kill(proc.pid, signal.SIGKILL)
        os.waitpid(proc.pid, 0)

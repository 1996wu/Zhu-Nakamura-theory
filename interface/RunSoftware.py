import os
import shutil
import subprocess
import time

from .prepare import soft, soft_path, EnterDir
from .prepare import output_suffix, input_suffix, input_prefix

__all__ = ["software_running"]


def temp_dic():
    # GaussianTemp OrcaTemp MolproTemp
    tempDic = soft + "Temp"
    try:
        os.mkdir(tempDic)
    except FileExistsError:
        pass


# soft = soft_type
temp_dic()
if soft == "Gaussian":
    shutil.copy('gauss.chk', './GaussianTemp')


class QuRun:
    """
    This is class for running quantum chemistry software
    """

    def __init__(self):
        self.soft = soft
        self.path = None
        self.temDic = soft + "Temp"
        self.filename = None
        self.output = None
        self.command = None
        self.suffix = output_suffix[self.soft]

    @staticmethod
    def current_time():
        return time.asctime(time.localtime(time.time()))

    def prepare(self):
        command_type = {"Gaussian": "which g16", "Orca": "which orca", "Molpro": "which molpro"}
        # GaussianPath
        # /home/wzb/software/g16/g16
        # OrcaPath
        # /home/wzb/software/orca_4_2_1_linux_x86-64_shared_openmpi314/orca
        try:
            self.path = soft_path
        except NameError:
            _path = subprocess.run(command_type[self.soft], shell=True, stdout=subprocess.PIPE, text=True)
            self.path = _path.stdout.replace("\n", "")

    def input_file(self, file=None):
        if file:
            self.filename = file
        else:
            self.filename = input_prefix[self.soft] + "." + input_suffix[self.soft]

    def output_file(self):
        prefix = self.filename.split(".")[0]
        self.output = prefix + "." + self.suffix

    def running_command(self):
        if self.soft == "Gaussian":
            # g16 gauss.gjf => gauss.out
            self.command = self.path + " " + self.filename
        elif self.soft == "Orca":
            # orcaPath  orca.inp &> orca.out
            self.command = self.path + " " + self.filename + " &> " + self.output
        elif self.soft == "Molpro":
            # save wfu file in current dir and forbid backup file
            self.command = self.path + " " + " -W ./  " + " --backup 1 " + self.filename
        else:
            # I do not know...how to run molpro molcas BDF
            pass

    def running(self):
        print("The %s begin running at %s" % (self.soft, self.current_time()), flush=True)
        shutil.copy(self.filename, self.temDic)
        with EnterDir(self.temDic):
            proc = subprocess.Popen(self.command, shell=True)
            proc.wait()
            shutil.copy(self.output, "../")
        print("The %s ended at %s" % (self.soft, self.current_time()), flush=True)

    def worker(self, inputfile=None):
        self.prepare()
        self.input_file(inputfile)
        self.output_file()
        self.running_command()
        self.running()


def software_running(filename=None):
    sr = QuRun()
    sr.worker(filename)


if __name__ == "__main__":
    software_running()
from abics.applications.latgas_abinitio_interface import aenet
import numpy as np
import os, pathlib, shutil, subprocess
import time


class aenet_trainer:
    def __init__(
        self,
        structures,
        energies,
        generate_inputdir,
        train_inputdir,
        predict_inputdir,
        generate_exe,
        train_exe,
    ):
        self.structures = structures
        self.energies = energies
        self.generate_inputdir = generate_inputdir
        self.train_inputdir = train_inputdir
        self.predict_inputdir = predict_inputdir
        self.generate_exe = generate_exe
        self.train_exe = train_exe
        assert len(self.structures) == len(self.energies)
        self.numdata = len(self.structures)
        self.is_prepared = False
        self.is_trained = False
        self.generate_outputdir = None

    def prepare(self, latgas_mode = True, st_dir = "aenetXSF"):
        rootdir = os.getcwd()
        xsfdir = os.path.join(rootdir, st_dir)
        
        # prepare XSF files for aenet
        os.makedirs(xsfdir, exist_ok=True)
        os.chdir(xsfdir)
        xsfdir = os.getcwd()
        if latgas_mode:
            for i, st in enumerate(self.structures):
                xsf_string = aenet.to_XSF(st, write_force_zero=False)
                xsf_string = (
                    "# total energy = {} eV\n\n".format(self.energies[i]) + xsf_string
                )
                with open("structure.{}.xsf".format(i), "w") as fi:
                    fi.write(xsf_string)
        else:
            for i, st in enumerate(self.structures):
                xsf_string = aenet.to_XSF(st, write_force_zero=False)
                xsf_string = (
                    "# total energy = {} eV\n\n".format(self.energies[i]) + xsf_string
                )
                with open("structure.{}.xsf".format(i), "w") as fi:
                    fi.write(xsf_string)

        os.chdir(rootdir)

    def generate_run(self, xsfdir="aenetXSF", generate_dir="generate"):
        # prepare generate
        xsfdir = str(pathlib.Path(xsfdir).resolve())
        if os.path.exists(generate_dir):
            shutil.rmtree(generate_dir)
        shutil.copytree(self.generate_inputdir, generate_dir)
        os.chdir(generate_dir)
        with open("generate.in.head", "r") as fi:
            generate_head = fi.read()
            xsf_paths = [
                os.path.join(xsfdir, "structure.{}.xsf".format(i))
                for i in range(self.numdata)
            ]
            generate = (
                generate_head
                + "\n"
                + "FILES\n"
                + str(self.numdata)
                + "\n"
                + "\n".join(xsf_paths)
                + "\n"
            )
            with open("generate.in", "w") as fi_in:
                fi_in.write(generate)

        command = self.generate_exe + " generate.in"
        with open(os.path.join(os.getcwd(), "stdout"), "w") as fi:
            #subprocess.run(
            self.gen_proc = subprocess.Popen(
                command, shell=True, stdout=fi, stderr=subprocess.STDOUT,#, check=True
                )
        self.generate_outputdir = os.getcwd()
        os.chdir(pathlib.Path(os.getcwd()).parent)
        #self.is_prepared = True
        
    def generate_wait(self):
        self.gen_proc.wait()
        self.is_prepared = True

    def train(self, train_dir = "train"):
        try:
            assert self.is_prepared
        except AssertionError as e:
            e.args += "you have to prepare the trainer before training!"
        if os.path.exists(train_dir):
            shutil.rmtree(train_dir)
        shutil.copytree(self.train_inputdir, train_dir)
        os.chdir(train_dir)
        os.rename(
            os.path.join(self.generate_outputdir, "aenet.train"),
            os.path.join(os.getcwd(), "aenet.train"),
        )
        command = self.train_exe + " train.in"

        while True:
            # Repeat until test set error begins to rise
            with open(os.path.join(os.getcwd(), "stdout"), "w") as fi:
                subprocess.run(
                    command, shell=True, stdout=fi, stderr=subprocess.STDOUT, check=True
                )
            with open("stdout", "r") as trainout:
                fullout = trainout.readlines()
                epoch_data = []
                for li in fullout:
                    if "<" in li:
                        epoch_data.append(li)
            with open("epochdat", "w") as epochdatfi:
                epochdat_str = "".join(epoch_data[1:]).replace("<", "")
                epochdatfi.write(epochdat_str)
            epoch_dat_arr = np.loadtxt("epochdat")
            # Find epoch id with minimum test set RMSE
            testRMSE = epoch_dat_arr[:, 4]
            minID = np.argmin(testRMSE)
            if minID == 0:
                minID = np.argmin(testRMSE[1:]) + 1
            num_epoch = len(testRMSE)
            if minID < num_epoch*0.7:  # this "0.7" is a heuristic
                break

        print("Best fit at epoch ID ", minID)
        self.train_outputdir = os.getcwd()
        self.train_minID = minID
        os.chdir(pathlib.Path(os.getcwd()).parent)
        self.is_trained = True

    def new_baseinput(self, baseinput_dir):
        try:
            assert self.is_trained
        except AssertionError as e:
            e.args += "you have to train before getting results!"

        # Some filesystems may delay making a directory due to cache
        # especially when mkdir just after rmdir, and hence 
        # we should make sure that the old directory is removed and the new one is made.
        # Since `os.rename` is an atomic operation,
        # `baseinput_dir` is removed after `os.rename`.
        if os.path.exists(baseinput_dir):
            os.rename(baseinput_dir, baseinput_dir + "_temporary")
            shutil.rmtree(baseinput_dir + "_temporary")
        os.makedirs(baseinput_dir, exist_ok=False)
        while not os.path.exists(baseinput_dir):
            time.sleep(0.1)

        iflg = False
        for name in ["predict.in", "in.lammps"]:
            if os.path.isfile(os.path.join(self.predict_inputdir), name):
                iflg = True
                shutil.copyfile(
                    os.path.join(self.predict_inputdir, name),
                    os.path.join(baseinput_dir, name),
                )
        if iflg is False:
            print("Warning: predict.in or in.lammps should be in the predict directory.")


        NNPid_str = "{:05d}".format(self.train_minID)
        NNPfiles = [fi for fi in os.listdir(self.train_outputdir) if NNPid_str in fi]
        for fi in NNPfiles:
            shutil.copyfile(
                os.path.join(self.train_outputdir, fi),
                os.path.join(baseinput_dir, fi),
            )
            # os.rename is guaranteed to be atomic
            os.rename(
                os.path.join(baseinput_dir, fi),
                os.path.join(baseinput_dir, fi[:-6]),
            )

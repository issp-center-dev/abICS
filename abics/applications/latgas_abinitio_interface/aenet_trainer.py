from abics.applications.latgas_abinitio_interface import aenet
import numpy as np
import os, pathlib, shutil, subprocess


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

    def prepare(self, latgas_mode = True):
        # prepare XSF files for aenet
        os.makedirs("aenetXSF", exist_ok=True)
        os.chdir("aenetXSF")
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

        xsfdir = os.getcwd()
        os.chdir(pathlib.Path(os.getcwd()).parent)
        # prepare generate
        if os.path.exists("generate"):
            shutil.rmtree("generate")
        shutil.copytree(self.generate_inputdir, "generate")
        os.chdir("generate")
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
            subprocess.run(
                command, shell=True, stdout=fi, stderr=subprocess.STDOUT, check=True
            )
        self.generate_outputdir = os.getcwd()
        os.chdir(pathlib.Path(os.getcwd()).parent)
        self.is_prepared = True

    def train(self):
        try:
            assert self.is_prepared
        except AssertionError as e:
            e.args += "you have to prepare the trainer before training!"
        if os.path.exists("train"):
            shutil.rmtree("train")
        shutil.copytree(self.train_inputdir, "train")
        os.chdir("train")
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
        os.makedirs(baseinput_dir, exist_ok=True)
        shutil.copyfile(
            os.path.join(self.predict_inputdir, "predict.in"),
            os.path.join(baseinput_dir, "predict.in"),
        )
        NNPid_str = "{:05d}".format(self.train_minID)
        NNPfiles = [fi for fi in os.listdir(self.train_outputdir) if NNPid_str in fi]
        for fi in NNPfiles:
            shutil.copyfile(
                os.path.join(self.train_outputdir, fi),
                os.path.join(baseinput_dir, fi[:-6]),
            )

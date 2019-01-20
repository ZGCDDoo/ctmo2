import os
import shutil
import unittest
import numpy as np
import subprocess
import json

# to run only one test: python -m unittest -v test.test_integration.TestIntegration.test_dca2x2


class TestIntegration(unittest.TestCase):

    mpi_cmd = "mpirun -np 4 ctmo params 1"
    mpi_dca_cmd = "mpirun -np 4 ctmo_DCA params 1"

    @classmethod
    def setUpClass(cls):
        print("\nIn TestIntegration.\n")
        # shutil.rmtree("./tmp", ignore_errors=True)

    def test_holstein_dmft(self):
        base_path = os.getcwd()
        holstein_path = os.path.join(base_path, "test/Simulations/holstein_DMFT_T0.3")
        tmp_path = os.path.join(base_path, "tmp/holstein_DMFT_T0.3")
        shutil.rmtree(tmp_path, ignore_errors=True)
        shutil.copytree(holstein_path, tmp_path)
        holstein_path = os.path.join(holstein_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_cmd, shell=True)
        good_result = np.loadtxt(
            os.path.join(tmp_path, "GoodReza/Region1_LocalGFn.dat")
        )
        result = np.loadtxt("greenUp1.dat")
        os.chdir(base_path)
        # shutil.rmtree(tmp_path)

        len_result = result.shape[0]
        len_high = len_result

        # test all the  frequencies
        np.testing.assert_allclose(
            result[:len_result],
            good_result[:len_result],
            rtol=1e-3,
            atol=2e-3,
            verbose=True,
        )

        # Test the first few and last frequencies - imaginary part
        np.testing.assert_allclose(
            result[0:5, 2], good_result[0:5, 2], rtol=1e-3, atol=5e-4
        )
        np.testing.assert_allclose(
            result[200:205, 2], good_result[200:205, 2], rtol=1e-3, atol=4e-4
        )

        np.testing.assert_allclose(
            result[500:520, 2], good_result[500:520, 2], rtol=1e-3, atol=3e-4
        )

    def test_triangle2x2(self):
        base_path = os.getcwd()
        triangle_path = os.path.join(
            base_path, "test/Simulations/2x2Triangle_tp04_b5_U3"
        )
        tmp_path = os.path.join(base_path, "tmp/2x2Triangle_tp04_b5_U3")
        shutil.rmtree(tmp_path, ignore_errors=True)

        shutil.copytree(triangle_path, tmp_path)
        triangle_path = os.path.join(triangle_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_cmd, shell=True)
        good_result = np.loadtxt(os.path.join(tmp_path, "GoodNew/greenUp1.dat"))
        good_result_old = np.loadtxt(os.path.join(tmp_path, "GoodOld/greenUp1.dat"))
        result = np.loadtxt("greenUp1.dat")
        os.chdir(base_path)
        # shutil.rmtree(tmp_path)

        np.testing.assert_allclose(result, good_result, rtol=5e-4, atol=4e-4)

        # The imaginary parts should be even closer
        np.testing.assert_allclose(
            result[:, 2], good_result[:, 2], rtol=5e-4, atol=1e-4
        )

        len_result = result.shape[0]

        # only test the first green function, because order is not the same.
        np.testing.assert_allclose(
            result[:, 0:3], good_result_old[:len_result, 0:3], rtol=5e-4, atol=5e-4
        )

        # now test the observables
        with open(os.path.join(tmp_path, "Obs.json")) as fin:
            jj_result = json.load(fin)

        with open(os.path.join(tmp_path, "GoodOld/Obs.json")) as fin:
            jj_good = json.load(fin)

        keys = ["docc", "n", "sign", "k"]
        for key in keys:
            np.testing.assert_allclose(
                jj_result[key][0], jj_good[key][0], rtol=1e-3, atol=1e-3
            )

    def test_dca2x2(self):
        base_path = os.getcwd()
        dca2x2_path = os.path.join(base_path, "test/Simulations/2x2_DCA_b5_U3")
        tmp_path = os.path.join(base_path, "tmp/2x2_DCA_b5_U3")
        shutil.rmtree(tmp_path, ignore_errors=True)

        shutil.copytree(dca2x2_path, tmp_path)
        dca2x2_path = os.path.join(dca2x2_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_dca_cmd, shell=True)
        good_result_old = np.loadtxt(os.path.join(tmp_path, "GoodOld/greenUp1.dat"))
        result = np.loadtxt("greenUp1.dat")
        os.chdir(base_path)
        # shutil.rmtree(tmp_path)

        len_result = result.shape[0]

        # test all
        np.testing.assert_allclose(
            result, good_result_old[:len_result], rtol=4e-3, atol=1e-3
        )

        # test the high frequencies
        np.testing.assert_allclose(
            result[70:75, 1::], good_result_old[70:75, 1::], rtol=1e-4, atol=1e-4
        )

        # the columns 3 and 5 should be zero (counting from zero)
        np.testing.assert_allclose(result[:, 3], result[:, 5], atol=1e-8)
        np.testing.assert_allclose(result[:, 3], np.zeros(len_result), atol=1e-9)

        # now test the observables
        with open(os.path.join(tmp_path, "Obs.json")) as fin:
            jj_result = json.load(fin)

        with open(os.path.join(tmp_path, "GoodOld/Obs.json")) as fin:
            jj_good = json.load(fin)

        keys = ["docc", "n", "sign", "k"]
        for key in keys:
            np.testing.assert_allclose(
                jj_result[key][0], jj_good[key][0], rtol=1e-3, atol=1e-3
            )

    def test_dca4x4(self):
        base_path = os.getcwd()
        dca4x4_path = os.path.join(base_path, "test/Simulations/4x4_DCA_b10_U3")
        tmp_path = os.path.join(base_path, "tmp/4x4_DCA_b10_U3")
        shutil.rmtree(tmp_path, ignore_errors=True)

        shutil.copytree(dca4x4_path, tmp_path)
        dca4x4_path = os.path.join(dca4x4_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_dca_cmd, shell=True)
        good_result_old = np.loadtxt(os.path.join(tmp_path, "Good/greenUp1.dat"))
        result = np.loadtxt("greenUp1.dat")
        os.chdir(base_path)
        # shutil.rmtree(tmp_path)

        len_result = result.shape[0]

        # test all
        np.testing.assert_allclose(
            result, good_result_old[:len_result], rtol=5e-3, atol=4e-3
        )

        # test the high frequencies
        np.testing.assert_allclose(
            result[130:140, 1::], good_result_old[130:140, 1::], rtol=1e-4, atol=1e-4
        )

        # now test the observables
        with open(os.path.join(tmp_path, "Obs.json")) as fin:
            jj_result = json.load(fin)

        with open(os.path.join(tmp_path, "Good/Obs.json")) as fin:
            jj_good = json.load(fin)

        keys = ["docc", "n", "sign", "k"]
        for key in keys:
            np.testing.assert_allclose(
                jj_result[key][0], jj_good[key][0], rtol=1e-3, atol=1e-3
            )

    def test_square2x2_tp_tpp(self):
        base_path = os.getcwd()
        square2x2_path = os.path.join(
            base_path, "test/Simulations/2x2Square_tp_tpp_beta50_U5"
        )
        tmp_path = os.path.join(base_path, "tmp/2x2Square_tp_tpp_beta50_U5")
        shutil.rmtree(tmp_path, ignore_errors=True)

        shutil.copytree(square2x2_path, tmp_path)
        square2x2_path = os.path.join(square2x2_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_cmd, shell=True)
        good_result = np.loadtxt(os.path.join(tmp_path, "StatsPat/green_moy.dat"))
        result = np.loadtxt("greenUp1.dat")
        os.chdir(base_path)

        # test the first few frequencies
        np.testing.assert_allclose(
            result[:100], good_result[:100], rtol=5e-3, atol=5e-3
        )

    def test_square2x2_2Orb(self):
        base_path = os.getcwd()
        square2x2_path = os.path.join(base_path, "test/Simulations/2x2Square_2Orb_U4.5")
        tmp_path = os.path.join(base_path, "tmp/2x2Square_2Orb_U4.5")
        shutil.rmtree(tmp_path, ignore_errors=True)

        shutil.copytree(square2x2_path, tmp_path)
        square2x2_path = os.path.join(square2x2_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_cmd, shell=True)
        good_result_old = np.loadtxt(
            os.path.join(tmp_path, "Stats_Good/greenUp112.dat")
        )
        result = np.loadtxt("greenUp1.dat")
        os.chdir(base_path)

        len_result = result.shape[0]

        # test all the green fct
        np.testing.assert_allclose(
            result, good_result_old[:len_result], rtol=5e-3, atol=5e-3
        )

        with open(os.path.join(tmp_path, "Obs.json")) as fin:
            jj_result = json.load(fin)

        with open(os.path.join(tmp_path, "Stats_Good/statsobs.json")) as fin:
            jj_good = json.load(fin)

        keys = ["docc", "n", "sign", "k"]
        for key in keys:
            np.testing.assert_allclose(
                jj_result[key][0], jj_good[key + ".dat"][0], rtol=1e-3, atol=1e-3
            )


if __name__ == "__main__":
    unittest.main()

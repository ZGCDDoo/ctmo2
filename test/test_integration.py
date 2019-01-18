import os
import shutil
import unittest
import numpy as np
import subprocess

# to run only one test: python -m unittest -v test.test_integration.TestIntegration.test_dca2x2


class TestIntegration(unittest.TestCase):

    mpi_cmd = "mpirun -np 4 ctmo params 1"
    mpi_dca_cmd = "mpirun -np 4 ctmo_DCA params 1"

    @classmethod
    def setUpClass(cls):
        print("\nIn TestIntegration.\n")
        shutil.rmtree("./tmp", ignore_errors=True)

    def test_holstein_dmft(self):
        base_path = os.getcwd()
        holstein_path = os.path.join(base_path, "test/Simulations/holstein_DMFT_T0.3")
        tmp_path = os.path.join(base_path, "tmp/holstein_DMFT_T0.3")
        shutil.copytree(holstein_path, tmp_path)
        holstein_path = os.path.join(holstein_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_cmd, shell=True)
        good_result = np.loadtxt(
            os.path.join(tmp_path, "GoodReza/Region1_LocalGFn.dat")
        )
        result = np.loadtxt("greenUp.dat")
        os.chdir(base_path)
        # shutil.rmtree(tmp_path)

        len_result = result.shape[0]
        len_high = len_result

        # test the low frequencies
        np.testing.assert_allclose(
            result[:len_high],
            good_result[:len_high],
            rtol=1e-3,
            atol=1e-3,
            verbose=True,
        )

        # Test the first few and last frequencies - imaginary part
        np.testing.assert_allclose(
            result[0:5, 2], good_result[0:5, 2], rtol=1e-3, atol=1e-3
        )
        np.testing.assert_allclose(
            result[200:205, 2], good_result[200:205, 2], rtol=1e-3, atol=1e-3
        )

        # test the high frequencies which should be more precise
        # np.testing.assert_allclose(
        #     result[len_result - 20 : len_result],
        #     good_result[len_result - 20 :],
        #     rtol=1e-2,
        #     atol=1e-2,
        # )

    def test_triangle2x2(self):
        base_path = os.getcwd()
        triangle_path = os.path.join(
            base_path, "test/Simulations/2x2Triangle_tp04_b5_U3"
        )
        tmp_path = os.path.join(base_path, "tmp/2x2Triangle_tp04_b5_U3")
        shutil.copytree(triangle_path, tmp_path)
        triangle_path = os.path.join(triangle_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_cmd, shell=True)
        good_result = np.loadtxt(os.path.join(tmp_path, "GoodNew/greenUp1.dat"))
        good_result_old = np.loadtxt(os.path.join(tmp_path, "GoodOld/greenUp1.dat"))
        result = np.loadtxt("greenUp.dat")
        os.chdir(base_path)
        # shutil.rmtree(tmp_path)

        np.testing.assert_allclose(result, good_result, rtol=5e-4, atol=5e-4)

        len_result = result.shape[0]

        # only test the first green function, because order is not the same.
        np.testing.assert_allclose(
            result[:, 0:3], good_result_old[:len_result, 0:3], rtol=1e-3, atol=1e-3
        )

    def test_dca2x2(self):
        base_path = os.getcwd()
        triangle_path = os.path.join(base_path, "test/Simulations/2x2_DCA_b5_U3")
        tmp_path = os.path.join(base_path, "tmp/2x2_DCA_b5_U3")
        shutil.copytree(triangle_path, tmp_path)
        triangle_path = os.path.join(triangle_path, tmp_path)
        os.chdir(tmp_path)

        subprocess.run(self.mpi_dca_cmd, shell=True)
        good_result_old = np.loadtxt(os.path.join(tmp_path, "GoodOld/greenUp1.dat"))
        result = np.loadtxt("greenUp.dat")
        os.chdir(base_path)
        # shutil.rmtree(tmp_path)

        len_result = result.shape[0]

        np.testing.assert_allclose(
            result, good_result_old[:len_result], rtol=5e-3, atol=5e-3
        )

        # np.testing.assert_allclose(
        #     result[len_result - 20 :],
        #     good_result_old[len_result - 20 : len_result],
        #     rtol=1e-2,
        #     atol=1e-2,
        # )

        # the columns 3 and 5 should be zero (counting from zero)
        np.testing.assert_allclose(result[:, 3], result[:, 5], atol=1e-4)
        np.testing.assert_allclose(result[:, 3], np.zeros(len_result), atol=1e-4)


if __name__ == "__main__":
    unittest.main()

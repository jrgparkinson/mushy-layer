import unittest
import mushyLayerRunUtils as util
import mock
import os

# Run with e.g.
# coverage run -m unittest test_utils.py; coverage html
# python -m unittest
class TestUtil(unittest.TestCase):

    def test_get_executable(self):

        with mock.patch("socket.gethostname", return_value="fake_host"):
            self.assertEqual(util.get_executable(), 'mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex')
            self.assertEqual(util.get_executable(dim=3), 'mushyLayer3d.Linux.64.mpiCC.gfortran.OPT.MPI.ex')
            self.assertEqual(util.get_executable(base_name='compare', dim=3), 'compare3d.Linux.64.mpiCC.gfortran.OPT.MPI.ex')

        with mock.patch("socket.gethostname", return_value="gyre3"):
            self.assertEqual(util.get_executable(), 'mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.GYRE.ex')


    def test_get_base_output_dir(self):

        with mock.patch("socket.gethostname", return_value="atmlxmaster"):
            self.assertEqual(util.get_base_output_dir(), '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/TestDevelopment-19Dec/')

        with mock.patch("socket.gethostname", return_value="atmlxlap005"):
            self.assertEqual(util.get_base_output_dir(), '/home/parkinsonjl/mushy-layer/test/output/')

        with mock.patch("socket.gethostname", return_value="somethingelse"):


            with mock.patch.dict('os.environ', {'MUSHY_LAYER_TEST_PATH': 'testing'}):
                self.assertEqual(util.get_base_output_dir(), 'testing')

            # Check that with no mushy layer path available we raise an error
            with mock.patch.dict('os.environ', {}):
                with self.assertRaises(SystemExit) as cm:
                    util.get_base_output_dir()()

                self.assertEqual(cm.exception.code, -1)

    def test_get_data_dir(self):

        with mock.patch("socket.gethostname", return_value="atmlxlap005"):
            self.assertEqual(util.get_data_dir(), '/home/parkinsonjl/mnt/sharedStorage/')

        with mock.patch("socket.gethostname", return_value="somethingelse"):
            with self.assertRaises(SystemExit) as cm:
                util.get_data_dir()

            self.assertEqual(cm.exception.code, -1)


    def test_get_matlab_base_command(self):

        with mock.patch("os.path.abspath", return_value='/path/to/parent/'):
            self.assertEqual(util.get_matlab_base_command(), 'cd /path/to/parent/matlab/MethodsPaper;\n'
                                                         'source /etc/profile.d/modules.sh\n'
                                                         'module load idl/870\n'
                                                         'module load matlab/R2018b'
                                                         '\n\n'
                                                         'matlab -nodisplay -nosplash -nodesktop -r' )

    def test_get_current_vcs_revision(self):
        with mock.patch("subprocess.check_output", return_value='version'.encode()):
            self.assertEqual(util.get_current_vcs_revision(), 'version')

    def test_get_mush_layer_dir(self):
        with mock.patch("os.path.abspath", return_value='/path/to/this/file/'):
            self.assertEqual(util.get_mushy_layer_dir(), '/path/to/this')

    def test_get_executable_name(self):
        pass

    def test_construct_run_name(self):
        pass


    def test_read_inputs(self):
        pass

    def test_write_inputs(self):
        pass

    def test_has_reached_steady_state(self):
        pass

    def test_time_since_folder_updated(self):
        pass

    def test_time_since_file_updated(self):
        pass

    def test_get_restart_file(self):
        pass

    def test_get_final_plot_file(self):
        pass

    def test_get_final_chk_file(self):
        pass


    def test_check_exec_exists(self):
        pass









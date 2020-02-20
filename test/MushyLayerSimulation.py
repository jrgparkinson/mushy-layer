import os
from mushyLayerRunUtils import write_inputs


class SimulationStatus:
    NONEXISTANT = 'Non-existant'
    UNSTARTED = 'Unstarted'
    RUNNING = 'Running'
    FINISHED = 'Finished'
    CRASHED = 'Crashed'


# This class describes a single, self contained, simulation
# which exists in it's own folder.
# It should make it easier to work with simulations from python:
#  creating them, running them, and analysing them
class MushyLayerSimulation:
    folder = ''
    status = SimulationStatus.NONEXISTANT
    num_proc = 0
    inputs_file_loc = ''
    inputs = None

    def __init__(self, folder):
        self.folder = folder
        self.inputs_file_loc = os.path.join(self.folder, 'inputs')

        self.determine_status()

    def run(self):
        if self.status == SimulationStatus.NONEXISTANT:
            self.make_simulation()

        if not self.status == SimulationStatus.UNSTARTED:
            print('Unable to run simulation in folder %s' % self.folder)
            print('Current status: %s' % self.status)
            return

    def make_simulation(self):
        if os.path.exists(self.folder):
            print('Unable to make simulation as folder exists: %s' % self.folder)
            return

        if self.inputs is None:
            print('Unable to make simulation as no inputs have been defined')
            return

        # Make folder
        os.makedirs(self.folder)

        # Write inputs file
        write_inputs(self.inputs_file_loc, self.inputs)

        self.status = SimulationStatus.UNSTARTED

    def determine_status(self):
        """ Work out what state the simulation is in """

        if not self.folder:
            self.status = SimulationStatus.NONEXISTANT
            return

        if not os.path.exists(self.folder):
            self.status = SimulationStatus.NONEXISTANT
            return

        output_files = os.listdir(self.folder)

        time_table_files = [f for f in output_files if 'time.table.' in f]
        pout_files = [f for f in output_files if 'pout.' in f]

        if len(pout_files) == 0:
            self.status = SimulationStatus.UNSTARTED
            return

        if len(time_table_files) == 0:
            self.status = SimulationStatus.RUNNING
            return

        # TODO: determine if simulation finished properly, or crashed

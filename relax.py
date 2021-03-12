import logging
import time

import pyrosetta as prs
from pyrosetta.rosetta.protocols.relax import ClassicRelaxCreator

import benchmark_prot_fetcher

prs.init(
    options="-ex1 -ex2", set_logging_handler=True, extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta"
)  # no output from the design process
prs.logging_support.set_logging_sink()
logger = logging.getLogger("rosetta")
logger.setLevel(logging.ERROR)


def relax_with_config(fa_reps):
    script_manager = prs.rosetta.protocols.relax.RelaxScriptManager()

    # get benchmark protein
    pdbs = benchmark_prot_fetcher.get_crystal()

    # with open('relax_script', 'wb') as f:
    #     f.write("repeat 5 \n\
    #             coord_cst_weight 1.0\n\
    #             scale:fa_rep {}\n\
    #             repack\n\
    #             scale:fa_rep {}\n\
    #             min 0.01\n\
    #             coord_cst_weight 0.5\n\
    #             scale:fa_rep {}\n\
    #             repack\n\
    #             scale:fa_rep {}\n\
    #             min 0.01\n\
    #             coord_cst_weight 0.0\n\
    #             scale:fa_rep {}\n\
    #             repack\n\
    #             scale:fa_rep {}\n\
    #             min 0.01\n\
    #             coord_cst_weight 0.0\n\
    #             scale:fa_rep {}\n\
    #             repack\n\
    #             min 0.00001\n\
    #             accept_to_best\n\
    #             endrepeat".format(*fa_reps))

        def minimize(self):
            ori_pose = pyrosetta.Pose()
            ori_pose.assign(self._pose)
            mut_pose = pyrosetta.Pose()
            mut_pose.assign(self._pose)
            #observer = pyrosetta.rosetta.protocols.moves.AddPyMOLObserver(mut_pose, True)
            score_before = self._score_fxn(ori_pose)
            # print(self._score_fxn(mut_pose))
            # Minimizing
            minimizer = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            movemap = pyrosetta.rosetta.core.kinematics.MoveMap()
            # movemap.set_bb(True)
            minimizer.movemap(movemap)
            #         minimizer.score_function(centroid_scorefxn)
            minimizer.score_function(self._score_fxn_weighted)
            minimizer.apply(mut_pose)
            score_after = self._score_fxn(mut_pose)
            # ------------- No Metropolis ------------------------
            pose_new = pyrosetta.Pose()
            pose_new.assign(mut_pose)
            energy_diff = -(score_after - score_before)
            # ----------------------------------------------------
            # print(energy_diff)

            return pose_new, energy_diff

            def relaxation(self):
                ori_pose = pyrosetta.Pose()
                ori_pose.assign(self._pose)
                mut_pose = pyrosetta.Pose()
                mut_pose.assign(self._pose)
                score_before = self._score_fxn(ori_pose)
                taskf = pyrosetta.rosetta.core.pack.task.TaskFactory()
                taskf.push_back(
                    pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
                taskf.push_back(
                    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
                packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
                    self._score_fxn_weighted)
                #         rlx.set_scorefxn(centroid_scorefxn)
                packer.task_factory(taskf)
                #         print(taskf.create_task_and_apply_taskoperations(mut_pose)) #print out the applied taks - useful for debug
                packer.apply(mut_pose)
                #         fa_switch = pyrosetta.rosetta.protocols.simple_moves.SwitchResidueTypeSetMover("fa_standard")
                #         fa_switch.apply(mut_pose)
                # print(self._score_fxn(mut_pose))
                score_after = self._score_fxn(mut_pose)

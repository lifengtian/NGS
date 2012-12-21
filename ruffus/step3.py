from ruffus import *
import sys

#---------------------------------------------------------------
#
#   first task
#
task1_param = [
                    [ None, 'job1.stage1'], # 1st job
                    [ None, 'job2.stage1'], # 2nd job
              ]

@files(task1_param)
def first_task(no_input_file, output_file):
    open(output_file, "w")
    #
    # pretend we have worked hard


#---------------------------------------------------------------
#
#   second task
#
task2_param = [
                    [ 'job1.stage1', "job1.stage2", "    1st_job"], # 1st job
                    [ 'job2.stage1', "job2.stage2", "    2nd_job"], # 2nd job
              ]

@follows(first_task)
@files(task2_param)
def second_task(input_file, output_file, extra_parameter):
    open(output_file, "w")
    print extra_parameter

pipeline_printout(sys.stdout, [second_task], verbose = 3)

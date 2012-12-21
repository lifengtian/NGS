from ruffus import *
import time

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

#---------------------------------------------------------------
#
#       Show flow chart and tasks before running the pipeline
#
print "Show flow chart and tasks before running the pipeline"
pipeline_printout_graph ( open("flowchart_before.png", "w"),
                         "png",
                         [second_task],
                         no_key_legend=True)

#---------------------------------------------------------------
#
#       Run
#
pipeline_run([second_task])


# modify job1.stage1
open("job1.stage1", "w").close()


#---------------------------------------------------------------
#
#       Show flow chart and tasks after running the pipeline
#
print "Show flow chart and tasks after running the pipeline"
pipeline_printout_graph ( open("flowchart_after.png", "w"),
                         "png",
                         [second_task],
                         no_key_legend=True)

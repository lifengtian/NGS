from ruffus import *

def first_task():
	print "hello"

@follows (first_task)
def second_task():
	print "world"

pipeline_run([second_task])


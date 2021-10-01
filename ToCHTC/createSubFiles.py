import argparse
import os, sys
import datetime

time = datetime.datetime.now()

subFile = 'block.sub'
dockerImage = 'nathanlazarus/oilgame:v2'

script = 'OilGame.py'
helper_functions = 'helper_functions.py'

WhichStage = 'WhichStage.txt'
Donefile = 'Done.txt'

# parser = argparse.ArgumentParser()
# parser.add_argument('Nplayers', type = int, help = 'number of players')
# parser.add_argument('Nblocks', type = int, help = 'number of blocks')
# parser.add_argument('entire_problem_max_DS', type = int, help = 'size of blocks')
# parser.add_argument('maxWorkers', type = int, help = 'maximum number of worker jobs')
# args = parser.parse_args()

if not os.path.exists(WhichStage):
    thisStage = 1
else:
    with open(WhichStage, 'r') as f:
    	thisStage = int(next(f))


nextStage = thisStage + 1

with open(WhichStage, 'w+') as f:
	f.write(str(nextStage))

with open('iter' + str(thisStage) + 'time.txt', 'w') as f:
    f.write(str(time) + ', ' + str(thisStage))

# Nstates = entire_problem_max_DS ** Nplayers

# Niter = Nstates

def states_summing_to_x(length, total_sum, min_vals, max_vals):
    if length == 1:
        yield (total_sum,)
    else:
        range_min = max(min_vals[len(min_vals) - length], total_sum - sum(max_vals[len(max_vals) - length + 1:]))
        range_max = min(total_sum - sum(min_vals[len(min_vals) - length + 1:]), max_vals[len(max_vals) - length])
        for value in range(range_min, range_max + 1):
            for permutation in states_summing_to_x(length - 1, total_sum - value, min_vals, max_vals):
                yield (value,) + permutation


parser = argparse.ArgumentParser()
parser.add_argument('Nplayers', type = int, help = 'number of players')
parser.add_argument('Nblocks', type = int, help = 'number of blocks')
parser.add_argument('entire_problem_max_DS', type = int, help = 'size of blocks')
parser.add_argument('maxWorkers', type = int, help = 'maximum number of worker jobs')

args = parser.parse_args()

globals().update(args.__dict__)



Nblocks_by_dim = [Nblocks] * Nplayers
grand_scheme_min_DS_list = [1] * Nplayers


thisStageTot = list(range(sum(Nblocks_by_dim), sum(grand_scheme_min_DS_list) - 1, -1))[thisStage - 1]
thisStageJobs = states_summing_to_x(Nplayers, thisStageTot, grand_scheme_min_DS_list, Nblocks_by_dim)
Njobs = len(list(thisStageJobs))

isFinalStage = int(thisStageTot == sum(grand_scheme_min_DS_list))
print('isFinalStage', isFinalStage)

with open(Donefile, 'w+') as f:
    f.write(str(isFinalStage))

if isFinalStage == True:
    os.remove(WhichStage)

if thisStage > 1:
    lastStageTot = list(range(sum(Nblocks_by_dim), sum(grand_scheme_min_DS_list) - 1, -1))[thisStage - 2]
    lastStageJobs = states_summing_to_x(Nplayers, lastStageTot, grand_scheme_min_DS_list, Nblocks_by_dim)
    past_state_files = ['valuearray_' + str(Nplayers) + '_' + '_'.join(map(str, job)) + '.nc' for job in lastStageJobs]
else:
    past_state_files = []

input_files = ', '.join([script, helper_functions] + past_state_files)


NworkersThisBlock = min(Njobs, maxWorkers)

sub_text = """\
# {0}

workerID                = $(Process) + 1
universe                = docker
docker_image            = {1}
executable              = python
transfer_executable     = False
arguments               = {2} {3} {4} {5} {6} --maxworkers {7} -i $INT(workerID)
transfer_input_files    = {8}
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
error                   = outfiles/worker{6}_$INT(workerID).err
output                  = outfiles/worker{6}_$INT(workerID).out
log                     = outfiles/worker{6}.log
max_retries             = 10
periodic_release        = (HoldReasonCode == 35) && (NumJobStarts < 10)

Rank                    = kflops
request_disk            = 300MB
request_memory          = 300MB
queue {7}
""".format(subFile, dockerImage, script, Nplayers, Nblocks, entire_problem_max_DS, thisStage, NworkersThisBlock, input_files)


with open(subFile, 'w') as sub_file:
    sub_file.write(sub_text)

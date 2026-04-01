#!/usr/bin/env python
##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Module to Analyse memory usage information from PBS and /usr/bin/time for
lfric_apps jobs run on the Met Office EX systems.

"""

import argparse
import datetime
import json
import os
import re
import statistics
import time

import numpy as np
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


def num_fmt(y, pos):
    """Helper function to allow y axis formatting of MiB to GiB."""
    label = f'{y/1024:,.0f}'
    return label


yfmt = tkr.FuncFormatter(num_fmt)

KB_TO_MIB = 1000/(1024**2)

run_env_keys = ["TARGET_PLATFORM",
                "RUN_METHOD",
                "HYPERTHREADS",
                "CORES_PER_NODE",
                "NUMA_REGIONS_PER_NODE",
                "OMP_NUM_THREADS",
                "TOTAL_RANKS",
                "XIOS_SERVER_MODE",
                "XIOS_SERVER_RANKS",
                "EXEC_NAME",
                "PAT_EXE_EXTEN"
                ]


def plot_run_job(run, out_filename):
    """
    Parse information from job, job.out & job.err.
    Produce a plot of the Memory profile from the information.
    `run` should be a path to a job folder, where job.err, job.out & job
    are available. It is expected that the selected job was run with
    MEMORY_PROFILE=1 on Met Office EX host machines for the collection &
    labelling of data from the executed process.

    """

    run_stats = {}
    # Parse run configuration information from the job & job.out PBS files.
    if not os.path.exists(run):
        raise FileNotFoundError(f'{run} does not exist')
    job_name = os.path.basename(os.path.dirname(run))
    with open(os.path.join(run, 'job'), encoding="utf-8") as jf:
        jfr = jf.read()
    # assume 1 node, unless PBS configures nodes
    nnodes = 1
    # re pattern to fetch number of nodes from PBS config
    nnodes_rep = re.compile('#PBS -l select=([0-9]+)')
    if len(nnodes_rep.findall(jfr)) == 1:
        nnodes, = nnodes_rep.findall(jfr)
    nnodes = int(nnodes)
    # re pattern to fetch number of nodes from PBS config
    cores_per_node_mpiproc_rep = re.compile(r'[\s\S]*#PBS -l select='
                                            r'([0-9]+):coretype=milan:'
                                            r'mpiprocs=([0-9]+).*[\s\S]*')
    matcher = cores_per_node_mpiproc_rep.match(jfr)
    mpiprocs = None
    if matcher:
        mpiprocs = int(matcher.group(2))

    if not os.path.exists(os.path.join(run, 'job.out')):
        time.sleep(5)
        if not os.path.exists(os.path.join(run, 'job.out')):
            raise FileNotFoundError(f'{run}/job.out does not exist')

    with open(os.path.join(run, 'job.out'), encoding="utf-8") as jof:
        jofr = jof.read()
        # re pattern to find max memory usage report by node
        mem_per_node_rep = re.compile(r'[\s\S]*Maximum memory usage per node'
                                      r' \(MiB\)((\n  nid.*)+)\n[\s\S]*')
        matcher = mem_per_node_rep.match(jofr)
        if matcher:
            mem_per_node = matcher.group(1).split('\n')[1:]
            mem_per_node.sort(key=lambda x: int(x.split()[1]))
            mem_per_node = [{'node': m.split()[0], 'ncount': m.split()[1],
                             'memMiB': m.split()[2]} for m in mem_per_node]
        else:
            raise ValueError('Failed to parse memory per node from job.out. '
                             'Check that the PBS script has run and output '
                             'and that both MEMORY_PROFILE is true & '
                             'TARGET_PLATFORM = "meto-ex1a" in the '
                             'environment variables passed to launch-exe.')

    run_env_vars = {}
    # parse environment variables from job.out
    for k, v in re.findall((f'({"|".join([e for e in run_env_keys])})='
                            r'([0-9a-zA-Z\-\_]+)'), jofr):
        run_env_vars[k] = v
    mps = ""
    # try to parse wallclock time from job.out
    # "walltime": "00:01:16" or
    # Elapsed Time:          0:31:49 (1909 seconds, 80% of total)
    wclock_rep = re.compile('(?:"walltime": "|Elapsed Time:[ ]+)'
                            '([0-9]+:[0-9]+:[0-9]+)')
    wclock = "missing"
    if len(wclock_rep.findall(jofr)) == 1:
        wclock, = wclock_rep.findall(jofr)
    else:
        # but the job.out can be read before PBS has reported the wall clock
        # time, so wait, then try again
        time.sleep(11)
        with open(os.path.join(run, 'job.out'), encoding="utf-8") as jof:
            jofr = jof.read()
        # but this can also fail, as PBS End of Job may not get to cylc job.out
        # so only update wclock if it is found
        if wclock_rep.match(jofr) is not None:
            wclock, = wclock_rep.findall(jofr)
    if wclock == "missing":
        # wclock = "missing"
        # 2025-06-04T13:50:31Z INFO - started
        # 2025-06-04T14:48:33Z INFO - succeeded
        start_pattern = re.compile(r'([0-9\-:T]+)Z INFO - started')
        astart, = start_pattern.findall(jofr)
        success_pattern = re.compile(r'([0-9\-:T]+)Z INFO - succeeded')
        asuccess, = success_pattern.findall(jofr)
        wt = datetime.datetime.fromisoformat(asuccess) \
            - datetime.datetime.fromisoformat(astart)
        wclock = str(wt)

    if mpiprocs is not None:
        mps = f"{mpiprocs} mpiprocs, "
    runtitle = (f"{job_name}: wallclock= {wclock}\n"
                f"{run_env_vars.get('TARGET_PLATFORM', '')} "
                f"{nnodes} nodes, "
                f"{mps}"
                f"{run_env_vars.get('EXEC_NAME', '')} "
                f"{run_env_vars.get('TOTAL_RANKS', '')} ranks, "
                f"{run_env_vars.get('OMP_NUM_THREADS', '')} OMP_threads, "
                f"{run_env_vars.get('', '')}")
    if run_env_vars.get('XIOS_SERVER_MODE') == 'True':
        rranks = run_env_vars.get('XIOS_SERVER_RANKS', '')
        runtitle = runtitle + f"xios_server {rranks} ranks, "

    run_stats['run'] = {}
    run_stats['run']['wallclock'] = wclock
    run_stats['run']['nnodes'] = nnodes
    tranks = int(run_env_vars.get('TOTAL_RANKS', '-1'))
    run_stats['run'][f"{run_env_vars.get('EXEC_NAME', '')} ranks"] = tranks
    run_stats['run']['OMP_threads'] = run_env_vars.get('OMP_NUM_THREADS', '')
    if run_env_vars.get('XIOS_SERVER_MODE') == 'True':
        run_stats['run']['xios_server ranks'] = \
            int(run_env_vars.get('XIOS_SERVER_RANKS', ''))
    # Run configuration information is now encoded into the title to be used
    # for plotting.
    # Next, parse run configuration information from the job.err PBS file;
    # this was encoded into job.err by the /usr/bin/time utility which wrapped
    # the run command for mpiexec.
    with open(os.path.join(run, 'job.err'), encoding="utf-8") as jef:
        jefr = jef.read()
        # re pattern for max_mem_{process} for mem per rank
        mem_per_rank_rep = re.compile('(([a-z]+[0-9]+) ([0-9]+): '
                                      'max_mem_([a-zA-Z0-9_]+)_([0-9]+)kb)+')
        mem_per_rank = mem_per_rank_rep.findall(jefr)
        mem_per_rank.sort(key=lambda x: (x[1], int(x[2])))
        mem_per_rank = [{'node': m[1], 'rank': int(m[2]), 'exec': m[3],
                         'memkB': int(m[4]), 'log': m[0]}
                        for m in mem_per_rank]

    for k, v in re.findall((r'(MPICH_RANK_REORDER_METHOD)[ ]*='
                            r'[ ]*([0-9a-zA-Z\-\_]+)'), jefr):
        run_env_vars[k] = v

    # Now there is a multi-line title string, with information captured
    # from the run Environment.
    # There are also data objects: mem-per_node & mem_per_rank
    # which are now suitable to create a plot with.

    nodes = [m['node'] for m in mem_per_node]

    width = 0.45
    offset = 0.5

    # Begin plot configuration, using a matplotlib figure.
    fig, ax = plt.subplots(figsize=(24, 9))
    # container for stacking bar plots by y value
    bottom = np.zeros(len(nodes))

    # Define X axis.
    x = np.arange(len(nodes))

    # Create a bar chart of maximum memory reported to PBS per node.
    nbar = ax.bar(x+offset, [int(m['memMiB']) for m in mem_per_node],
                  tick_label=nodes, width=width, color='purple', label='node')

    # Iterate through each element of the memory per rank list.
    # Assign a label and colour, then plot this element as a stacked bar
    # updating the stacking container and adding labels.
    # handle minimal legend for multi-bar stack with these flags
    lfirst = True
    xfirst = True
    lfh = None
    xfh = None
    for m in mem_per_rank:
        label = "executable"
        col = "black"
        xval = nodes.index(m['node'])
        y = m['memkB']*KB_TO_MIB
        if m['exec'] == 'lfric' or m['exec'] == 'lfric_atm' or \
           m['exec'] == 'lfric2um' or m['exec'] == 'um2lfric':
            label = m['exec']
            col = 'green'
            if lfirst:
                lfirst = False
                lfh = ax.bar(xval, y, width=width, bottom=bottom[xval],
                             color=col, label=label, edgecolor='black')
            else:
                ax.bar(xval, y, width=width, bottom=bottom[xval], color=col,
                       label=label, edgecolor='black')
        elif m['exec'] == 'xios' or m['exec'] == 'xios_server':
            label = 'xios_server'
            col = 'blue'
            if xfirst:
                xfirst = False
                xfh = ax.bar(xval, y, width=width, bottom=bottom[xval],
                             color=col, label=label, edgecolor='black')
            else:
                ax.bar(xval, y, width=width, bottom=bottom[xval], color=col,
                       label=label, edgecolor='black')
        bottom[xval] += y

    # Rotate x axis labels for readability.
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    # Set y axis labels to be in GiB for readability.
    ax.yaxis.set_major_formatter(yfmt)
    ax.set_ylabel('Max Memory GiB')

    # Calculate summary statistics to append to the title.
    all_lf = [m['memkB']*KB_TO_MIB for m in mem_per_rank
              if m['exec'] == 'lfric_atm' or m['exec'] == 'lfric' or
              m['exec'] == 'lfric2um' or m['exec'] == 'um2lfric']
    all_x = [m['memkB']*KB_TO_MIB for m in mem_per_rank
             if m['exec'] == 'xios_server' or m['exec'] == 'xios']
    all_n = [int(int(m['memMiB'])/1024) for m in mem_per_node]

    tstats = (f'\n{int(min(all_lf))} ≤ lfric_atm (max mem MiB): '
              f'{int(statistics.mean(all_lf))}'
              f' ± σ={int(statistics.stdev(all_lf))} '
              f'≤ {int(max(all_lf))}')

    run_stats[run_env_vars.get('EXEC_NAME', '')] = {}
    run_stats[run_env_vars.get('EXEC_NAME', '')]['mean'] = \
        statistics.mean(all_lf)
    run_stats[run_env_vars.get('EXEC_NAME', '')]['max'] = max(all_lf)
    run_stats[run_env_vars.get('EXEC_NAME', '')]['min'] = min(all_lf)
    run_stats[run_env_vars.get('EXEC_NAME', '')]['std_dev'] = \
        statistics.stdev(all_lf)
    run_stats[run_env_vars.get('EXEC_NAME', '')]['units'] = 'max mem MiB'
    if all_x:
        tstats = tstats + (f'\n{int(min(all_x))} ≤ xios_server (max mem MiB): '
                           f'{int(statistics.mean(all_x))} ± '
                           f'σ={int(statistics.stdev(all_x))}'
                           f' ≤ {int(max(all_x))}')
        ax.legend(handles=[nbar, lfh, xfh])
        run_stats['xios_server'] = {}
        run_stats['xios_server']['mean'] = statistics.mean(all_x)
        run_stats['xios_server']['max'] = max(all_x)
        run_stats['xios_server']['min'] = min(all_x)
        run_stats['xios_server']['std_dev'] = statistics.stdev(all_x)
        run_stats['xios_server']['units'] = 'max mem MiB'
    else:
        ax.legend(handles=[nbar, lfh])
    if len(all_n) > 1:
        tstats = tstats + (f'\n{int(min(all_n))} ≤ node (max mem GiB): '
                           f'{int(statistics.mean(all_n))} ± '
                           f'σ={int(statistics.stdev(all_n))}'
                           f' ≤ {int(max(all_n))}')
    elif len(all_n) == 1:
        tstats = tstats + (f'\nSingle node (max mem GiB): '
                           f'{int(statistics.mean(all_n))}')

    run_stats['node'] = {}
    run_stats['node']['mean'] = statistics.mean(all_n)
    run_stats['node']['max'] = max(all_n)
    run_stats['node']['min'] = min(all_n)
    if len(all_n) > 1:
        run_stats['node']['std_dev'] = statistics.stdev(all_n)
    run_stats['node']['units'] = 'max mem MiB'

    # Add a horizontal line for Milan max memory
    ax.axhline(y=237*1024, color='r', linestyle='--')
    ax.text(0.5, 240*1024, 'EX Milan Estimated Usable Memory Threshold',
            color='black', va='center', ha='left',
            transform=ax.get_yaxis_transform())
    plt.title(runtitle + tstats)
    fig.savefig(out_filename + '.png')
    with open(out_filename + '.json', 'w', encoding="utf-8") as jsonout:
        jsonout.write(json.dumps(run_stats))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_job")

    args = parser.parse_args()
    this_run = args.run_job
    plot_filename = os.path.join(os.environ['CYLC_TASK_LOG_DIR'],
                                 "memory_profile")
    plot_run_job(this_run, plot_filename)

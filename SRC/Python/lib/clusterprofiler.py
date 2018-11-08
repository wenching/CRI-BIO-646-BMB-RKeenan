"""clusterProfiler"""
import os
import sys
import logging
import datetime

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def clusterprofiler(args, sw_cfg, task_cfg):
    """
    clusterProfiler

    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {clusterprofiler}
    :parm task_cfg: a dictionary of corresponding task configureation {case_vs_ctrl}
    :returns: returns corresponding code snippets, in which will be written to shell_script_path if provided
    :raises keyError: NA
    """
    # logging
    if args.log_file is None:
        log_file_path = '{}.{}.log'.format(
            SELF_FILE_PATH,
            datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    else:
        log_file_path = args.log_file

    formatter = "%(asctime)-15s %(levelname)-8s %(message)s"
    logging.basicConfig(
        level=[logging.NOTSET, logging.DEBUG, logging.INFO,
               logging.WARNING, logging.ERROR, logging.CRITICAL][1],
        format=formatter,
        filename=log_file_path)

    logging.info("[ {} ] Make Shell Script\n".format(SELF_FILE))

    sw_cfg_tool = sw_cfg['clusterprofiler']

    task_cfg_proj = task_cfg['clusterprofiler']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_proj else '''&>> {}'''.format(task_cfg_proj['log_file_path'])

    target_genome = next(iter(task_cfg['references']))
    task_cfg_reference = task_cfg['references'][target_genome]

    shell_script_path = '' if 'shell_script_path' not in task_cfg_proj else task_cfg_proj['shell_script_path']
    
    cmd_list = []

    if sw_cfg_tool['module'] is not None:
        cmd_list.append(sw_cfg_tool['module'])

    cmd_list.append('''{} \\
{} \\
--inPath="{}" \\
--compPair="{}" \\
--outPath="{}" \\
--assembly="{}" \\
--logPath="{}" \\
{}
\n'''.format(
         sw_cfg_tool['exe'],
         sw_cfg_tool['script'],
         ','.join(task_cfg_proj['in_file_path_dict'].values()),
         ','.join(task_cfg_proj['in_file_path_dict'].keys()),
         task_cfg_proj['out_file_path_list'][0],
         target_genome,
         '' if 'log_file_path' not in task_cfg_proj else task_cfg_proj['log_file_path'],
         log_file_path_cmd
     ))

    if shell_script_path:
        try:
            with open(shell_script_path, "w") as outfile:
                outfile.write(
                    '''\n{}\n'''.format(
                        '\n\n'.join(cmd_list)
                    )
                )
                
                outfile.write("\nexit 0\n")
                outfile.close()
        except IOError as exc:
            print(exc)

    logging.debug("[ {} ] Make Shell Script - DONE\n".format(SELF_FILE))

    return cmd_list

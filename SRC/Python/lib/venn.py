"""Venn Diagram"""
import os
import sys
import logging
import datetime

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def venn(args, sw_cfg, task_cfg):
    """
    Venn Diagram

    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {venn}
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

    sw_cfg_venn = sw_cfg['venn']

    task_cfg_comp = task_cfg[next(iter(task_cfg))]
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_comp else '''&>> {}'''.format(task_cfg_comp['log_file_path'])

    shell_script_path = '' if 'shell_script_path' not in task_cfg_comp else task_cfg_comp['shell_script_path']
    
    cmd_list = []

    if sw_cfg_venn['module'] is not None:
        cmd_list.append(sw_cfg_venn['module'])

    cmd_list.append('''{} \\
{} \\
--inPath="{}" \\
--outPath="{}" \\
--anchorPath="{}" \\
--logPath="{}" \\
{}
\n'''.format(
         sw_cfg_venn['exe'],
         sw_cfg_venn['script'],
         ','.join(task_cfg_comp['in_file_path_list']),
         task_cfg_comp['out_file_path_list'][0],
         task_cfg_comp['anchor_file_path'],
         '' if 'log_file_path' not in task_cfg_comp else task_cfg_comp['log_file_path'],
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

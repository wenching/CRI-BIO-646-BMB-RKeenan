"""sambamba"""
import os
import sys
import logging
import datetime
import re

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def index(args, sw_cfg, task_cfg):
    """
    sambamba::index
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {sambamba}
    :parm task_cfg: a dictionary of corresponding task configureation {in_file_path_list, .log_file_path, .shell_script_path}
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

    logging.info("[ {}::index ] Make Shell Command\n".format(SELF_FILE))

    sw_cfg_sambamba = sw_cfg['sambamba']
    nthreads = 1 if 'nthreads' not in sw_cfg_sambamba else min(1, sw_cfg_sambamba['threads'])
  
    task_cfg_sambamba = task_cfg['sambamba']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_sambamba else '''&>> {}'''.format(task_cfg_sambamba['log_file_path'])
    
    shell_script_path = '' if "shell_script_path" not in task_cfg_sambamba else task_cfg_sambamba['shell_script_path']
    
    cmd_list = []

    for i,v in enumerate(task_cfg_sambamba['in_file_path_list']):
        cmd_list.append('''{} \\
index \\
--nthreads={} \\
{} \\
{} \\
{}'''.format(
             sw_cfg_sambamba['exe'],
             nthreads,
             v,
             re.sub(r""".bam$""", ".bai", v),
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


    logging.debug("[ {}::index ] Make Shell Command - DONE\n".format(SELF_FILE))

    return cmd_list


def merge(args, sw_cfg, task_cfg):
    """
    sambamba::merge
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {sambamba}
    :parm task_cfg: a dictionary of corresponding task configureation {in_file_path_list, .log_file_path, .shell_script_path, out_file_path_list}
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

    logging.info("[ {}::merge ] Make Shell Command\n".format(SELF_FILE))

    sw_cfg_sambamba = sw_cfg['sambamba']
    nthreads = 1 if not 'nthreads' in sw_cfg_sambamba else sw_cfg_sambamba['nthreads']
    
    task_cfg_sambamba = task_cfg['sambamba']
    in_file_path_list = task_cfg_sambamba['in_file_path_list']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_sambamba else '''&>> {}'''.format(task_cfg_sambamba['log_file_path'])
    
    shell_script_path = "" if "shell_script_path" not in task_cfg_sambamba else task_cfg_sambamba['shell_script_path']

    cmd_list = []

    if len(in_file_path_list) == 1:
        cmd_list.append('''ln -sf \\
{} \\
{} \\
{}'''.format(
         in_file_path_list[0],
         task_cfg_sambamba['out_file_path_list'][0],
         log_file_path_cmd         
     ))
    else:
        if sw_cfg_sambamba['module'] is not None:
            cmd_list.append(sw_cfg_sambamba['module'])

        for i,v in enumerate(in_file_path_list):
            cmd_list.append('''{} \\
merge \\
--nthreads={} \\
{} \\
{} \\
{}'''.format(
         sw_cfg_sambamba['exe'],
         nthreads,
         task_cfg_sambamba['out_file_path_list'][0],
         ' '.join(task_cfg_sambamba['in_file_path_list']),
         log_file_path_cmd         
     ))
    
    task_cfg_sambamba_index = \
        {'sambamba':
            {
                'in_file_path_list': task_cfg_sambamba['out_file_path_list'],
                'log_file_path': None if 'log_file_path' not in task_cfg_sambamba else task_cfg_sambamba['log_file_path']
            }            
        }

    cmd_list.extend(
        index(
            args,
            sw_cfg,
            task_cfg_sambamba_index
        )
    )    
    
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


    logging.debug("[ {}::merge ] Make Shell Command - DONE\n".format(SELF_FILE))

    return cmd_list
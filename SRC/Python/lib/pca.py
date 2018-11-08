"""PCA"""
import os
import sys
import logging
import datetime

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def pca(args, sw_cfg, task_cfg):
    """
    PCA

    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {pca}
    :parm task_cfg: a dictionary of corresponding task configureation {pca, library_info}
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

    sw_cfg_pca = sw_cfg['pca']

    task_cfg_pca = task_cfg['pca']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_pca else '''&>> {}'''.format(task_cfg_pca['log_file_path'])

    shell_script_path = '' if 'shell_script_path' not in task_cfg_pca else task_cfg_pca['shell_script_path']
    
    library_2_group_key = 'NULL' if 'library_info' not in task_cfg_pca else ','.join(task_cfg_pca['library_info'].keys())
    library_2_group_value = 'NULL' if 'library_info' not in task_cfg_pca else ','.join(task_cfg_pca['library_info'].values())

    cmd_list = []

    if sw_cfg_pca['module'] is not None:
        cmd_list.append(sw_cfg_pca['module'])

    cmd_list.append('''{} \\
{} \\
--inPath="{}" \\
--outPath="{}" \\
--ids="{}" \\
--groups="{}" \\
--logPath="{}" \\
{}
\n'''.format(
         sw_cfg_pca['exe'],
         sw_cfg_pca['script'],
         task_cfg_pca['in_file_path_list'][0],
         task_cfg_pca['out_file_path_list'][0],
         library_2_group_key,
         library_2_group_value,
         '' if 'log_file_path' not in task_cfg_pca else task_cfg_pca['log_file_path'],
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

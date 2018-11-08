"""FastQC."""
import os
import sys
import logging
import datetime

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def fastqc(args, sw_cfg, task_cfg):
    """
    FastQC
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {fastqc}
    :parm task_cfg: a dictionary of corresponding task configureation {in_file_path_list, out_dir_path, shell_script_path, log_file_path}
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

    sw_cfg_fastqc = sw_cfg['fastqc']
    
    task_cfg_fastqc = task_cfg['fastqc']
    shell_script_path = task_cfg_fastqc['shell_script_path']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_fastqc else '''&>> {}'''.format(task_cfg_fastqc['log_file_path'])

    cmd_list = []
    
    if sw_cfg_fastqc['module'] is not None:
        cmd_list.append(sw_cfg_fastqc['module'])
    
    cmd_list.append('''{} \\
--extract \\
-o {} \\
-t {} \\
--nogroup \\
{} \\
{}'''.format(
         sw_cfg_fastqc['exe'],
         task_cfg_fastqc['out_dir_path'],
         len(task_cfg_fastqc['in_file_path_list']),
         ','.join(task_cfg_fastqc['in_file_path_list']),
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

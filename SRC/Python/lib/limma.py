"""limma"""
import os
import sys
import logging
import datetime

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def limma(args, sw_cfg, task_cfg):
    """
    limma
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {limma}
    :parm task_cfg: a dictionary of corresponding task configureation {limma, references}
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

    sw_cfg_limma = sw_cfg['limma']
    
    task_cfg_limma = task_cfg['limma']
    shell_script_path = task_cfg_limma['shell_script_path']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_limma else '''&>> {}'''.format(task_cfg_limma['log_file_path'])

    task_cfg_reference = task_cfg['references'][next(iter(task_cfg['references']))]

    cmd_list = []
    
    if sw_cfg_limma['module'] is not None:
        cmd_list.append(sw_cfg_limma['module'])
    
    cmd_list.append('''{} \\
{} \\
--metadata="{}" \\
--GTF="{}" \\
--inPath="{}" \\
--outPath="{}" \\
--filterLowExpr="{}" \\
--criteria="{}" \\
--foldChange={} \\
--pval={} \\
--Case2CtrlPairs="{}" \\
--application="{}" \\
--logPath="{}" \\
{}
\n'''.format(
         sw_cfg_limma['exe'],
         sw_cfg_limma['script'],
         task_cfg_limma['metadata_path'],
         task_cfg_reference['anno_gtf'],
         task_cfg_limma['in_dir_path'],
         task_cfg_limma['out_file_path_list'][0],
         "TRUE" if task_cfg_limma['filter_4_low_expr'] else "FALSE",
         task_cfg_limma['criteria'],
         task_cfg_limma['fold_change'],
         task_cfg_limma['pval'],
         '~'.join(task_cfg_limma['comparison'].keys()),
         task_cfg_limma['application'],
         '' if 'log_file_path' not in task_cfg_limma else task_cfg_limma['log_file_path'],
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

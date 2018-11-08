"""DESeq2"""
import os
import sys
import logging
import datetime

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def deseq2(args, sw_cfg, task_cfg):
    """
    DESeq2
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {deseq2}
    :parm task_cfg: a dictionary of corresponding task configureation {deseq2, references}
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

    sw_cfg_deseq2 = sw_cfg['deseq2']
    
    task_cfg_deseq2 = task_cfg['deseq2']
    shell_script_path = task_cfg_deseq2['shell_script_path']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_deseq2 else '''&>> {}'''.format(task_cfg_deseq2['log_file_path'])

    task_cfg_reference = task_cfg['references'][next(iter(task_cfg['references']))]

    cmd_list = []
    
    if sw_cfg_deseq2['module'] is not None:
        cmd_list.append(sw_cfg_deseq2['module'])
    
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
         sw_cfg_deseq2['exe'],
         sw_cfg_deseq2['script'],
         task_cfg_deseq2['metadata_path'],
         task_cfg_reference['anno_gtf'],
         task_cfg_deseq2['in_dir_path'],
         task_cfg_deseq2['out_file_path_list'][0],
         "TRUE" if task_cfg_deseq2['filter_4_low_expr'] else "FALSE",
         task_cfg_deseq2['criteria'],
         task_cfg_deseq2['fold_change'],
         task_cfg_deseq2['pval'],
         '~'.join(task_cfg_deseq2['comparison'].keys()),
         task_cfg_deseq2['application'],
         '' if 'log_file_path' not in task_cfg_deseq2 else task_cfg_deseq2['log_file_path'],
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

"""Subread"""
import os
import sys
import logging
import datetime

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def featurecounts(args, sw_cfg, task_cfg):
    """
    Subread::featureCounts

    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {featurecounts}
    :parm task_cfg: a dictionary of corresponding task configureation {featurecounts, references}
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

    sw_cfg_featurecounts = sw_cfg['featurecounts']

    task_cfg_main = None if 'main' not in task_cfg.keys() else task_cfg['main']

    task_cfg_featurecounts = task_cfg['featurecounts']
    out_file_path = task_cfg_featurecounts['out_file_path_list'][0]
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_featurecounts else '''&>> {}'''.format(task_cfg_featurecounts['log_file_path'])
    
    task_cfg_ref = task_cfg['references']

    shell_script_path = '' if 'shell_script_path' not in task_cfg_featurecounts else task_cfg_featurecounts['shell_script_path']


    args_list = []
    
    annot_ext = sw_cfg_featurecounts['annot_ext']
    isGTFAnnotationFile = sw_cfg_featurecounts['isGTFAnnotationFile']
    annot_inbuilt = sw_cfg_featurecounts['annot_inbuilt']
    
    if annot_ext is not None:
        if annot_ext not in ['GTF', 'SAF']:
            logging.error("[ {} ] not supported annotation extension: [ {} ]\n".format(SELF_FILE, annot_ext))

        args_list.append('''-F {}'''.format(annot_ext))
        if annot_ext == "GTF":
            args_list.append('''-a {}'''.format(task_cfg_ref[next(iter(task_cfg_ref))]['anno_gtf']))
        else:
            args_list.append('''-a {}'''.format(sw_cfg_featurecounts['annot_file']))  
    else:
        if annot_inbuilt in ['hg19', 'hg38', 'mm9', 'mm10']:
            args_list.append('''-a {}'''.format(annot_inbuilt))
        else:
            logging.error("[ {} ] not supported inbuilt annotation: [ {} ]\n".format(SELF_FILE, annot_inbuilt))

    args_list.append('''-t {}'''.format(sw_cfg_featurecounts['GTF_featureType']))
    args_list.append('''-g {}'''.format(sw_cfg_featurecounts['GTF_attrType']))

    if sw_cfg_featurecounts['primaryOnly']: 
        args_list.append('--primary')
    else:
        if sw_cfg_featurecounts['countMultiMappingReads']: args_list.append('-M')
        if sw_cfg_featurecounts['fraction']: args_list.append('--fraction')
    
    if task_cfg_featurecounts['strandSpecific'] is not None: args_list.append('''-s {}'''.format(task_cfg_featurecounts['strandSpecific']))
    if sw_cfg_featurecounts['juncCounts']: args_list.append('-J')
    if task_cfg_featurecounts['isPairedEnd']:
        args_list.append('-p')
        
        if sw_cfg_featurecounts['requireBothEndsMapped']: args_list.append('-B')
        if sw_cfg_featurecounts['countChimericFragments']: args_list.append('-C')


    cmd_list = []

    if sw_cfg_featurecounts['module'] is not None:
        cmd_list.append(sw_cfg_featurecounts['module'])

    cmd_list.append('''{} \\
{} \\
-o {} \\
{} \\
{}'''.format(
             sw_cfg_featurecounts['exe'],
             ' \\\n'.join(args_list),
             out_file_path,
             ' '.join(task_cfg_featurecounts['in_file_path_list']),
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

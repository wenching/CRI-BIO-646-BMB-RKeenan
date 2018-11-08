"""RSeQC"""
import os
import sys
import logging
import datetime

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def clipping_profile(args, sw_cfg, task_cfg):
    """
    RSeQC::clipping_profile
    
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

    sw_cfg_rseqc = sw_cfg['rseqc']

    module = "clipping_profile.py"
    task_cfg_module = task_cfg[module]
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_module else '''&>> {}'''.format(task_cfg_module['log_file_path'])

    task_cfg_reference = task_cfg['references'][next(iter(task_cfg['references']))]

    shell_script_path = '' if 'shell_script_path' not in task_cfg_module else task_cfg_module['shell_script_path']

    cmd_list = []

    if sw_cfg_rseqc['module'] is not None:
        cmd_list.append(sw_cfg_rseqc['module'])

    cmd_list.append('''{} \\
--input-file={} \\
--out-prefix={} \\
--mapq={} \\
--sequencing="{}" \\
{}
\n'''.format(
         module,
         task_cfg_module['in_file_path_list'][0],
         task_cfg_module['out_file_base'],
         30,
         task_cfg_module['sequencing_layout'],
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


def geneBody_coverage(args, sw_cfg, task_cfg):
    """
    RSeQC::geneBody_coverage
    
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

    sw_cfg_rseqc = sw_cfg['rseqc']

    module = "geneBody_coverage.py"
    task_cfg_module = task_cfg[module]
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_module else '''&>> {}'''.format(task_cfg_module['log_file_path'])

    task_cfg_reference = task_cfg['references'][next(iter(task_cfg['references']))]

    shell_script_path = '' if 'shell_script_path' not in task_cfg_module else task_cfg_module['shell_script_path']

    cmd_list = []

    if sw_cfg_rseqc['module'] is not None:
        cmd_list.append(sw_cfg_rseqc['module'])

    cmd_list.append('''{} \\
--input={} \\
--refgene={} \\
--minimum_length={} \\
--format={} \\
--out-prefix={} \\
{}
\n'''.format(
         module,
         task_cfg_module['in_file_path_list'][0],
         task_cfg_reference['anno_bed'],
         100,
         ['pdf', 'png', 'jpeg'][0],
         task_cfg_module['out_file_base'],
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



def infer_experiment(args, sw_cfg, task_cfg):
    """
    RSeQC::infer_experiment
    
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

    sw_cfg_rseqc = sw_cfg['rseqc']

    module = "infer_experiment.py"
    task_cfg_module = task_cfg[module]
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_module else '''&>> {}'''.format(task_cfg_module['log_file_path'])

    task_cfg_reference = task_cfg['references'][next(iter(task_cfg['references']))]

    shell_script_path = '' if 'shell_script_path' not in task_cfg_module else task_cfg_module['shell_script_path']

    cmd_list = []

    if sw_cfg_rseqc['module'] is not None:
        cmd_list.append(sw_cfg_rseqc['module'])

    cmd_list.append('''{} \\
--input-file={} \\
--refgene={} \\
--sample-size={} \\
--mapq={} \\
> {}
\n'''.format(
         module,
         task_cfg_module['in_file_path_list'][0],
         task_cfg_reference['anno_bed'],
         2000,
         30,
         task_cfg_module['out_file_path_list'][0]
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



def RPKM_saturation(args, sw_cfg, task_cfg):
    """
    RSeQC::RPKM_saturation
    
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

    sw_cfg_rseqc = sw_cfg['rseqc']

    module = "RPKM_saturation.py"
    task_cfg_module = task_cfg[module]
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_module else '''&>> {}'''.format(task_cfg_module['log_file_path'])
    strand_rule_cmd = '' if task_cfg_module['strand_specificity'] is None else """--strand='{}'""".format(task_cfg_module['strand_specificity'])
    shell_script_path = '' if 'shell_script_path' not in task_cfg_module else task_cfg_module['shell_script_path']

    task_cfg_reference = task_cfg['references'][next(iter(task_cfg['references']))]


    cmd_list = []

    if sw_cfg_rseqc['module'] is not None:
        cmd_list.append(sw_cfg_rseqc['module'])

    cmd_list.append('''{} \\
--input-file={} \\
--out-prefix={} \\
--refgene={} \\
{} \\
--percentile-floor={} \\
--percentile-ceiling={} \\
--percentile-step={} \\
--rpkm-cutoff={} \\
--mapq={} \\
{}
\n'''.format(
         module,
         task_cfg_module['in_file_path_list'][0],
         task_cfg_module['out_file_base'],
         task_cfg_reference['anno_bed'],
         strand_rule_cmd,
         5,
         100,
         5,
         0.01,
         30,
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


"""STAR"""
import os
import sys
import logging
import datetime
import re

import lib.sambamba
import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def star(args, sw_cfg, task_cfg):
    """
    STAR
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {picard, rseqc}
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

    sw_cfg_star = sw_cfg['star']

    task_cfg_main = task_cfg['main']

    task_cfg_star = task_cfg['star']
    out_file_path = task_cfg_star['out_file_path_list'][0]
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_star else '''&>> {}'''.format(task_cfg_star['log_file_path'])

    task_cfg_ref = task_cfg['references']

    shell_script_path = '' if 'shell_script_path' not in task_cfg_star else task_cfg_star['shell_script_path']

    cmd_list = []

    if sw_cfg_star['module'] is not None:
        cmd_list.append(sw_cfg_star['module'])

    cmd_list.append('''{} \\
--runMode alignReads \\
--genomeLoad {} \\
--outFileNamePrefix {} \\
--readFilesCommand {} \\
--genomeDir {} \\
--readFilesIn {} \\
--runThreadN {} \\
--outSAMstrandField {} \\
--outFilterIntronMotifs {} \\
--outSAMtype {} \\
--outReadsUnmapped {} \\
--outFilterMultimapNmax {} \\
--twopassMode {} \\
{}'''.format(
             sw_cfg_star['exe'],
             sw_cfg_star['genome_load'],
             re.sub("bam$", "", out_file_path),
             sw_cfg_star['readFiles_command'],
             task_cfg_ref[next(iter(task_cfg_ref))]['star_index'],
             ' '.join(task_cfg_main['in_file_path_list']),
             args.threads,
             sw_cfg_star['out_sam_strand_field'],
             sw_cfg_star['out_filter_intron_motifs'],
             sw_cfg_star['out_sam_type'],
             sw_cfg_star['out_reads_unmapped'],
             sw_cfg_star['max_multiple_hit'],
             sw_cfg_star['two_pass_mode'],
             log_file_path_cmd             
         ))


    out_file_substr = "Aligned.sortedByCoord.out.bam" if 'SortedByCoordinate' in sw_cfg_star['out_sam_type'] else "Aligned.out.bam"

    cmd_list.append('''ln -sf \\
{} \\
{} \\
{}'''.format(
             os.path.basename(re.sub("bam$", out_file_substr, out_file_path)),
             out_file_path,
             log_file_path_cmd             
         ))

    sw_cfg_sambamba = util.ddictfunc.subset(sw_cfg, ['sambamba'])

    task_cfg_sambamba = \
        {'sambamba':
            {
                'in_file_path_list': [out_file_path],
                'out_file_path_list': ["""{}.bai""".format(out_file_path)],
                'log_file_path': None if 'log_file_path' not in task_cfg_star else task_cfg_star['log_file_path']
            }            
        }

    cmd_list.extend(
        lib.sambamba.index(
            args,
            sw_cfg_sambamba,
            task_cfg_sambamba
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


    logging.debug("[ {} ] Make Shell Script - DONE\n".format(SELF_FILE))

    return cmd_list

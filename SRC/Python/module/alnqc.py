"""Alignment QC"""
import os
import sys
import logging
import datetime

import lib.picard
import lib.rseqc

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def aln_qc(args, sw_cfg, task_cfg):
    """
    Alignment QC
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {picard|resqc}
    :parm task_cfg: a dictionary of corresponding task configureation {main, picard|resqc, references}
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
    
    for alnqctool in task_cfg.keys():
        if alnqctool in ["main", "references"]: continue
        
        logging.info("[ {} ] Aligment QC Tool: {}\n".format(SELF_FILE, alnqctool))
        
        module_list = task_cfg[alnqctool]['main']['module']
        
        logging.info(
            "[ {} ] Alignment QC Tool: {} with module(s) {}".format(
                SELF_FILE,
                alnqctool,
                ', '.join(module_list)
            )
        )
        
        for aligner in task_cfg[alnqctool].keys():
            if aligner in ["main"]: continue
            
            logging.info("[ {} ] Aligner: {}\n".format(SELF_FILE, aligner))
    
            for library in task_cfg[alnqctool][aligner].keys():
                if library in ['main']: continue
                
                logging.info("[ {} ] Library: {}\n".format(SELF_FILE, library))

                out_dir_t_path = task_cfg[alnqctool][aligner][library]['main']['out_dir_path']

                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ $out_dir_t_path ]")
                    os.makedirs(out_dir_t_path, exist_ok=True)

                out_dir_t_path = task_cfg[alnqctool][aligner][library]['main']['tmp_dir_path']

                if not os.path.exists(out_dir_t_path):
                    logging.debug("MKDIR: [ $out_dir_t_path ]")
                    os.makedirs(out_dir_t_path, exist_ok=True)                    

                for module in module_list:
                    logging.info(
                        "[ {} ] Alignment QC Tool::Module: {}::{}".format(
                            SELF_FILE,
                            alnqctool,
                            module
                        )
                    )

                    task_cfg_lb = {
                        'references': task_cfg['references'],
                        module: task_cfg['main'][aligner][library]
                    }                    
                    task_cfg_lb[module].update(task_cfg[alnqctool][aligner][library]['main'])
                    task_cfg_lb[module].update(task_cfg[alnqctool][aligner][library][module])
                    
                    #util.ddictfunc.pprint_ddicts(task_cfg_lb)
                    '''
                    'CollectRnaSeqMetrics': {'in_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Aln/star/KO01/KO01.star.bam'],
                                              'log_file_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/run.alnQC.picard.star.KO01.CollectRnaSeqMetrics.log',
                                              'out_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01',
                                              'out_file_base': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard',
                                              'out_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics',
                                                                     '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics.pdf'],
                                              'shell_script_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/Shell/run.alnQC.picard.star.KO01.CollectRnaSeqMetrics.sh',
                                              'strand_specificity': 'NONE',
                                              'tmp_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/AlnQC/picard/star/KO01/tmp'},
                     'references': {'grch38': {'anno_bed': None,
                                               'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.bed12',
                                               'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.gtf',
                                               'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.refFlat.txt',
                                               'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.bed',
                                               'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.interval_list',
                                               'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                               'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                               'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.fa',
                                               'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.dict',
                                               'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11'}}}
                    
                    '''

                    sw_cfg_tool = util.ddictfunc.subset(sw_cfg, [alnqctool])

                    if(alnqctool == "picard"):
                        if(module == "CollectRnaSeqMetrics"):
                            lib.picard.collect_rna_seq_metrics(args, sw_cfg_tool, task_cfg_lb)
                    elif(alnqctool == "rseqc"):
                        if(module == "clipping_profile.py"):
                            lib.rseqc.clipping_profile(args, sw_cfg_tool, task_cfg_lb)
                        elif(module == "geneBody_coverage.py"):
                            lib.rseqc.geneBody_coverage(args, sw_cfg_tool, task_cfg_lb)
                        elif(module == "infer_experiment.py"):
                            lib.rseqc.infer_experiment(args, sw_cfg_tool, task_cfg_lb)
                        elif(module == "RPKM_saturation.py"):
                            lib.rseqc.RPKM_saturation(args, sw_cfg_tool, task_cfg_lb)


                        """
                        RSeQC::bam_stat.py:               Summarizing mapping statistics of a BAM or SAM file
                        [Y] RSeQC::clipping_profile.py:   Calculate the distributions of clipped nucleotides across reads
                        [Y] RSeQC::geneBody_coverage.py:  Calculate the RNA-seq reads coverage over gene body
                        [Y] RSeQC::infer_experiment.py:   Used to “guess” how RNA-seq sequencing were configured, particulary how reads were stranded for strand-specific RNA-seq data, through comparing the “strandness of reads” with the “standness of transcripts”
                        RSeQC::junction_annotation.py:    Compare detected splice junctions to reference gene model
                        RSeQC::junction_saturation.py:    Check for saturation by resampling total alignments and then detects splice junctions from each subset and compares them to reference gene model.
                        RSeQC::read_distribution.py:      Calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions)
                        RSeQC::read_GC.py:                GC content distribution of reads
                        RSeQC::read_NVC.py:               Check the nucleotide composition bias
                        RSeQC::read_quality.py:           Box plot and heat map of read quality per base
                        [Y] RSeQC::RPKM_saturation.py:    Resample a series of subsets from total RNA reads and then calculate RPKM value using each subset
                        """

                    logging.debug(
                        "[ {} ] Alignment QC Tool::Module: {}::{} - DONE\n".format(
                            SELF_FILE,
                            alnqctool,
                            module
                        )
                    )
    
                logging.debug("[ {} ] Library: {} - DONE\n".format(SELF_FILE, library))
            logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))
        logging.debug("[ {} ] Aligment QC Tool: {} - DONE\n".format(SELF_FILE, alnqctool))

    return

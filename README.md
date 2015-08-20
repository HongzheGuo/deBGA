###Introduction
deBGA is a seed-and-extension-based read alignment tool. It is suitable for aligning various kinds of high-throughput sequencing reads to multiple similar reference genomes.

deBGA indexes the genome through de Bruijn graph framework, which captures the genomic substrings implied by the unipaths and indexes them with a hash table-based index. With the index, deBGA adopts a seed-and-extension strategy. In seeding phase, a series of seeds from the reads are matched to the reference genome inferring a set of putative read positions (PRPs) from the matched positions. In the extension phase, the reads are aligned with the local sequences around the PRPs to compose full read alignments.

deBGA has outstanding throughput on aligning reads from various prokaryote and eukaryote genomes. A series of benchmarks on model organism genomes, e.g., E. coli, H. sapiens, etc., demonstrated that it can simultaneously achieve good throughput, sensitivity and accuracy in various kinds of read alignment tasks. deBGA is open source and free for non-commercial use.

deBGA is mainly designed by Bo Liu and developed by Hongzhe Guo in Center for Bioinformatics, Harbin Institute of Technology, China.

###Memory requirement
The memory usage of deBGA can fit the configurations of most modern servers and workstations. Its peak memory footprint depends on the length of reference genome, i.e., 25.30 Gigabytes and 5.20 Gigabytes respectively for the real H. Sapiens and E. Coli genomes, on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04.

###Installation
Current version of deBGA needs to be run on Linux operating system.  
The source code is written in C, and can be directly download from: https://github.com/HongzheGuo/deBGA  
The makefile is attached. Use the make command for generating the executable file.  

###Synopsis
deBGA index [options] reference.fasta \<index_route\>  
Index reference in RdBG-Index format  

deBGA aln [options] \<index_route\> \<single_end_read.fastq [pair_end_read1.fastq pair_end_read2.fastq]\> \<result_file.sam\>  
Align read to its primitive location in Reference  

###Parameters (could be updated in the future for adding new functions)
```
deBGA index   
-k,                     INT the k-mer size of the vertices of Reference de Bruijn Graph (RdBG)[22]. 
    
deBGA aln 
-k,                     INT the minimum length of a valid Uni-MEM seed. For current version of deBGA, this setting should be equal to the k-mer size of the RdBG-index[22].    
-s,                     INT the number of rounds of re-seeding. deBGA iteratively aligns a read in at most (-s + 1) rounds with various set of seeds.   
                            This option works combining with the minimum interval of seeding (the -i option) and the maximum allowed number of hits per seed (the -n option).  
                            That is, in the r-th round (r = 1 ,…, -s), deBGA tries to generate seeds at every (-s – r +1)*(-i) bp along the read.   
                            If the read still cannot be successfully aligned after -s rounds, deBGA would ignore -n option to handle very repetitive reads in the (-s+1)-th round[4].    
-i,                     INT the minimum interval of seeding. This parameter determines the density of seeds, which is related to the sensitivity and efficiency of alignment.   
                            Configuring this parameter lower can make deBGA generate seed more densely, which could benefit the sensitivity, but at the expense of throughput.  
                            The default value of this parameter is determined by a series of benchmarking on various datasets[5].    
-n,                     INT the maximum allowed number of hits per seed. In the first -s rounds of the alignment process, the seeds with more than -n hits would be discarded for improving the efficiency.  
                            deBGA would ignore this restriction to introduce repetitive seeds if the read still cannot be successfully aligned after -s rounds[300].    
-c,                     NUM the threshold of confident alignment. In each round, deBGA would check the edit distance of the best obtained alignment.  
                            If the ratio ED_best/RL < -c, where ED_best and RL respectively indicate the edit distance of the best obtained alignment and the length of the read,  
                            deBGA considers that the read is confidently aligned and early-stops the alignment process[0.05].    
--cl,                   NUM the adjusted threshold of confident alignment. For a given round, if there is Uni-MEM seed available for extension, but no successful alignment is obtained,  
                            the threshold of confident alignment can be dynamically adjusted in next rounds. This is a heuristic for acceleratingthe alignment of reads with many changes,   
                            e.g., reads having many sequencing errors. If --cl is not set, there will be no change on the -c option during the process (default: following the setting of -c).    
--local,                    the local alignment option for confident alignment. If --local option is set, deBGA will perform local alignment instead of end-to-end alignment in next rounds   
                            when the read cannot be successfully aligned with valid Uni-MEM seeds. The best obtained local alignment will be output as the result.   
                            It is also worthnoting that the --cl option and --local option should not be simultaneously set (default: not set).    
--local-match,          INT the score for a matched base in the local alignment.This option will take effect only if --local option is set[1].     
--local-mismatch,       INT the penalty for a mismatched base in the local alignment. This option will take effect only if --local option is set[4].    
--local-gap-open,       INT the penalty for a gap open in the local alignment.This option will take effect only if --local option is set[6].    
--local-gap-extension,  INT the penalty for gap extension in the local alignment.This option will take effect only if --local option is set[1].     
-u,                     INT the upper limit of insert size. For a pair-end read, deBGA pairs the alignments of the two ends according to the upper (-u option) and lower (-f option) limits of the insert size.   
                            deBGA will consider it as a suitable pair-end alignment only if the distance between the aligned positions of the two read ends are within the range [-f, -u][1000].    
-f,                     INT ower limit of insert size. For a pair-end read, deBGA pairs the alignments of the two ends according to the upper (-u option) and lower (-f option) limits of the insert size.  
                            deBGA will consider it as a suitable pair-end alignment only if the distance between the aligned positions of the two read ends are within the range [-f, -u][50].    
-e,                     INT budget for single-end alignment. In single-end read alignment, deBGA sets a budget on the computation resourcein advancefor balancing the efficiency and the sensitivity.  
                            More precisely, in the extension phase, deBGA subsequentlyextend the candidate seedsin order oftheir coverage lengths,until more than -e extension operations have been totally   
                             executed after extending some of the seeds, or all the seeds are extended[100].    
-o,                     INT the maximum number of alignment output.deBGA outputs at most -o alignments for the read. This is except for the pair-end reads which are handled with the anchoring alignment strategy.    
                            For such reads, the number of output alignments is determined by the -x option[20].  
-x,                     INT the maximum number of alignment output for anchoring alignment. For the pair-end reads aligned with the anchoring alignment strategy, deBGA will output at most -x alignments[150].    
-l,                     INT the maximum allowed read length. For the current version of deBGA, reads shorter than -l bp will be normally processed, andfor reads longer than -l bp, only the first -l bp will be aligned,  
                            and the other parts will be trimmed. Set -l option with a larger number may slightly increasethe RAM usage. For most nowadays next generation sequencing reads, e.g., reads from Illumina platforms,   
                            the default setting is long enough to work withoutthe trimming. Moreover, the current version ofdeBGA can support reads upto 4096 bp (setting -l to 4096)[512].     
-p,                     INT the number of threads. The current version of deBGA supports upto 32 threads in read alignment[1].    
```

###Quick start
Genome indexing:
deBGA index Reference Index_Dir
Read alignment:
deBGA aln Index_Dir Fastq_File Sam_file

###Simulation benchmarking
We simulated a series of datasets from various genomes, i.e., human genome build GRCh37/hg19 and the reference sequences of the 62 E. Coli strains, through Mason Simulator (version0.1.2). For human genome, five read lengths from 100 bp to 250 bp were used and the mean and standard deviation of the insert size are respectively 500 bp and 25 bp for all the datasets. For the E.coli strains, 100 bp Illumina-like pair-end reads were simulated for evaluation. These datasets helped us to evaluate the performance of deBGA. The datasets have been uploaded to Google Drive, and can be downloaded through the following link:


###Reference
deBGA: Read Alignment with de Bruijn Graph-based Seed and Extension. Manuscript in preparation.

###Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn


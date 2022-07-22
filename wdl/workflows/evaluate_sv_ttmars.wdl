version 1.0

import "../tasks/ttmars.wdl" as ttmars_t

workflow evaluateSVsWithTTMars {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Evaluate structural variants against truth assemblies using TT-Mars (https://github.com/ChaissonLab/TT-Mars). If necessary, lift-over files can be made from the reference genome and the two assembled haplotypes. "
    }

    input {
        File referenceFile           # reference fasta (.fa or .fa.gz)
        File svsFile                 # VCF with the SVs to evaluate (.vcf or .vcf.gz)
        File hap1File                # truth assembly for haplotype 1 (.fa or .fa.gz)
        File hap2File                # truth assembly for haplotype 1 (.fa or .fa.gz)
        File trfFile                 # BED file listing all the regions in the reference with simple repeats
        File centromereFile          # centromere location of the reference (see https://github.com/ChaissonLab/TT-Mars/blob/main/centromere_hg38.txt)
		Int threads=16               # how many threads to use for the LRA alignment (when lifting over is needed)
        Int nbXchr=2                 # basically specify if sample is male (nbXchr=1) or female (nbXchr=2)
        File? nonCovReg1File         # lift-over coverage of haplotype 1 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh
        File? nonCovReg2File         # lift-over coverage of haplotype 2 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh
        File? loPosAssem1File        # lift-over information of reference on haplotype 1. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh
        File? loPosAssem2File        # lift-over information of reference on haplotype 2. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh
        File? loPosAssem10File       # lift-over information of haplotype 1 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh
        File? loPosAssem20File       # lift-over information of haplotype 2 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh
    }

    ## if any of those input files are missing, the assemblies have to be lifted to the reference sequence
    if(!defined(nonCovReg1File) || !defined(nonCovReg2File) || !defined(loPosAssem1File) || !defined(loPosAssem2File) || !defined(loPosAssem10File) || !defined(loPosAssem20File)){
        ## align each haplotype to the reference with LRA
        call ttmars_t.lra_t as lra1 {
            input:
            reference_file=referenceFile,
            hap_file=hap1File,
            threads=threads
        }
        call ttmars_t.lra_t as lra2 {
            input:
            reference_file=referenceFile,
            hap_file=hap2File,
            threads=threads
        }
        ## make liftover files for TTmars
        call ttmars_t.liftover_t {
            input:
            hap1_lra_file=lra1.bamOutput,
            hap2_lra_file=lra2.bamOutput
        }
    }

    File non_cov_reg_1_file = select_first([nonCovReg1File, liftover_t.non_cov_reg_1_file])
    File non_cov_reg_2_file = select_first([nonCovReg2File, liftover_t.non_cov_reg_2_file])
    File lo_pos_assem1_file = select_first([loPosAssem1File, liftover_t.lo_pos_assem1_file])
    File lo_pos_assem2_file = select_first([loPosAssem2File, liftover_t.lo_pos_assem2_file])
    File lo_pos_assem1_0_file = select_first([loPosAssem10File, liftover_t.lo_pos_assem1_0_file])
    File lo_pos_assem2_0_file = select_first([loPosAssem20File, liftover_t.lo_pos_assem2_0_file])

    ## call TT-Mars
    call ttmars_t.ttmars_t {
        input:
        reference_file=referenceFile,
        svs_file=svsFile,
        hap1_file=hap1File,
        hap2_file=hap2File,
        centromere_file=centromereFile,
        trf_file=trfFile,
        non_cov_reg_1_file=non_cov_reg_1_file,
        non_cov_reg_2_file=non_cov_reg_2_file,
        lo_pos_assem1_file=lo_pos_assem1_file,
        lo_pos_assem2_file=lo_pos_assem2_file,
        lo_pos_assem1_0_file=lo_pos_assem1_0_file,
        lo_pos_assem2_0_file=lo_pos_assem2_0_file,
        nb_x_chr=nbXchr
    }

	output {
        File tsvOutput = ttmars_t.tsvOutput
        File? liftoverNonCovReg1File = liftover_t.non_cov_reg_1_file
        File? liftoverNonCovReg2File = liftover_t.non_cov_reg_2_file
        File? liftoverLoPosAssem1File = liftover_t.lo_pos_assem1_file
        File? liftoverLoPosAssem2File = liftover_t.lo_pos_assem2_file
        File? liftoverLoPosAssem10File = liftover_t.lo_pos_assem1_0_file
        File? liftoverLoPosAssem20File = liftover_t.lo_pos_assem2_0_file
	}
}

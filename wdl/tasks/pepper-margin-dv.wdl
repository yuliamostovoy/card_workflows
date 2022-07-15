version 1.0

task pepper_margin_dv_t {
  input {
    Int threads
    File reference
	File bamAlignment
	File bamAlignmentIndex
	String mapMode = "ont"
	Int memSizeGb = 256
	Int diskSizeGb = 1024
  }

  String pepperMode = if mapMode == "ont" then "--ont_r9_guppy5_sup" else "--hifi"

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    ## unzip reference fasta is necessary
    REF=~{reference}
    if [[ $REF == *.gz ]]
    then
        gunzip $REF -c > ref.fa
        REF=ref.fa
    fi
    ## index reference
    samtools faidx $REF

    ## make sure both bam and bai are "in" the same directory
    ln -s ~{bamAlignment}
    ln -s ~{bamAlignmentIndex}
    BAM_FILE=$(basename ~{bamAlignment})

    ## run PEPPER-Margin-DeepVariant
    run_pepper_margin_deepvariant call_variant -b $BAM_FILE -f $REF -o `pwd` -t ~{threads} ~{pepperMode} -p PMDV_FINAL
  >>>

  output {
	File pepperVcf = "PMDV_FINAL.vcf.gz"
    #File haplotaggedBam = "intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam"
  }

  runtime {
    docker: "kishwars/pepper_deepvariant:r0.8"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}

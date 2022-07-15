version 1.0

task pepper_margin_dv_t {
  input {
    Int threads
    File reference
	File bamAlignment
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

    ## index reference and bam
    samtools faidx $REF
    samtools index -@ ~{threads} ~{bamAlignment}

    ## run PEPPER-Margin-DeepVariant
    run_pepper_margin_deepvariant call_variant -b ~{bamAlignment} -f $REF -o `pwd` -t ~{threads} ~{pepperMode} -p PMDV_FINAL
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

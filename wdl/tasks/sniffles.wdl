version 1.0

task sniffles_t {
  input {
	File bamAlignment
    File reference
    Int threads = 8
	Int memSizeGb = 128
	Int diskSizeGb = 256
	String trfAnnotations = ""
  }
  
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    TRF_STRING=""
    if [ ! -z ~{trfAnnotations} ]
    then
       TRF_STRING="--tandem-repeats /opt/trf_annotations/~{trfAnnotations}"
    fi
    echo $TRF_STRING

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
    
    sniffles -i ~{bamAlignment} -v sniffles.vcf.gz --reference $REF -t ~{threads} ${TRF_STRING}
  >>>

  output {
	File snifflesVcf = "sniffles.vcf.gz"
  }

  runtime {
    docker: "mkolmogo/card_sniffles:2.0.3"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}

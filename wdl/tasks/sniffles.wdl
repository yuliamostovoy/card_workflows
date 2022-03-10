version 1.0

task sniffles_t {
  input {
    Int threads
	File bamAlignment
	Int memSizeGb = 128
	Int diskSizeGb = 256
  }
  
  File snifflesOut = "sniffles.vcf"
  
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    samtools index -@ 10 ~{bamAlignment}
    sniffles -i ~{bamAlignment} -v ~{snifflesOut} -t ~{threads}
  >>>

  output {
	File snifflesVcf = snifflesOut
  }

  runtime {
    docker: "mkolmogo/card_sniffles:2.0.3"
    cpu: threads
	memory: memSizeGb + " GB"
	disks: "local-disk " + diskSizeGb + " SSD"
  }
}

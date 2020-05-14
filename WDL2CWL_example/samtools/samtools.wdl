workflow samtools {

 
  File referenceGenome
  String name
  String samtoolsPath
  String bcftoolsPath
 
  call Alignment {
    input:
         ReferenceGenome=referenceGenome,
         sampleName=name,
         index=name
}

 call samtoolsView {
    input:
         inputSAM=Alignment.rawSAM,
         samtoolsPath = samtoolsPath,
         sampleName=name
}

 call samtoolsSort {
    input:
         inputBAM=samtoolsView.rawBAM,
         samtoolsPath = samtoolsPath,
         sampleName=name
    
}

 call samtoolsIndex {
    input:
        inputBAM=samtoolsSort.rawBAM,
        samtoolsPath = samtoolsPath,
        sampleName=name
}


 call bcftools {
   input:
        inputBAM=samtoolsSort.rawBAM,
        sampleName=name,
        bcftoolsPath = bcftoolsPath,
        ReferenceGenome=referenceGenome
       
} 
        }

task Alignment {
 
  File leftFastq
  File rightFastq
  File ReferenceGenome
  String sampleName
  String index


  command {
           
          export PATH=$PATH:/home/ngsap2/Downloads/bowtie2-2.4.1

          bowtie2-build ${ReferenceGenome} ${index} 
          
          bowtie2 -q -x ${index} -1 ${leftFastq} -2 ${rightFastq} -S ${sampleName}.sam           
}


output {
    
    File rawSAM = "${sampleName}.sam"
  }


}

task samtoolsView {
 
 File inputSAM
 String sampleName
 String samtoolsPath
 
 command {
            samtools view -bS ${inputSAM} > ${sampleName}.bam
}


output {
 
 File rawBAM = "${sampleName}.bam"
}


}

task samtoolsSort {

 File inputBAM
 String sampleName
 String samtoolsPath

command {

  samtools sort -l 0 -o ${sampleName}.sorted.bam ${inputBAM}

}

output {
 
 File rawBAM = "${sampleName}.sorted.bam"

}

}

task samtoolsIndex {

 File inputBAM
 String sampleName
 String samtoolsPath

command {

  samtools index ${inputBAM} > ${sampleName}.bam.bai
  
 }


output {
 
 File rawBAM = "${sampleName}.bam.bai"

}

}

task bcftools
{

File inputBAM
String sampleName
File ReferenceGenome
String bcftoolsPath

command {

bcftools mpileup -f ${ReferenceGenome} ${inputBAM} | bcftools call -mv -Ob -o ${sampleName}.vcf
 
}

output{

File rawVCF = "${sampleName}.vcf"

}
 
}



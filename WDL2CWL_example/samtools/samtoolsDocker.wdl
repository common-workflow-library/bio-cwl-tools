workflow samtoolsDocker     
{

  String bowtieDocker
  String samtoolsDocker
  String bcftoolsDocker 
  File referenceGenome
  String name
 
  call Alignment {
    input:
         ReferenceGenome=referenceGenome,
         sampleName=name,
         docker = bowtieDocker,
         index=name
}

 call samtoolsView {
    input:
         inputSAM=Alignment.rawSAM,
         docker = samtoolsDocker,
         sampleName=name
}

 call samtoolsSort {
    input:
         inputBAM=samtoolsView.rawBAM,
         docker = samtoolsDocker,
         sampleName=name
    
}

 call samtoolsIndex {
    input:
        inputBAM=samtoolsSort.rawBAM,
        docker = samtoolsDocker,
        sampleName=name
}


 call bcftools {
   input:
        inputBAM=samtoolsSort.rawBAM,
        sampleName=name,
        ReferenceGenome=referenceGenome,
        docker=bcftoolsDocker
} 
        }

task Alignment {
 
  File leftFastq
  File rightFastq
  File ReferenceGenome
  String sampleName
  String index
  String docker

  command {
          
          bowtie2-build ${ReferenceGenome} ${index} 
          
          bowtie2 -q -x ${index} -1 ${leftFastq} -2 ${rightFastq} -S ${sampleName}.sam           
}

runtime
{
  docker:docker
}

output {
    
    File rawSAM = "${sampleName}.sam"
  }


}

task samtoolsView {
 
 File inputSAM
 String sampleName
 String docker
 command {
            samtools view -bS ${inputSAM} > ${sampleName}.bam
}

runtime
{
  docker:docker
}

output {
 
 File rawBAM = "${sampleName}.bam"
}


}

task samtoolsSort {

 File inputBAM
 String sampleName
 String docker

command {

  samtools sort -l 0 -o ${sampleName}.sorted.bam ${inputBAM}

}

runtime
{
  docker:docker
}

output {
 
 File rawBAM = "${sampleName}.sorted.bam"

}

}

task samtoolsIndex {

 File inputBAM
 String sampleName
 String docker

command {

  samtools index ${inputBAM} > ${sampleName}.bam.bai
  
 }

runtime
{
  docker:docker
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
String docker


command {

bcftools mpileup -f ${ReferenceGenome} ${inputBAM} | bcftools call -mv -Ob -o ${sampleName}.vcf
 
}

runtime
{
  docker:docker
}

output{

File rawVCF = "${sampleName}.vcf"

}
 
}



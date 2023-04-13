#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: ['java', '-jar', '/usr/src/myapp/MetFragCommandLine-2.5.0.jar']

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/maw-metfrag_2.5.0:1.0.4
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - entryname: param_file.txt
      entry: |-
        PeakListPath = $(inputs.PeakList.path)
        IonizedPrecursorMass =  $(inputs.IonizedPrecursorMass)
        PrecursorIonMode = $(inputs.PrecursorIonMode)
        IsPositiveIonMode = $(inputs.IsPositiveIonMode)
        MetFragDatabaseType = LocalCSV
        LocalDatabasePath = $(inputs.LocalDatabase.path)
        DatabaseSearchRelativeMassDeviation = 5
        FragmentPeakMatchAbsoluteMassDeviation = 0.001
        FragmentPeakMatchRelativeMassDeviation = 15
        MetFragCandidateWriter = CSV
        SampleName = $(inputs.SampleName)
        ResultsPath = .
        MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter
        MetFragPostProcessingCandidateFilter = InChIKeyFilter
        MaximumTreeDepth = 2
        NumberThreads = $(runtime.cores)


inputs: # additional inputs for all files; make them to show certain paths
  PeakList: File
  IonizedPrecursorMass: string
  PrecursorIonMode: int
  IsPositiveIonMode: boolean
  LocalDatabase: File
  SampleName: string
  #MetFragDatabaseType: string
  #DatabaseSearchRelativeMassDeviation: int
  #FragmentPeakMatchAbsoluteMassDeviation: float
  #FragmentPeakMatchRelativeMassDeviation: float
  #MetFragCandidateWriter: string
  #MetFragPreProcessingCandidateFilter: string
  #MetFragPostProcessingCandidateFilter: string
  #MaximumTreeDepth: int
  #NumberThreads: int

arguments:
  - ParameterFile=param_file.txt

outputs:
  candidate_list:
    type: File
    outputBinding:
        glob: "*.csv"

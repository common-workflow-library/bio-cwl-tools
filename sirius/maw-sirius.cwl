cwlVersion: v1.0 
class: CommandLineTool

$namespaces:
  cwltool: http://commonwl.org/cwltool#
hints:
  "cwltool:Secrets":
    secrets:
      - sirius_user
      - sirius_password

requirements:
  DockerRequirement:
    dockerPull: docker.io/zmahnoor/maw-sirius5:dev
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  EnvVarRequirement:
    envDef:
      SIRIUS_USER: $(inputs.sirius_user)
      SIRIUS_PASSWORD: $(inputs.sirius_password)

inputs: 
  sirius_user:
    type: string
  sirius_password:
    type: string
  spectrum:
    type: File
#    inputBinding: 
#      prefix: --input
  profile:
    type: string
    default: "orbitrap"
  candidates:
    type: int
    default: 30
  ppm_max:
    type: int
    default: 5
  ppm_max_ms2:
    type: int
    default: 15
  database:
    type: string
    default: "coconut"
  isotope:
    type: boolean
    default: False


arguments:
    - sirius
    - login
    - --user-env 
    - SIRIUS_USER
    - --password-env 
    - SIRIUS_PASSWORD 
    - shellQuote: False
      valueFrom: ";"
    - sirius
    - --input
    - $(inputs.spectrum.path)
    - --output
    - $(inputs.spectrum.nameroot).json
    - formula
    - --profile
    - $(inputs.profile)
    - | 
      ${
        if (inputs.isotope) { 
          return ["--no-isotope-filter",
            "--no-isotope-score"];
        } else {
          return null;
        }
        
      }
    - --candidates
    - $(inputs.candidates)
    - --ppm-max
    - $(inputs.ppm_max)
    - --ppm-max-ms2
    - $(inputs.ppm_max_ms2)
    - fingerprint
    - structure
    - --database
    - $(inputs.database)
    - compound-classes
    - write-summaries
    - --output
    - $(inputs.spectrum.nameroot).json

outputs:
  results: 
    type: Directory
    outputBinding:
       glob: $(inputs.spectrum.nameroot).json
  
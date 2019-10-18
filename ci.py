from subprocess import check_output, check_call
import sys

has_failed = False
changed_files = check_output("git --no-pager diff --name-status origin/release..$(git branch | grep \* | cut -d ' ' -f2)", shell=True)
for line in changed_files.decode('utf-8').rstrip().split('\n'):
    change = line.split('\t')[0]

    if change not in ['M', 'A', 'C', 'R', 'X']:
        continue
    fs = line.split('\t')[1]
    print(fs)

    if not fs.lower().endswith(".cwl"):
        continue

    file_validation_status = check_call("cwltool --validate ")
    if file_validation_status != 0:
        print(f'Tool Failed: {fs}')
        has_failed = True
if has_failed == True:
    sys.exit(1)

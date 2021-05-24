import crayons
import os
from subprocess import check_output, check_call
import sys

has_failed = False
tool_failed = False
# changed_files = check_output("git --no-pager diff --name-status release..$(git branch | grep \* | cut -d ' ' -f2)", shell=True)
try:
    changed_files = check_output(
        "git --no-pager diff --name-status --diff-filter d release..FETCH_HEAD",
        shell=True,
    )
    print("Using Head")
except:
    print("Using Master")
    changed_files = check_output("git --no-pager diff --name-status HEAD", shell=True)

for line in changed_files.decode("utf-8").rstrip().split("\n"):
    tool_failed = False
    change = line.split("\t")[0]

    if change not in ["M", "A", "C", "R", "X"]:
        continue
    fs = line.split("\t")[1]

    if not fs.lower().endswith(".cwl"):
        continue
    print(crayons.blue(f"Testing CWL Validation: {fs} \n"))
    file_validation_status = check_call(f"cwltool --validate {fs}", shell=True)
    if file_validation_status != 0:
        print(f"Tool Failed Validation: {fs}")
        tool_failed = True

    if file_validation_status == 0:
        print(crayons.green(f"Tool Passed Validation: {fs}\n"))

    print(crayons.blue(f"Testing Repo Requirements\n"))
    with open(fs) as f:
        lines_expected = [
            ["#!/usr/bin/env cwl-runner"],
            [
                "cwlVersion:v1.0",
                "cwlVersion: v1.0",
                "cwlVersion:v1.1",
                "cwlVersion: v1.1",
            ],
            ["class: CommandLineTool", "class: Workflow", "class: ExpressionTool"],
        ]
        for index, line in enumerate(f.readlines()):
            line = line.rstrip()
            if index == 3:
                break
            patterns = lines_expected[index]
            for line_expected in patterns:
                if line == line_expected:
                    break
            else:
                # if we get to the end of the loop without breaking, it means the line has not been matched
                print(
                    crayons.red(
                        f"Tool Failed Requirements: Line {index + 1} : {fs} : got: '{line}', expected one of {patterns}\n"
                    )
                )
                tool_failed = True
    if tool_failed == False:
        print(crayons.green(f"Tool Passed Repo Requirements: {fs}\n"))
    else:
        print(crayons.red(f"Tool Failed Repo Requirements: {fs}\n"))
        has_failed = True

    if os.access(fs, os.X_OK):
        print(crayons.green(f"Tool Passed Executable Test: {fs}"))
    else:
        print(crayons.red(f"File Failed Executable Requirement: {fs}"))
        has_failed = True
    print("=====================")

if has_failed == True:
    sys.exit(1)

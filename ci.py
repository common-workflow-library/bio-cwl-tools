from subprocess import call

changed_files = call("git --no-pager diff --name-status origin/release..$(git branch | grep \* | cut -d ' ' -f2)", shell=True)
print(changed_files)

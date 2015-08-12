#!/bin/bash

#---automatically push current code to github
make docs clean
rm -r repo
timestamp=$(date +%Y.%m.%d.%H%M)
du -h --max-depth=1
echo "pushing code to github via"
echo
echo "git commit -a -m \""$timestamp" : "${@:2}"\""
echo
read -p "okay to commit and push y/N? " -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
	commit_message=$timestamp" : "${@:2}
	git add . --all
	find . -name "*.pyc" -exec git rm -f "{}" \;
  	git commit -a -m "$commit_message"
    git push
fi
echo


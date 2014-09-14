==Orientation==

You are in one of three:
* Currently, the "origin" github repository https://github.com/yanailab/CEL-Seq-pipeline.git
        Only houses master
* The dev directory /data/tools/pipeline_dev
        Houses dev and master and sub-branches
* The stable directory /data/tools/pipeline_stable
        Only houses master

==Synching across==

After done making changes at pipeline_dev the stable/master branch;
* Sync to github:
        cd /data/tools/pipeline_dev
        git push origin master
* Sync to stable pipeline:
        cd /data/tools/pipeline_stable
        git pull origin master

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
        git push origin master
* Sync to stable pipeline:
        git push /data/tools/pipeline_stable/ master

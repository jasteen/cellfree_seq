# cellfree_seq


```
/usr/local/python/2.7.12_static/bin/virtualenv --system-site-packages pipeline_env
module load drmaa
export DRMAA_LIBRARY_PATH=/usr/local/drmaa/1.0.7/lib/libdrmaa.so
source pipeline_env/bin/activate
pip install git+https://github.com/SoutheyLab/cellfree_seq
cellfree_seq --config pipeline.config --use_threads --log_file pipeline.log --jobs 10 --verbose 3 --just_print
```

The following lines need to be added to your .bash_rofile so that VEP will function correctly.

```
PERL5LIB=$PERL5LIB:/projects/vh83/reference/annotation_databases/VEP_Modules/.vep/Plugins
export PERL5LIB
```




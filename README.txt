DeSNP Dependencies

samtools
pysam
numpy

if you have a modulefile,  like "module load desnp/version" you can modify PYTHONPATH as needed in the modulefile.  I would do something like /opt/compsci/desnp/1.0/ for the scripts and /opt/compsci/desnp/1.0/lib for your python libraries.  the module file will add /opt/compsci/desnp/1.0/ to PATH and /opt/compsci/desnp/1.0/lib to PYTHONPATH

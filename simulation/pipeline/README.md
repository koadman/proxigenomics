First commit of pipeline and ancillaries.

This may not work out of the box since I had a number of paths relative
to my home directory and have also moved the code to a subdir.


##Nestly pipeline

1) run `hic_nestly.py` to prepare the sweep folder.
2) run overlord.sh to compute all the tasks defined in the sweep folder.

The current sweep defintion will create around 80GB of data.

Matt DeMaere


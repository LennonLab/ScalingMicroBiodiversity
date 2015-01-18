import skbio
print skbio.__path__
for i in dir(skbio): print i


from skbio import diversity
print diversity.__path__
for i in dir(diversity): print i

from skbio import stats
print stats.__path__
for i in dir(diversity): print i

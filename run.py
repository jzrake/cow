
from cow import *
import os, sys

print sys.argv

d = cow_domain_new()
f = cow_dfield_new(d, "primitive")

cow_domain_commit(d)
cow_dfield_commit(f)

print cow_dfield_getname(f)

cow_dfield_del(f)
cow_domain_del(d)


import cow
from cow import *
import os, sys

#print sys.argv
#print sys.path

d = cow_domain_new()
f = cow_dfield_new(d, "primitive")

cow_domain_commit(d)
cow_dfield_commit(f)

print test_trans, TEST_TRANS
print cow_dfield_getname(f)

print test_trans(None, None, None, None)
cow_dfield_loop(f, TEST_TRANS, None)
cow_dfield_del(f)
cow_domain_del(d)

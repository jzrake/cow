

from _cow import *
import os, sys
import atexit
import numpy as np

atexit.register(cow_finalize)

cow_init()
domain = cow_domain_new()
dfield = cow_dfield_new(domain, "prim")
cow_dfield_addmember(dfield, "vx")
cow_dfield_addmember(dfield, "vy")
cow_dfield_addmember(dfield, "vz")

cow_domain_setndim(domain, 2);
cow_domain_setguard(domain, 3);
cow_domain_setsize(domain, 0, 10);
cow_domain_setsize(domain, 1, 12);
cow_domain_commit(domain)
cow_dfield_commit(dfield)

print test_trans, TEST_TRANS
print cow_dfield_getname(dfield)
print cow_dfield_getbuffer(dfield)

A = getarray3(dfield)
A[:,:,0] = 1.0
A[:,:,1] = 2.0
A[:,:,2] = 3.0
B = A[3:-3,3:-3,:]
print A.shape, B.shape
print A.flags, B.flags
print A.max(), A.min()

#test_trans(None, None, None, None)
#cow_dfield_loop(dfield, TEST_TRANS, None)
cow_dfield_del(dfield)
cow_domain_del(domain)

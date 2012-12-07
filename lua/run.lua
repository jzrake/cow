
local cow = require "cow"

local modes = 0--bit32.bor(cow.COW_DISABLE_MPI, cow.COW_NOREOPEN_STDOUT)
cow.cow_init(0, nil, modes)


local function DistributedDomain(kwargs)
   local self = { }
   local domain = cow.cow_domain_new()
   local size = kwargs.size
   local guard = kwargs.guard
   cow.cow_domain_setndim(domain, #size)
   for i=1,#size do
      cow.cow_domain_setsize(domain, i-1, size[i])
   end
   cow.cow_domain_setguard(domain, guard)
   cow.cow_domain_commit(domain)
   self._c = domain

   function self:ndim()
      return cow.cow_domain_getndim(domain)
   end
   function self:size()
      local ret = { }
      for i=1,#size do
	 ret[i] = cow.cow_domain_getsize(self._c, i-1)
      end
      return ret
   end

   local meta = { }
   function meta:__gc()
      cow.cow_domain_del(self._c)
   end
   function meta:__tostring()
      return "<DistributedDomain: ["..table.concat(self:size(), ", ").."]>"
   end
   setmetatable(self, meta)
   return self
end

local function test_lowlevel()
   local domain = cow.cow_domain_new()
   local dfield = cow.cow_dfield_new()

   cow.cow_domain_setndim(domain, 3)
   cow.cow_domain_setsize(domain, 0, 10)
   cow.cow_domain_setsize(domain, 1, 10)
   cow.cow_domain_setsize(domain, 2, 10)
   cow.cow_domain_commit(domain)

   assert(cow.cow_domain_getndim(domain) == 3)

   cow.cow_dfield_addmember(dfield, "vx")
   cow.cow_dfield_addmember(dfield, "vy")
   cow.cow_dfield_addmember(dfield, "vz")
   cow.cow_dfield_setname(dfield, "velocity")
   cow.cow_dfield_setdomain(dfield, domain)
   cow.cow_dfield_commit(dfield)

   assert(cow.cow_dfield_getname(dfield) == "velocity")

   cow.cow_domain_setchunk(domain, 1)
   cow.cow_domain_setcollective(domain, 1)
   cow.cow_dfield_write(dfield, "outfile.h5")
   cow.cow_dfield_read(dfield, "outfile.h5")
   cow.cow_dfield_del(dfield)
   cow.cow_domain_del(domain)
end

local function test_highlevel()
   local D = DistributedDomain{size={10,10}, guard=2}
   print(D)
end

test_lowlevel()
test_highlevel()

cow.cow_finalize()

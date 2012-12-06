


LUA_A = lua-5.2.1/lib/liblua.a
COW_A = ../lib/libcow.a
LUA_I = -Ilua-5.2.1/include
LIB = $(HDF5_LIB) $(FFTW_LIB)

default : ../bin/cow

cow_wrap.c : ../include/cow.h
	swig -module cow -o cow_wrap.c -lua ../include/cow.h
	patch -i lua52-swig.patch $@

../bin/cow : cowlua.c cow_wrap.c
	$(CC) -Wall -o $@ $^ $(LUA_I) $(LUA_A) $(COW_A) -I../include $(LIB)

clean :
	rm -f *.o cow cow_wrap.c cow_wrap.c.orig

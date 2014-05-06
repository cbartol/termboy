include nall/Makefile

ananke := ananke
gb  := gb

profile := accuracy
target  := ethos

# options += debugger
# arch := x86
# console := true

# compiler
flags   += -I. -O3 -fomit-frame-pointer
#flags   += -I. -ggdb -O0
link    += -lncursesw
link += -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-s
objects := libco

# profile-guided optimization mode
# pgo := instrument
# pgo := optimize

ifeq ($(pgo),instrument)
  flags += -fprofile-generate
  link += -lgcov
else ifeq ($(pgo),optimize)
  flags += -fprofile-use
endif

# platform
ifeq ($(platform),x)
  flags += -march=native
  link += #-Wl,-export-dynamic -ldl -lX11 -lXext
else ifeq ($(platform),osx)
  flags += -march=native
else ifeq ($(platform),win)
  ifeq ($(arch),win32)
    flags += -m32
    link += -m32
  endif
  ifeq ($(console),true)
    link += -mconsole
  else
    link += -mwindows
  endif
  link += -mthreads -luuid -lkernel32 -luser32 -lgdi32 -lcomctl32 -lcomdlg32 -lshell32 -lole32 -lws2_32
  link += -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc
else
  $(error unsupported platform.)
endif

ui := target-$(target)

# implicit rules
compile = \
    $(if $(filter %.c,$<), \
      $(compiler) $(cflags) $(flags) $1 -c $< -o $@, \
      $(if $(filter %.cpp,$<), \
        $(compiler) $(cppflags) $(flags) $1 -c $< -o $@ \
      ) \
    ) \

%.o: $<; $(call compile)

all: build;	

obj/libco.o: libco/libco.c libco/*

include $(ui)/Makefile
flags := $(flags) $(foreach o,$(call strupper,$(options)),-D$o)

# targets
clean:
	-@$(call delete,obj/*.o)
	-@$(call delete,obj/*.a)
	-@$(call delete,obj/*.so)
	-@$(call delete,obj/*.dylib)
	-@$(call delete,obj/*.dll)
	-@$(call delete,*.res)
	-@$(call delete,*.manifest)
	-@$(call delete,/tmp/benchmarks.*)

archive:
	if [ -f termboy.tar.xz ]; then rm termboy.tar.xz; fi
	tar -cJf termboy.tar.xz `ls`
sync:
ifeq ($(shell id -un),byuu)
	if [ -d ./libco ]; then rm -r ./libco; fi
	if [ -d ./nall ]; then rm -r ./nall; fi
	if [ -d ./ruby ]; then rm -r ./ruby; fi
	if [ -d ./phoenix ]; then rm -r ./phoenix; fi
	cp -r ../libco ./libco
	cp -r ../nall ./nall
	cp -r ../ruby ./ruby
	cp -r ../phoenix ./phoenix
	rm -r libco/doc
	rm -r libco/.test
	rm -r nall/.test
	rm -r ruby/.test
	rm -r phoenix/.test
endif

help:;

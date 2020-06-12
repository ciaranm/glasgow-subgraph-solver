TARGET := sip_to_lad

SOURCES := \
    sip_to_lad.cc

TGT_PREREQS := libcommon.a
ifeq ($(shell uname -s), Linux)
TGT_LDLIBS := libcommon.a $(boost_ldlibs) -lstdc++fs
else
TGT_LDLIBS := libcommon.a $(boost_ldlibs)
endif


TARGET := sip_to_opb

SOURCES := \
    sip_to_opb.cc

TGT_PREREQS := libcommon.a
ifeq ($(shell uname -s), Linux)
TGT_LDLIBS := libcommon.a $(boost_ldlibs) -lstdc++fs
else
TGT_LDLIBS := libcommon.a $(boost_ldlibs)
endif


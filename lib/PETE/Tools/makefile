
PROJECT_ROOT = ../../../
include $(PROJECT_ROOT)/config/head.mk

#CXX=g++-3.4
#CXXINC=-I..
#CXXARG=-O2 -g

default:: MakeOperators

MakeOperators: MakeOperators.cpp ClassDescriptor.h Header.h \
	OperatorDescriptor.h Parser.h PrintFunctions.h PrintOperators.h \
	PrintList.h Options.cpp Options.h PeteOps.cpp PeteOps.h \
	PAssert.cpp PAssert.h DescriptorBase.h Join.h Join.cpp
	$(CXX) -o MakeOperators -I.. $(CXX_OPT_APP_ARGS) \
	MakeOperators.cpp Options.cpp PAssert.cpp PeteOps.cpp Join.cpp

newpeteops: MakeOperators
	MakeOperators --operators PeteOps.in --guard "" --insert-op \
	--lanl-boilerplate --o PeteOps.cpp

peteheaders: OperatorTags.h

OperatorTags.h: MakeOperators
	MakeOperators --guard PETE_PETE_OPERATORTAGS_H --op-tags \
	--operators PeteOps.in --lanl-boilerplate --o OperatorTags.h

installheaders: OperatorTags.h
	cp OperatorTags.h ../PETE

clean:
	rm -f MakeOperators
	rm -f OperatorTags.h
	rm -f *.o
	rm -rf ti_files
	rm -rf ii_files

# ACL:rcsinfo
#  ----------------------------------------------------------------------
#  $RCSfile: makefile,v $   $Author: richard $
#  $Revision: 1.2 $   $Date: 2004/11/01 18:16:57 $
#  ----------------------------------------------------------------------
# ACL:rcsinfo

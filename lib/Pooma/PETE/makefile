
MAKEOPERATORS = ../../PETE/Tools/MakeOperators

#
# the files we will need to build
#

COMMONFILES = \
	PoomaOperatorTags.h

ARRAYFILES = \
	ArrayOperators.h \
	PoomaArrayOperators.h \
	VectorArrayOperators.h

DYNAMICARRAYFILES = \
	DynamicArrayOperators.h \
	PoomaDynamicArrayOperators.h \
	VectorDynamicArrayOperators.h

FIELDFILES = \
	FieldOperators.h \
	PoomaFieldOperators.h \
	VectorFieldOperators.h

ALLFILES = $(COMMONFILES) $(ARRAYFILES) $(DYNAMICARRAYFILES) $(FIELDFILES)


#
# the standard build targets
#

default: $(ALLFILES)

clean:
	rm -f $(ALLFILES)

install: $(ALLFILES)
	cp -f PoomaOperatorTags.h ..
	cp -f ArrayOperators.h ../../Array
	cp -f PoomaArrayOperators.h ../../Array
	cp -f VectorArrayOperators.h ../../Array
	cp -f PoomaDynamicArrayOperators.h ../../DynamicArray
	cp -f VectorDynamicArrayOperators.h ../../DynamicArray
	cp -f DynamicArrayOperators.h ../../DynamicArray
	cp -f FieldOperators.h ../../Field
	cp -f PoomaFieldOperators.h ../../Field
	cp -f VectorFieldOperators.h ../../Field

diff: $(ALLFILES)
	-diff -u ../PoomaOperatorTags.h PoomaOperatorTags.h
	-diff -u ../../Array/ArrayOperators.h ArrayOperators.h
	-diff -u ../../Array/PoomaArrayOperators.h PoomaArrayOperators.h
	-diff -u ../../Array/VectorArrayOperators.h VectorArrayOperators.h
	-diff -u ../../DynamicArray/PoomaDynamicArrayOperators.h PoomaDynamicArrayOperators.h
	-diff -u ../../DynamicArray/VectorDynamicArrayOperators.h VectorDynamicArrayOperators.h
	-diff -u ../../DynamicArray/DynamicArrayOperators.h DynamicArrayOperators.h
	-diff -u ../../Field/FieldOperators.h FieldOperators.h
	-diff -u ../../Field/PoomaFieldOperators.h PoomaFieldOperators.h
	-diff -u ../../Field/VectorFieldOperators.h VectorFieldOperators.h

#
# Common operations
#

PoomaOperatorTags.h:
	$(MAKEOPERATORS) --operators PoomaOps.in \
	--guard POOMA_POOMA_POOMAOPERATORTAGS_H \
	--no-expression --op-tags \
	--o PoomaOperatorTags.h

#
# Array-related operations
#

ArrayOperators.h:
	$(MAKEOPERATORS) --classes PoomaClass.in \
	--guard POOMA_ARRAY_ARRAYOPERATORS_H \
	--no-expression \
	--o ArrayOperators.h

PoomaArrayOperators.h:
	$(MAKEOPERATORS) --classes PoomaClass.in --operators PoomaOps.in \
	--guard POOMA_POOMA_POOMAARRAYOPERATORS_H \
	--no-expression \
	--o PoomaArrayOperators.h

VectorArrayOperators.h:
	$(MAKEOPERATORS) --classes PoomaVectorClass.in \
	--guard POOMA_POOMA_VECTORARRAYOPERATORS_H \
	--operators PoomaVectorOps.in \
	--no-expression --scalars \
	--o VectorArrayOperators.h

#
# DynamicArray-related operations
#

DynamicArrayOperators.h:
	$(MAKEOPERATORS) --classes PoomaDynamicArray.in \
	--guard POOMA_DYNAMIC_ARRAY_DYNAMICARRAYOPERATORS_H \
	--no-expression --extra-classes \
	--o DynamicArrayOperators.h

PoomaDynamicArrayOperators.h:
	$(MAKEOPERATORS) --classes PoomaDynamicArray.in \
	--guard POOMA_POOMA_POOMADYNAMICARRAYOPERATORS_H \
	--operators PoomaOps.in \
	--no-expression --extra-classes \
	--o PoomaDynamicArrayOperators.h

VectorDynamicArrayOperators.h:
	$(MAKEOPERATORS) --classes PoomaDynamicArrayVector.in \
	--guard POOMA_POOMA_VECTORDYNAMICARRAYOPERATORS_H \
	--operators PoomaVectorOps.in \
	--no-expression --scalars \
	--o VectorDynamicArrayOperators.h

#
# Field-related operations
#

FieldOperators.h:
	$(MAKEOPERATORS) --classes PoomaField.in \
	--guard POOMA_FIELD_FIELDOPERATORS_H \
	--no-expression \
	--o FieldOperators.h

PoomaFieldOperators.h:
	$(MAKEOPERATORS) --classes PoomaField.in \
	--guard POOMA_FIELD_POOMAFIELDOPERATORS_H \
	--operators PoomaOps.in \
	--no-expression \
	--o PoomaFieldOperators.h

VectorFieldOperators.h:
	$(MAKEOPERATORS) --classes PoomaFieldVector.in \
	--guard POOMA_POOMA_VECTORFIELDOPERATORS_H \
	--operators PoomaVectorOps.in \
	--no-expression --scalars \
	--o VectorFieldOperators.h

# ACL:rcsinfo
#  ----------------------------------------------------------------------
#  $RCSfile: makefile,v $   $Author: richard $
#  $Revision: 1.10 $   $Date: 2004/11/01 18:17:05 $
#  ----------------------------------------------------------------------
# ACL:rcsinfo

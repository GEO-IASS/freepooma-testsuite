set(DIFF  "PS_cmp.py -T @ -u ?")
qa_set_run_time_limit(150)
qa_set_compile_time_limit(400)
qa_add_compiler_flags(-I${test_suite_path}/lib)

set(lib_sources
    Connect/Lux/LuxAppPointer.cmpl.C
    Connect/Paws/PawsAppPointer.cmpl.C
    Connect/Connection.cmpl.C
    DataBrowser/DataBrowser.cmpl.C
    Domain/DomainCalculus.cmpl.C
    IO/DiskLayout.cmpl.C
    IO/DiskMeta.cmpl.C
    IO/FileSetReader.cmpl.C
    IO/MetaTokenIterator.cmpl.C
    Layout/DynamicLayout.cmpl.C
    Layout/GlobalIDDataBase.cmpl.C
    Field/Relations/RelationGroups.cmpl.C
    Field/FieldCentering.cmpl.C
    Particles/AttributeList.cmpl.C
    Particles/ParticleBCList.cmpl.C
    Partition/UniformMapper.cmpl.C
    Pooma/Pooma.cmpl.C
    Threads/IterateSchedulers/SerialAsync.cmpl.C
    Tulip/Messaging.cmpl.C
    Tulip/PatchSizeSyncer.cmpl.C
    Utilities/Benchmark.cmpl.C
    Utilities/Inform.cmpl.C
    Utilities/Options.cmpl.C
    Utilities/PAssert.cmpl.C
    Utilities/Pool.cmpl.C
    Utilities/Statistics.cmpl.C
    Utilities/Tester.cmpl.C
    Utilities/Unique.cmpl.C
)


set(real_lib_sources )
foreach(src ${lib_sources})
    set(real_lib_sources ${real_lib_sources} ${test_suite_path}/lib/${src})
endforeach()

qa_make_library(${qa_work_path}/libpooma-gcc.a ${real_lib_sources})


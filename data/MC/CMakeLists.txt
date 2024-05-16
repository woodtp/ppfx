file(GLOB root_files *.root)
install( FILES ${root_files} DESTINATION ${product}/${version}/data/MC )

add_subdirectory(QGSPv9.2)
add_subdirectory(FTFP)
add_subdirectory(FTFPv9.2)

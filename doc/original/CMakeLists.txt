file(GLOB docs *.pdf *.docx *.txt)
install(FILES ${docs} DESTINATION doc)
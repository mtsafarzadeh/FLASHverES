
IO += IO_init.o \
	IO_output.o \
	IO_setScalar.o \
	IO_readCheckpoint.o \
	IO_getScalar.o \
	IO_outputInitial.o \
	IO_outputFinal.o \
	IO_sendOutputData.o \
	IO_getPrevScalar.o \
	IO_writeIntegralQuantities.o \
	IO_writeParticles.o \
	IO_readParticles.o \
	IO_writeCheckpoint.o \
	IO_updateScalars.o \
	IO_writePlotfile.o \
	IO_finalize.o \
	IO_writeUserArray.o \
	IO_readUserArray.o \
	IO_initRPsFromCheckpoint.o \
        IO_interface.o \
	IO_checkForPlot.o \
        IO_startRayWrite.o \
	IO_writeRays.o \
        IO_endRayWrite.o

IO += IO_data.o \
	io_getAllScalars.o io_getNumScalars.o \
	io_prepareListsWrite.o io_prepareListsRead.o \
	io_finalizeListsRead.o io_finalizeListsWrite.o \
	io_bcastScalars.o io_getOutputName.o \
	io_prepareSimInfo.o \
	io_isPlotVar.o \
	io_setPrevScalar.o \
	io_getVarExtrema.o \
	io_memoryReport.o \
	io_restrictBeforeWrite.o \
	io_mpi_type.o \
	io_xfer_cont_slab.o \
	io_xfer_tree_data.o \
	io_attribute.o \
	io_create_dataset.o \
	io_c_interface.o \
	io_repack_data.o

io_flash.h : Flash.h; test -f $@ && ${TOUCH} $@

IO_output.o : IO_data.o Grid_sendOutputData.o
IO_init.o : RuntimeParameters_get.o IO_data.o
IO_data.o : nameValueLL_data.o
IO_setScalar.o : IO_data.o
IO_getScalar.o : IO_data.o
io_getNumScalars.o : IO_data.o
io_getAllScalars.o : IO_data.o
io_getVarExtrema.o : IO_data.o
io_bcastScalars.o : IO_data.o
IO_writePlotfile.o : io_getOutputName.o
IO_writeCheckpoint.o : Logfile_stamp.o
IO += io_h5write_header.o\
	io_h5write_lists.o \
	io_h5read_header.o \
	io_h5read_lists.o \
	io_h5read_unknowns.o \
	io_h5write_unknowns.o \
	io_h5write_unknowns_sp.o \
	io_h5write_generic_int_arr.o \
	io_h5write_generic_real_arr.o \
	io_h5read_generic_int_arr.o \
	io_h5read_generic_real_arr.o \
	io_comp_decomp.o \
	io_h5_xfer_wrapper.o \
	io_h5_xfer.o \
	io_h5create_dataset.o \
	io_h5_type.o \
	io_h5_attribute.o \
	io_h5_type_matched_xfer.o \
	io_h5create_raydset.o \
	io_h5write_raydata.o \
	io_h5_read_present_dims.o \
	io_h5_report_xfer_method.o \
	io_createDatasets.o

io_h5_attribute.o: constants.h Flash.h io_flash.h io_h5_type.h
#	Makefile for FLASH HDF5 output (serial version)

IO += io_h5file_interface.o \
	io_initFile.o io_closeFile.o




io_h5file_interface.o : hdf5_flash.h constants.h
IO += io_writeData.o \
	io_readData.o



IO += io_ptCorrectNextPartTime.o \
      io_ptResetNextFile.o \
      io_ptSendOutputData.o \
      io_ptInit.o \
      io_readRPsFromCheckpoint.o \
      io_ptInterface.o \
      io_typeInterface.o \
      io_rescaleCellBoxes.o


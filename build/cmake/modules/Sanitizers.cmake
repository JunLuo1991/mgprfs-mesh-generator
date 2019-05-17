# Define the ENABLE_ASAN and ENABLE_UBSAN options.
option(ENABLE_ASAN "Enable the Address Sanitizer (ASAN)" false)
option(ENABLE_UBSAN "Enable the Undefined Behavior Sanitizer (UBSAN)" false)

# Enable ASAN, if appropriate.
if (ENABLE_ASAN)
	if (CMAKE_CXX_COMPILER_ID MATCHES GNU OR
	  CMAKE_CXX_COMPILER_ID MATCHES Clang)
		# Handle the cases of GCC and Clang.
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
		set(CMAKE_EXE_LINKER_FLAGS
		  "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
	else()
		# Handle the case of unsupported compilers.
		message(FATAL_ERROR "ASAN support is unavailable.")
	endif()
endif()

# Enable UBSAN, if appropriate.
if (ENABLE_UBSAN)
	if (CMAKE_CXX_COMPILER_ID MATCHES GNU OR
	  CMAKE_CXX_COMPILER_ID MATCHES Clang)
		# Handle the cases of GCC and Clang.
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")
		set(CMAKE_EXE_LINKER_FLAGS
		  "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=undefined")
	else()
		# Handle the case of unsupported compilers.
		message(FATAL_ERROR "UBSAN support is unavailable.")
	endif()
endif()

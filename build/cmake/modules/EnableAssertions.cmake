# Define the ALWAYS_ENABLE_ASSERTIONS option.
option(ALWAYS_ENABLE_ASSERTIONS
  "Always enable assertions (even for release builds)" true)

# Delete any compiler option that defines the preprocessor symbol NDEBUG.
if (ALWAYS_ENABLE_ASSERTIONS)
	foreach (var
		CMAKE_CXX_FLAGS_DEBUG
		CMAKE_CXX_FLAGS_RELEASE
		CMAKE_CXX_FLAGS_RELWITHDEBINFO
		CMAKE_CXX_FLAGS_MINSIZEREL)
		# Print the original value of the variable.
		message(STATUS "Original value of ${var}: ${${var}}")
		if (CMAKE_CXX_COMPILER_ID MATCHES GNU OR
		  CMAKE_CXX_COMPILER_ID MATCHES Clang)
			# Handle the case of GCC and Clang
			# Remove -DNDEBUG or -D NDEBUG.
			string(REGEX REPLACE "(^| )-D *NDEBUG($| )" " " ${var} "${${var}}")
		elseif (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
			# Handle the case of MSVC.
			# Remove /DNDEBUG or /D NDEBUG.
			string(REGEX REPLACE "(^| )/D *NDEBUG($| )" " " ${var} "${${var}}")
		endif()
		# Print the new value of the variable.
		message(STATUS "New value of ${var}: ${${var}}")
	endforeach()
endif()

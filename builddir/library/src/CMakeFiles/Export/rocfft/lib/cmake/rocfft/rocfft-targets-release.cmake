#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "roc::rocfft-device-misc" for configuration "Release"
set_property(TARGET roc::rocfft-device-misc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(roc::rocfft-device-misc PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/rocfft/lib/librocfft-device-misc.so.0.1"
  IMPORTED_SONAME_RELEASE "librocfft-device-misc.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS roc::rocfft-device-misc )
list(APPEND _IMPORT_CHECK_FILES_FOR_roc::rocfft-device-misc "${_IMPORT_PREFIX}/rocfft/lib/librocfft-device-misc.so.0.1" )

# Import target "roc::rocfft-device-single" for configuration "Release"
set_property(TARGET roc::rocfft-device-single APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(roc::rocfft-device-single PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/rocfft/lib/librocfft-device-single.so.0.1"
  IMPORTED_SONAME_RELEASE "librocfft-device-single.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS roc::rocfft-device-single )
list(APPEND _IMPORT_CHECK_FILES_FOR_roc::rocfft-device-single "${_IMPORT_PREFIX}/rocfft/lib/librocfft-device-single.so.0.1" )

# Import target "roc::rocfft-device-double" for configuration "Release"
set_property(TARGET roc::rocfft-device-double APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(roc::rocfft-device-double PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/rocfft/lib/librocfft-device-double.so.0.1"
  IMPORTED_SONAME_RELEASE "librocfft-device-double.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS roc::rocfft-device-double )
list(APPEND _IMPORT_CHECK_FILES_FOR_roc::rocfft-device-double "${_IMPORT_PREFIX}/rocfft/lib/librocfft-device-double.so.0.1" )

# Import target "roc::rocfft" for configuration "Release"
set_property(TARGET roc::rocfft APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(roc::rocfft PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "roc::rocfft-device-single;roc::rocfft-device-double;roc::rocfft-device-misc"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/rocfft/lib/librocfft.so.0.1"
  IMPORTED_SONAME_RELEASE "librocfft.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS roc::rocfft )
list(APPEND _IMPORT_CHECK_FILES_FOR_roc::rocfft "${_IMPORT_PREFIX}/rocfft/lib/librocfft.so.0.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

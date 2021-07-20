# Install script for directory: /home/ugrad/ugrad003/roctosycl/rocFFT/library/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/opt/rocm")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/rocfft/include" TYPE DIRECTORY FILES "/home/ugrad/ugrad003/roctosycl/rocFFT/library/include/" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$" REGEX "/[^/]*\\.hh$" REGEX "/[^/]*\\.hxx$" REGEX "/[^/]*\\.inl$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/rocfft/include" TYPE DIRECTORY FILES "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/include/" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$" REGEX "/[^/]*\\.hh$" REGEX "/[^/]*\\.hxx$" REGEX "/[^/]*\\.inl$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so.0"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/rocfft/lib" TYPE SHARED_LIBRARY FILES
    "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/librocfft.so.0.1"
    "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/librocfft.so.0"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so.0"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/device:"
           NEW_RPATH "")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/rocfft/lib" TYPE SHARED_LIBRARY FILES "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/librocfft.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so"
         OLD_RPATH "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/device:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/librocfft.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/cmake/rocfft/rocfft-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/cmake/rocfft/rocfft-targets.cmake"
         "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/CMakeFiles/Export/rocfft/lib/cmake/rocfft/rocfft-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/cmake/rocfft/rocfft-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft/lib/cmake/rocfft/rocfft-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/rocfft/lib/cmake/rocfft" TYPE FILE FILES "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/CMakeFiles/Export/rocfft/lib/cmake/rocfft/rocfft-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/rocfft/lib/cmake/rocfft" TYPE FILE FILES "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/CMakeFiles/Export/rocfft/lib/cmake/rocfft/rocfft-targets-release.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/rocfft/lib/cmake/rocfft" TYPE FILE FILES
    "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/rocfft-config.cmake"
    "/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/rocfft-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  
        set(SUBDIR $ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/rocfft)
        file(GLOB_RECURSE FILES RELATIVE ${SUBDIR} ${SUBDIR}/*)
        foreach(FILE ${FILES})
            set(SRC ${SUBDIR}/${FILE})
            set(DEST $ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${FILE})
            get_filename_component(DEST_DIR ${DEST} DIRECTORY)
            file(MAKE_DIRECTORY ${DEST_DIR})
            file(RELATIVE_PATH SRC_REL ${DEST_DIR} ${SRC})
            message(STATUS "symlink: ${SRC_REL} -> ${DEST}")
            execute_process(COMMAND ln -sf ${SRC_REL} ${DEST})
        endforeach()
    
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/ugrad/ugrad003/roctosycl/rocFFT/builddir/library/src/device/cmake_install.cmake")

endif()


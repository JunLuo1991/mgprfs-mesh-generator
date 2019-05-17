file(GLOB USER_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/../util.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/../image.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/../scan_triangle.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/../error_diffusion.cpp"
)

find_package(OpenCV REQUIRED)

include_directories(
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
	${SPL_INCLUDE_DIRS}
	${SPLEL_INCLUDE_DIRS}
	${CGAL_INCLUDE_DIR}
)

link_libraries(
        ${SPL_LIBRARIES}
        ${SPLEL_LIBRARIES}
        ${CGAL_LIBRARY}
        ${GMP_LIBRARIES}
        ${Boost_LIBRARIES}
)

# TIP: UNCOMMENT FOLLOWING LINES TO COMPILE image_compare PROGRAM
#add_executable(image_compare image_compare.cpp)
#target_include_directories(image_compare PUBLIC ${OpenCV_INCLUDE_DIRS})
#target_link_libraries(image_compare LINK_PUBLIC ${OpenCV_LIBS})


add_executable(output_mesh_size output_mesh_size.cpp)

add_executable(output_image_attributes output_image_attributes.cpp
               ../image.cpp ../util.cpp)

add_executable(data_compare data_compare.cpp)

add_executable(off_to_triangulation off_to_triangulation.cpp)

add_executable(mirror_image mirror_image.cpp
               ../image.cpp ../util.cpp)

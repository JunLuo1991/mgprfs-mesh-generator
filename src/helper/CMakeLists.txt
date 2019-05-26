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


add_executable(output_mesh_size output_mesh_size.cpp)

add_executable(output_image_attributes output_image_attributes.cpp
               ../image.cpp ../util.cpp)

add_executable(data_compare data_compare.cpp)

add_executable(off_to_triangulation off_to_triangulation.cpp)

add_executable(mirror_image mirror_image.cpp
               ../image.cpp ../util.cpp)

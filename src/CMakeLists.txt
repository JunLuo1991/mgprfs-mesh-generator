include_directories(
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

add_executable(make_model make_model.cpp mesh_generator.cpp
               image.cpp scan_triangle.cpp error_diffusion.cpp
               util.cpp pixel.hpp)

add_executable(rasterize_model rasterize_model.cpp scan_triangle.cpp 
               image.cpp util.cpp pixel.hpp)


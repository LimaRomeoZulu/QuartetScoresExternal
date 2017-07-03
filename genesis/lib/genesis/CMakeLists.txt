file (GLOB_RECURSE genesis_lib_sources   ${PROJECT_SOURCE_DIR}/lib/genesis/*.cpp)
include_directories (${PROJECT_SOURCE_DIR}/lib/genesis)
add_library           (genesis_lib_static SHARED ${genesis_lib_sources})
set_target_properties (genesis_lib_static PROPERTIES OUTPUT_NAME genesis)
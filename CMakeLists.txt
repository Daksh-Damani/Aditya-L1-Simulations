cmake_minimum_required(VERSION 3.25)

project(AdityaL1
    VERSION 1.0.0
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# GLFW
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

# GLM
FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG 1.0.1
)
FetchContent_MakeAvailable(glm)

# Dear ImGui 1.92.9 (current release as of this build guide)
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.92.9
)
FetchContent_MakeAvailable(imgui)

add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)

target_include_directories(imgui PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui PUBLIC glfw)

# Local GLAD 2 generated from gen.glad.sh
add_library(glad STATIC
    external/glad/src/gl.c
)
target_include_directories(glad PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include
)

add_executable(AdityaL1
    src/main.cpp
)

target_link_libraries(AdityaL1 PRIVATE
    glfw
    glm::glm
    glad
    imgui
)

target_include_directories(AdityaL1 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)

if(MSVC)
    target_compile_options(AdityaL1 PRIVATE /W4 /permissive-)
    target_compile_options(imgui PRIVATE /W2)
endif()

if(WIN32)
    target_link_libraries(AdityaL1 PRIVATE opengl32)
elseif(APPLE)
    find_library(OPENGL_FRAMEWORK OpenGL)
    target_link_libraries(AdityaL1 PRIVATE ${OPENGL_FRAMEWORK})
elseif(UNIX)
    find_package(OpenGL REQUIRED)
    target_link_libraries(AdityaL1 PRIVATE OpenGL::GL)
endif()

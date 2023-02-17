cmake_minimum_required(VERSION 3.14)
project (3dy4_test)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -O3)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
FetchContent_MakeAvailable(googletest)

include_directories( ./include)
file(GLOB SOURCES "src/*.c*" "test/*.c*")

set(PROJECT_LINK_LIBS pthread)
set(PROJECT_SOURCES ${SOURCES})
list(FILTER PROJECT_SOURCES EXCLUDE REGEX ".*[Tt]est\\.c.*$")
add_executable(project ${PROJECT_SOURCES})
target_link_libraries(project ${PROJECT_LINK_LIBS})

set(TEST_LINK_LIBS gtest_main)
set(TEST_SOURCES ${SOURCES})
list(FILTER TEST_SOURCES EXCLUDE REGEX ".*[Pp]roject\\.c.*$")
add_executable(test_fourier ${TEST_SOURCES})
target_link_libraries(test_fourier ${PROJECT_LINK_LIBS} ${TEST_LINK_LIBS})
include(GoogleTest)
gtest_discover_tests(test_fourier)

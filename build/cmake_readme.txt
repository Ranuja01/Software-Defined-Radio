
The use of "cmake" is optional, however it is highly recommended if you wish to unit test your new features (or the refactored code) incrementally. Follow the following few simple steps:

1) Keep your C++ source files in "src" sub-folder; your header files in the "include" sub-folder and your unit tests in the "test" sub-folder. It is assumed the top-level file of the project (with the "main" function) is in a file called "project.cpp" and all (and only) the files from the "test" sub-folder include "test" in the name.

2) Change to the "build" sub-folder and type in:

cmake ..

3) A Makefile is produced using the default configuration from CMakeLists.txt from the root directory of the project. After the Makefile has been produced, when you type in "make" you will generate two separate executables: one for the project and one for the unit tests.

For more details about Google's test framework check:

https://github.com/google/googletest/blob/main/docs/index.md

4) To generate a debug version of your executable, from the root directory of your project folder do the following:

mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

The executables from the "debug" sub-folder are now annotated with debug information.

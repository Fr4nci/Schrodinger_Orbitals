cmake_minimum_required(VERSION 3.10)
project(SchrodingerAtomo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova GLFW nel sistema
find_package(glfw3 REQUIRED)

# Aggiungi GLAD
add_library(glad external/glad/src/glad.c)
target_include_directories(glad PUBLIC external/glad/include)

# Eseguibile principale
add_executable(orbitals src/main.cpp)

# Linka le librerie
target_link_libraries(orbitals PRIVATE glad glfw)

# Copia la cartella shaders nella directory di build accanto all'eseguibile
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

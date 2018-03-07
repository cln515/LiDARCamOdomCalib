#pragma once

#include <Windows.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

std::vector<std::string> getFileWithNumber(std::string folder, std::string filename, std::string extension);

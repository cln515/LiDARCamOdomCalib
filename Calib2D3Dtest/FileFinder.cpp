#include "FileFinder.h"

std::vector<std::string> getFileWithNumber(std::string folder, std::string filename, std::string extension)
{
	HANDLE hFind;
	WIN32_FIND_DATAA win32fd;//defined at Windwos.h
	std::vector<std::string> filelist;

	int idx = 1;
	while (true) {
		std::stringstream ss;
		ss<< folder << "\\" << filename <<idx<< "." << extension;
		std::string search_name = ss.str();
		hFind = FindFirstFileA(search_name.c_str(), &win32fd);
		if (hFind == INVALID_HANDLE_VALUE) {
			break;
		}
		if (win32fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
		}
		else {
			filelist.push_back(search_name);
		}
		FindClose(hFind);
		idx++;
	}
	

	return filelist;
}
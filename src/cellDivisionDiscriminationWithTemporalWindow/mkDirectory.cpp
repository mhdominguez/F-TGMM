#include <iostream>
#include <string>

#if defined(_WIN32) || defined(_WIN64)
    #include <windows.h>
#else
	#include <unistd.h>
	#include "sys/stat.h"
#endif

using namespace std;

int mkDirectory(const std::string& folder)
{

#if defined(_WIN32) || defined(_WIN64)
		if (GetFileAttributes(folder.c_str()) == INVALID_FILE_ATTRIBUTES)//check if folder exists
		{

			bool errB = CreateDirectory(folder.c_str(), NULL);
			if (errB == false)
			{
				cout << "ERROR: mkDirectory" << folder << ": generating folder " << folder << endl;
				return 1;
			}
		}

#else
		struct stat St;
		if (stat(folder.c_str(), &St) != 0)//check if folder exists
		{
			string cmd("mkdir " + folder);
			int error = system(cmd.c_str());
			if (error>0)
			{
				cout << "ERROR (" << error << "): generating path " << folder << endl;
				cout << "Wtih command " << cmd << endl;
				return error;
			}
		}
#endif
	return 0;
}
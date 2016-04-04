#include <iostream>
#include <string>

#include "FD_Sn_Transport.hh"

using namespace std;

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        cerr << "usage: hare [file.xml]" << endl;
        return 1;
    }

    string filename = argv[1];
    
    FD_Sn_Transport transport(filename);
}

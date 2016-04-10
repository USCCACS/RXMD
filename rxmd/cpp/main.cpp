#include <iostream>
#include <fstream>
#include <string>

#include "params.h"

int main(int argc, char *argv[])
{
    for(int i=0; i<argc; ++i)
    {
        std::cout << i << " " << argv[i] << std::endl;
    }

    Params(std::string(argv[1]));

}

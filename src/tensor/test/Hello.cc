#include <iostream>
#include "Hello.h"

using namespace std;

Hello::Hello()
{
    name_ = "Tom";

}

Hello::Hello(std::string name)
{
    name_ = name;

}

void Hello::hello()
{
    cout << "Hello" << name_ << endl;
}

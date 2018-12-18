#ifndef HELLO_H
#define HELLO_H

using namespace std;

class Hello
{
    private:
    std::string name_;
    public:
    Hello();
    Hello(std::string name);
    void hello();
};

#endif
